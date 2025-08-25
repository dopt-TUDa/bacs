/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is part of the program BACS                                  */
/*                                                                           */
/*    an implementation of a Branch-And-Cut algorithm to solve the           */
/*    stable set problem.                                                    */
/*                                                                           */
/*    Copyright (C) 2024-  Discrete Optimization Group, TU Darmstadt         */
/*                                                                           */
/*                                                                           */
/*    Based on SCIP  --- Solving Constraint Integer Programs                 */
/*    Copyright (C) 2002-2025 Zuse Institute Berlin (ZIB)                    */
/*                                                                           */
/*    Both are licensed under the Apache License, Version 2.0.               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_neigh.cpp
 * @brief  Separator for neighborhood inequalities
 * @author Marc Pfetsch
 *
 * For each node \f$v\f$, we compute the stability number the neighborhood \f$S\f$. We add the corresponding
 * neighborhood inequality
 * \f[
 *   \sum_{i \in S} x_{ij} + \alpha(S) x_{vj} \leq \alpha(S),
 * \f]
 * for all \f$j\f$. We initially generate all neighborhood sets \f$S\f$ an check these for violated
 * cuts in each separation round.
 *
 * The inequalities here are investigated in the paper:@p
 * The Stable Set Problem: Clique and Nodal Inequalities Revisited@p
 * Adam N. Letchford, Fabrizio Rossi, Stefano Smriglio
 */

#include "sepa_neigh.h"
#include "probdata_bacs.h"
#include "vectorpool.hpp"
#include "prop_neighborhoods.h"
#include "heur_greedylp.h"
#include "heur_greedydeg.h"
#include "heur_greedyrounding.h"
#include "heur_tabu.h"
#include "heur_dynamicdeg.h"

#include <cassert>
#include <sys/times.h>

#define SEPA_NAME                 "neigh"
#define SEPA_DESC                 "neighborhood inequalities"
#define SEPA_PRIORITY                 -20
#define SEPA_FREQ                      -1
#define SEPA_MAXBOUNDDIST             0.0
#define SEPA_USESUBSCIP             FALSE
#define SEPA_DELAY                  FALSE


// default values:
#define DEFAULT_MAXNEIGHBORS          500
#define DEFAULT_MAXNEIGHFRACT         0.1
#define DEFAULT_NODELIMIT            5000
#define MAXLIMITREACHED                10


//! data for neighborhood separator
struct SCIP_SepaData
{
   vectorpool*           pool;               //!< pool
   int                   maxneighbors;       //!< maximal number of nodes in neighborhood to compute stable set for
   SCIP_Real             maxneighfract;      //!< maximal fraction of number of nodes present in neighborhood
   int                   nodelimit;          //!< nodelimit for branch-and-bound in subscip
};



//--------------------------------------------------------------------------------------------
//------------------------------------- local functions --------------------------------------
//--------------------------------------------------------------------------------------------

//! compute stability number for the subgraph induced by the neighbors of @a node
static
SCIP_RETCODE computeStabilityNumber(
   SCIP*                 scip,               //!< SCIP instance
   int                   nodelimit,          //!< nodelimit for branch-and-bound
   const Graph&          G,                  //!< original graph
   Vertex                node,               //!< current node
   size_t                vsize,              //!< number of nodes in the subgraph
   const std::vector<int>& ind,              //!< indices of nodes in V
   int&                  alpha               //!< resulting stability number, -1 if computation was not optimal
   )
{
   assert( scip != nullptr );

   // init computation as invalid
   alpha = -1;

   // create a new SCIP instance
   SCIP* subscip;
   SCIP_CALL( SCIPcreate(&subscip) );

   // copy plugins, we omit pricers (because we do not run if there are active pricers) and dialogs
   SCIP_Bool success;
#ifdef SCIP_MORE_DEBUG // we print statistics later, so we need to copy statistics tables
   SCIP_CALL( SCIPcopyPlugins(scip, subscip, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, &success) );
#else
   SCIP_CALL( SCIPcopyPlugins(scip, subscip, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
#endif

   if ( ! success )
   {
      SCIP_CALL( SCIPfree(&subscip) );
      return SCIP_OKAY;
   }

   SCIP_CALL( BACSincludePropNeighborhoods(subscip) );
   SCIP_CALL( BACSincludeHeurGreedyLP(subscip) );
   SCIP_CALL( BACSincludeHeurGreedyDeg(subscip) );
   SCIP_CALL( BACSincludeHeurGreedyRounding(subscip) );
   SCIP_CALL( BACSincludeHeurTabu(subscip) );
   SCIP_CALL( BACSincludeHeurDynamicdeg(subscip) );

   // copy parameter settings
   SCIP_CALL( SCIPcopyParamSettings(scip, subscip) );

   // disable presolving
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   // disable output, unless in extended debug mode
#ifndef SCIP_MORE_DEBUG
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
#else
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
#endif

   // create data of subproblem
   SCIP_PROBDATA* subprobdata;
   SCIP_CALL( BACSinitProblem(subscip, &subprobdata) );

   // construct graph
   Graph* SG = new Graph();

   Vertex* vertexmap;
   size_t n = boost::num_vertices(G);
   SCIP_CALL( SCIPallocBufferArray(scip, &vertexmap, n) );

   // create nodes
   subprobdata->unweighted = TRUE;
   for (size_t i = 0; i < n; ++i)
   {
      if ( ind[i] >= 0 )
      {
         Vertex v = boost::add_vertex(*SG);

         // treat weights
         SCIP_Real weight = boost::get(vertex_weight_t(), G, i);
         if ( ! SCIPisEQ(scip, weight, 1.0) )
            subprobdata->unweighted = FALSE;
         boost::put(vertex_weight_t(), *SG, v, weight);

         vertexmap[i] = v;
         assert( v <= vsize );
      }
      else
         vertexmap[i] = Graph::null_vertex();
   }

   // create edges
   EdgeIterator eit, eend;
   for (boost::tie(eit, eend) = boost::edges(G); eit != eend; ++eit)
   {
      Vertex s = boost::source(*eit, G);
      Vertex t = boost::target(*eit, G);

      if ( ind[s] < 0 || ind[t] < 0 )
         continue;

      Vertex subs = vertexmap[s];
      Vertex subt = vertexmap[t];
      assert( subs != Graph::null_vertex() && subt != Graph::null_vertex() );

      // create edge
      (void) boost::add_edge(subs, subt, *SG);
   }
   SCIPfreeBufferArray(scip, &vertexmap);

   subprobdata->G = SG;
   subprobdata->n = boost::num_vertices(*SG);
   subprobdata->m = boost::num_edges(*SG);

   SCIP_CALL( BACSinitDatastructures(subscip, subprobdata) );
   char name[SCIP_MAXSTRLEN];
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "neigh%u", (unsigned int) node);
   SCIP_CALL( BACSsetupProblem(subscip, name, subprobdata) );

   // update subscip depth
   SCIPsetSubscipDepth(subscip, SCIPgetSubscipDepth(scip) + 1);

   /* SCIPcopyLimits will set wrong time limits since it does not take into account time spent already in the sub-SCIP;
    * nevertheless, we call it to set the memory limit and unset all other limits, if set in the main SCIP. */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );

   // set time and memory limit for the subproblem
   //SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );

   // set node limit
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );

   // solve the subproblem
   SCIP_CALL( SCIPsolve(subscip) );

   if ( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
   {
      alpha = SCIPconvertRealToInt(scip, SCIPgetPrimalbound(subscip));
   }
   else
   {
      SCIPwarningMessage(scip, "Stable set computation not optimal.\n");
   }

   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}


/** Generate neighborhood inequalities
 *
 *  For each node \f$v\f$, we compute the stability number the neighborhood \f$S\f$. We add the corresponding
 *  information for the neighborhood inequality
 *  \f[
 *    \sum_{i \in S} x_i + \alpha(S) x_v \leq \alpha(S),
 *  \f]
 *  into the pool.
 */
static
SCIP_RETCODE generateNeighborhoodInequalities(
   SCIP*                 scip,               //!< SCIP pointer
   int                   maxneighbors,       //!< maximal number of nodes in neighborhood to compute stable set for
   SCIP_Real             maxneighfract,      //!< maximal fraction of number of nodes present in neighborhood
   vectorpool*           pool,               //!< pool
   int                   nodelimit,          //!< nodelimit in subscip
   unsigned int&         ngen,               //!< number of generated inequalities
   unsigned int&         ntested             //!< number of nodes tested
   )
{
   ngen = 0;
   ntested = 0;

   assert( scip != nullptr );
   assert( pool != nullptr );

   // get problem data
   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   // get sizes and variables
   const Graph* G = probdata->G;
   assert( G != nullptr );
   size_t n = probdata->n;
   assert( n == boost::num_vertices(*G) );

   // loop through all nodes
   unsigned int nlimitreached = 0;
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;

      // avoid sets that are too small or too large
      int degree = probdata->localdegrees[v];
      if ( degree < 4 || degree > maxneighbors || degree > (int) (maxneighfract * n))
         continue;

      assert( SCIPvarGetUbLocal(probdata->vars[v]) > 0.5 && SCIPvarGetLbLocal(probdata->vars[v]) < 0.5 );

      // collect neighbors
      std::vector<int> ind(n, -1);

      size_t vsize = 0;
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait!= aend; ++ait)
         ind[*ait] = (int) vsize++;

      assert( (int) vsize == degree );

      // compute stability number of subgraph induced by neighbors
      int alpha = 0;
      SCIPdebugMsg(scip, "Computing stability number for node %zu with neighborhood of size %zu ...\n", v, vsize);
      SCIP_CALL( computeStabilityNumber(scip, nodelimit, *G, v, vsize, ind, alpha) );
      ++ntested;

      // only generate set if computation was successfull and the complete set is not stable
      // (the cut is meaningless in this case)
      if ( alpha > 0 && alpha < (int) vsize )
      {
         std::vector<unsigned int> S;
         S.push_back((unsigned) alpha);

         for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait!= aend; ++ait)
         {
            S.push_back(1);                     // coefficient
            S.push_back((unsigned int) *ait);   // node
         }
         // add central node
         S.push_back((unsigned) alpha);
         S.push_back((unsigned int) v);

         assert( S.size() % 2 == 1 );
         (void) pool->insert(S);
         ++ngen;
      }
      else if ( alpha == -1 )
         ++nlimitreached;

      if ( nlimitreached >= MAXLIMITREACHED )
      {
         SCIPwarningMessage(scip, "Stopping computation, because %u runs reached the node limit.\n", nlimitreached);
         break;
      }
   }

   return SCIP_OKAY;
}


/** Separate neighborhood inequalities from the pool
 *
 *  We check for all sets in the pool whether they lead to violated inequalities.
 */
static
SCIP_RETCODE separatePools(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< solution to be separated */
   vectorpool*           pool,               /**< pool */
   unsigned int&         ngen,               /**< number of separated inequalities on exit */
   unsigned int          maxgen,             /**< maximal number of cuts to generate */
   SCIP_Bool&            cutoff              /**< whether a cutoff has been detected */
   )
{
   assert( scip != nullptr );
   ngen = 0;
   cutoff = FALSE;

   SCIPdebugMsg(scip, "Separate pool of separator <%s>...\n", SEPA_NAME);

   // get problem data
   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   // get sizes and variables
   size_t n = probdata->n;
   SCIP_VAR** vars = probdata->vars;

   SCIP_Real* tmpvals = nullptr;
   SCIP_VAR** tmpvars = nullptr;
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvals, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvars, n) );

   // loop through all sets
   vectorpool::const_iterator pit, pend = pool->end();
   for (pit = pool->begin(); pit != pend; ++pit)
   {
      std::vector<unsigned int>* S = *pit;
      std::vector<unsigned int>::const_iterator vit, vend = S->end();

      // the first entry gives the stable set number
      std::vector<unsigned int>::const_iterator beg = S->begin();
      unsigned alpha = *beg++;
      assert( S->size() % 2 == 1 );

      SCIP_Real sum = 0.0;
      for (vit = beg; vit != vend; ++vit)
      {
         SCIP_Real weight = (SCIP_Real) *vit++; /*lint !e850*/
         sum += weight * SCIPgetSolVal(scip, sol, vars[*vit]);
      }

      // if the weight of the left hand side is larger than the corresponding rhs
      if ( SCIPisEfficacious(scip, sum - (SCIP_Real) alpha) )
      {
         int nvars = 0;
         for (vit = beg; vit != vend; ++vit)
         {
            tmpvals[nvars] = (SCIP_Real) *vit++; /*lint !e850*/
            tmpvars[nvars++] = vars[*vit];
         }

         SCIP_ROW* row;
#ifndef NDEBUG
         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, "neigh", -SCIPinfinity(scip), (SCIP_Real) alpha, FALSE, FALSE, TRUE) );
#else
         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepa, "", -SCIPinfinity(scip), (SCIP_Real) alpha, FALSE, FALSE, TRUE) );
#endif
         SCIP_CALL( SCIPaddVarsToRow(scip, row, nvars, tmpvars, tmpvals) );
#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPprintRow(scip, row, nullptr) );
#endif

         // note: cut may not be violated if we are separating a solution different from the LP solution!
         if ( SCIPisFeasNegative(scip, SCIPgetRowLPFeasibility(scip, row)) )
         {
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &cutoff) );
            if ( cutoff )
               break;
            ++ngen;
         }

         SCIP_CALL( SCIPreleaseRow(scip, &row) );

         if ( ngen >= maxgen )
            break;
      }

      if ( cutoff )
         break;
   }

   // free storage
   SCIPfreeBufferArray(scip, &tmpvars);
   SCIPfreeBufferArray(scip, &tmpvals);

   return SCIP_OKAY;
}



//--------------------------------------------------------------------------------------------
//--------------------------------- SCIP callback functions ----------------------------------
//--------------------------------------------------------------------------------------------

//! initialization method of cut separator
SCIP_DECL_SEPAINIT(BACSinitSepaNeigh)
{
   assert( scip != nullptr );
   assert( sepa != nullptr );

   SCIPdebugMsg(scip, "Initialize separator <%s>.\n", SEPA_NAME);

   // get separator data
   SCIP_SEPADATA* sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != nullptr );

   // get problem data
   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   assert( sepadata->pool == nullptr );
   sepadata->pool = new vectorpool(probdata->n + 1000, 10000);

   // generate neighborhood inequalities if separation is turned on
   int freq = -1;
   SCIP_CALL( SCIPgetIntParam(scip, "separating/neigh/freq", &freq) );

   if ( freq >= 0 )
   {
      unsigned int ngen;
      unsigned int ntested;
      struct tms timer_beg;
      struct tms timer_end;

      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Generating neighborhoods ...\n");

      (void) times(&timer_beg);
      SCIP_CALL( generateNeighborhoodInequalities(scip, sepadata->maxneighbors, sepadata->maxneighfract, sepadata->pool, sepadata->nodelimit, ngen, ntested) );
      (void) times(&timer_end);

      double t = (timer_end.tms_utime - timer_beg.tms_utime) / (double)sysconf(_SC_CLK_TCK);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Generated neighborhoods:\t\t%u\t(%u tested from %zu)\n", ngen, ntested, probdata->n);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Time for generating neighborhoods:\t%4.2f\n\n", t);
   }

   return SCIP_OKAY;
}


//! deinitialization method of cut separator
SCIP_DECL_SEPAEXIT(BACSexitSepaNeigh)
{
   assert( scip != nullptr );
   assert( sepa != nullptr );

   SCIPdebugMsg(scip, "exit method of separator <%s>\n", SEPA_NAME);

   // get separator data
   SCIP_SEPADATA* sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != nullptr );

   assert( sepadata->pool != nullptr );
   if ( SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_MINIMAL && sepadata->pool->size() > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\nTotal number of neighborhood subsets found:\t%8zu\n", sepadata->pool->size());

      // get problem data
      SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
      assert( probdata != nullptr );

      // get sizes and variables
      size_t n = probdata->n;

      // get size statistics
      size_t* nsize = nullptr;
      SCIP_CALL( SCIPallocBufferArray(scip, &nsize, n) );

      // print statistics
      if ( sepadata->pool->size() > 0 )
      {
         for (unsigned int i = 0; i < n; ++i)
            nsize[i] = 0;
         size_t maxsize = 0;
         vectorpool::const_iterator vit, vend = sepadata->pool->end();
         for (vit = sepadata->pool->begin(); vit != vend; ++vit)
         {
            size_t size = ((*vit)->size() - 1)/2;
            assert( size <= n );
            ++nsize[size];
            if ( size > maxsize )
               maxsize = size;
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Neighborhood set sizes:\n");
         for (size_t i = 4; i <= maxsize; ++i)
         {
            if ( nsize[i] > 0 )
            {
               int w = (int) ceil(log10(i+0.1));
               w = MAX(w, (int) ceil(log10(nsize[i]+0.1)));
               SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "%*zu ", w, i);
            }
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\n");
         for (size_t i = 4; i <= maxsize; ++i)
         {
            if ( nsize[i] > 0 )
            {
               int w = (int) ceil(log10(i+0.1));
               w = MAX(w, (int) ceil(log10(nsize[i]+0.1)));
               SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "%*zu ", w, nsize[i]);
            }
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\n");
      }

      SCIPfreeBufferArray(scip, &nsize);
   }

   delete sepadata->pool;
   sepadata->pool = nullptr;

   // free sepadata
   SCIPfreeBlockMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, nullptr);

   return SCIP_OKAY;
}


//! separate inequalities
SCIP_DECL_SEPAEXECLP(BACSexeclpSepaNeigh)
{
   assert( scip != nullptr );
   assert( sepa != nullptr );
   assert( result != nullptr );

   *result = SCIP_DIDNOTRUN;

   // only call separator, if we are not close to terminating
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   // only call separator, if an optimal LP solution is at hand
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   // only start if there are fractional values for integer variables
   if ( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   // get separator data
   SCIP_SEPADATA* sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != nullptr );

   // check pool
   assert( sepadata->pool != nullptr );
   if ( sepadata->pool->size() > 0 )
   {
      unsigned int ngen = 0;
      SCIP_Bool cutoff;

      *result = SCIP_DIDNOTFIND;

      SCIP_CALL( separatePools(scip, sepa, nullptr, sepadata->pool, ngen, 100, cutoff) );
      SCIPdebugMsg(scip, "Separated neighborhood cuts from the pool: %u [%zu].\n", ngen, sepadata->pool->size());

      if ( cutoff )
         *result = SCIP_CUTOFF;
      else if ( ngen > 0 )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}


//! arbitrary primal solution separation method of separator
SCIP_DECL_SEPAEXECSOL(BACSexecsolSepaNeigh)
{
   assert( scip != nullptr );
   assert( sepa != nullptr );
   assert( result != nullptr );

   *result = SCIP_DIDNOTRUN;

   // only call separator, if we are not close to terminating
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   // only start if there are fractional values for integer variables
   if ( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   // get separator data
   SCIP_SEPADATA* sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != nullptr );

   // check pool
   assert( sepadata->pool != nullptr );
   if ( sepadata->pool->size() > 0 )
   {
      unsigned int ngen = 0;
      SCIP_Bool cutoff;

      *result = SCIP_DIDNOTFIND;

      SCIP_CALL( separatePools(scip, sepa, sol, sepadata->pool, ngen, 100, cutoff) );
      SCIPdebugMsg(scip, "Separated neighborhood cuts from the pool: %u [%zu].\n", ngen, sepadata->pool->size());

      if ( cutoff )
         *result = SCIP_CUTOFF;
      else if ( ngen > 0 )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}


//! creates the separator for neighborhood constraints and includes it in SCIP
SCIP_RETCODE BACSincludeSepaNeigh(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata = nullptr;
   SCIP_SEPA* sepa;

   // create separator data
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   sepadata->pool = nullptr;

   // include separator
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESUBSCIP, SEPA_DELAY, BACSexeclpSepaNeigh, BACSexecsolSepaNeigh, sepadata) );
   assert( sepa != nullptr );

   // set non-nullptr pointers to callback methods
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, BACSinitSepaNeigh) );
   SCIP_CALL( SCIPsetSepaExit(scip, sepa, BACSexitSepaNeigh) );

   // add parameters
   SCIP_CALL_ABORT( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxneighbors",
         "maximal number of nodes in neighborhood to compute stable set for",
         &sepadata->maxneighbors, TRUE, DEFAULT_MAXNEIGHBORS, 0, INT_MAX, nullptr, nullptr) );

   SCIP_CALL_ABORT( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/maxneighfract",
         "maximal fraction of number of nodes present in neighborhood",
         &sepadata->maxneighfract, TRUE, DEFAULT_MAXNEIGHFRACT, 0.0, 1.0, nullptr, nullptr) );

   SCIP_CALL_ABORT( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/nodelimit",
         "nodelimit for branch-and-bound in subscip",
         &sepadata->nodelimit, TRUE, DEFAULT_NODELIMIT, 0, INT_MAX, nullptr, nullptr) );

   return SCIP_OKAY;
}
