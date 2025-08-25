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

/**@file   prop_cliquefixing.cpp
 * @ingroup PROPAGATORS
 * @brief  propagator to fix nodes based on a cliquepartition
 * @author Annika Jaeger
 * @author Jonas Alker
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "prop_cliquefixing.h"
#include "graph.h"
#include "struct_probdata.h"
#include "cliquepartition.h"

#include <scip/prop_symmetry.h>


#define PROP_NAME                  "cliquefixing"
#define PROP_DESC                  "propagator to fix nodes based on a cliquepartition"
#define PROP_TIMING                 SCIP_PROPTIMING_BEFORELP
#define PROP_PRIORITY              -1000000  //!< priority of the propagator (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PROP_FREQ                         0  //!< propagation frequency
#define PROP_DELAY                     TRUE  //!< should propagation method be delayed, if other propagators found reductions?
#define PROP_PRESOL_PRIORITY       -9999999  //!< priority of the propagator (>= 0: before, < 0: after constraint handlers)
#define PROP_PRESOL_MAXROUNDS             0  //!< maximal number of presolving rounds the propagator participates in (-1: no limit)
#define PROP_PRESOLTIMING           SCIP_PRESOLTIMING_FAST //!< timing of the presolving method (fast, medium, or exhaustive)

#define DEFAULT_PROBINGONLY           FALSE  //!< whether the propagator should only run in probing
#define DEFAULT_CONFLICTANALYSIS      FALSE  //!< whether conflict analysis should be applied

#define DEFAULT_MINDENSITY              0.0  //!< minimum density the graph should have such that propagator runs
#define DEFAULT_MAXIMPROVEMENTGAP         3  //!< maximum gap between partition and lowerbound such that for an improvement is searched
#define DEFAULT_MAXFIXINGGAP              5  //!< maximum gap between partition and lowerbound such that for fixings is searched
#define DEFAULT_ONLYCUTOFF            FALSE  //!< whether only for cutoff should be searched

#define DEFAULT_PARTITIONGREEDY        TRUE  //!< use greedy partition
#define DEFAULT_PARTITIONLINEAR       FALSE  //!< use linear partition
#define DEFAULT_PARTITIONCOLORING     FALSE  //!< use coloring partition

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_Bool             probingonly;        //!< whether the propagator should only run in probing
   SCIP_Bool             conflictanalysis;   //!< whether conflict analysis should be applied

   SCIP_Real             mindensity;         //!< minimum density the graph should have such that propagator runs
   int                   maximprovementgap;  //!< maximum gap between partition and lowerbound such that for an improvement is searched
   int                   maxfixinggap;       //!< maximum gap between partition and lowerbound such that for fixings is searched
   SCIP_Bool             onlycutoff;         //!< whether only for cutoff should be searched

   SCIP_Bool             partitiongreedy;    //!< whether greedy partitioning should be used
   SCIP_Bool             partitionlinear;    //!< whether linear partitioning should be used
   SCIP_Bool             partitioncoloring;  //!< whether coloring partitioning should be used
};


/*
 * Local methods
 */


// fix high degree nodes to zero and isolated nodes to one
SCIP_RETCODE fixDegreeNodes(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Bool&            cutoff,             /**< whether we detected a cutoff */
   int&                  nfixed,             /**< number of variables fixed */
   size_t&               nfreenodes,         /**< number of locally unfixed nodes */
   size_t&               lowerbound          //*< local lower bound */
   )
{
   assert( scip != nullptr );
   assert( prop != nullptr );
   assert( probdata != nullptr );
   assert( ! cutoff );
   assert( nfixed == 0 );

   const Graph* G = probdata->G;

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      assert( probdata->localdegrees[v] < 0 || SCIPisGT(scip, SCIPvarGetObj(probdata->vars[v]), 0.0) );

      // skip fixed nodes
      if ( probdata->localdegrees[v] == -1 )
         continue;

      // fixed isolated nodes to 1
      if ( probdata->localdegrees[v] == 0 )
      {
         assert( SCIPvarGetUbLocal(probdata->vars[v]) > 0.5 );
         assert( SCIPvarGetLbLocal(probdata->vars[v]) < 0.5 );

#ifndef NDEBUG
         // all neighbors should be fixed to zero
         AdjacencyIterator ait, aend;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            Vertex w = *ait;
            assert( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 );
         }
#endif

         // node is isolated
         SCIP_Bool infeasible;
         SCIP_Bool fixed;
         if ( propdata->conflictanalysis )
            SCIP_CALL( SCIPinferBinvarProp(scip, probdata->vars[v], TRUE, prop, -1, &infeasible, &fixed) );
         else
            SCIP_CALL( SCIPfixVar(scip, probdata->vars[v], 1.0, &infeasible, &fixed) );

         assert( ! infeasible );
         assert( fixed );
         ++nfixed;
         --lowerbound;
         --nfreenodes;

         continue;
      }

      // fix v to 0.0 due to its local degree
      assert( probdata->localdegrees[v] >= 1 );
      if ( nfreenodes - (size_t) probdata->localdegrees[v] + 1 <= lowerbound )
      {
         assert( probdata->localdegrees[v] >= 2 );
         assert( SCIPvarGetLbLocal(probdata->vars[v]) < 0.5 && SCIPvarGetUbLocal(probdata->vars[v]) > 0.5 );

         SCIP_Bool infeasible;
         SCIP_Bool fixed;
         if ( propdata->conflictanalysis )
            SCIP_CALL( SCIPinferBinvarProp(scip, probdata->vars[v], FALSE, prop, 0, &infeasible, &fixed) );
         else
            SCIP_CALL( SCIPfixVar(scip, probdata->vars[v], 0.0, &infeasible, &fixed) );

         assert( ! infeasible );
         assert( fixed );
         ++nfixed;
         --nfreenodes;

         // maybe a cutoff due to fixings has been found
         if ( nfreenodes <= lowerbound )
         {
            cutoff = TRUE;
            SCIPdebugMsg(scip, "Node at depth: %d, cut off, nfreenodes: %lu, lowerbound: %lu\n", SCIPgetDepth(scip), nfreenodes, lowerbound);

            return SCIP_OKAY;
         }
      }
   }
   return SCIP_OKAY;
}


//! perform propagation
static
SCIP_RETCODE propagateCliquefixing(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   PARTITIONING_METHOD   partitioning,       /**< partitioning method */
   SCIP_Real&            upperbound,         /**< valid upperbound for the current node */
   SCIP_Bool&            cutoff,             /**< whether we detected a cutoff */
   int&                  nfixed              /**< number of variables fixed */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   cutoff = FALSE;
   nfixed = 0;
   upperbound = SCIPinfinity(scip);

   size_t nfixedzero = 0;

   size_t lowerbound;
   SCIP_CALL( BACScomputeLocalLowerbound(scip, probdata, lowerbound) );
   if ( lowerbound == 0 )
      return SCIP_OKAY;

   size_t ngen = 0;
   const Graph* G = probdata->G;
   size_t n = probdata->n;

   size_t nfreenodes = n - probdata->nlocalones - probdata->nlocalzeros;
   if ( nfreenodes <= lowerbound )
   {
      cutoff = TRUE;
      SCIPdebugMsg(scip, "Node at depth: %d, cut off, nfreenodes: %lu, lowerbound: %lu\n", SCIPgetDepth(scip), nfreenodes, lowerbound);

      return SCIP_OKAY;
   }

   // fix high degree nodes to zero and isolated nodes to one
   SCIP_CALL( fixDegreeNodes(scip, prop, probdata, cutoff, nfixed, nfreenodes, lowerbound) );

   if ( cutoff )
      return SCIP_OKAY;

   size_t* cliquesizes = probdata->cliquesizes;
   assert( cliquesizes != nullptr );

   // prepare loop
   int* cand = probdata->cliqueindex;
   assert( cand != nullptr );

   if ( partitioning == PARTITIONING_COLOR )
   {
      SCIP_CALL( computeCliquePartitionColoring(scip, probdata, cand, cliquesizes, ngen) );
   }
   else if ( partitioning == PARTITIONING_LINEAR )
   {
      SCIP_CALL( computeCliquePartitionLinear(scip, probdata, cand, cliquesizes, ngen) );
   }
   else if ( partitioning == PARTITIONING_GREEDY )
   {
      SCIP_CALL( computeCliquePartitionGreedy(scip, probdata, cand, cliquesizes, ngen) );
   }

   probdata->maxcliqueindex = ngen;
   probdata->ncliques = ngen;
   probdata->cliquenode = SCIPgetCurrentNode(scip);

   upperbound = probdata->nlocalones + probdata->ncliques + probdata->objoffset;
#ifndef NDEBUG
   BACScheckCliquesizes(scip, cand, cliquesizes, ngen, n);
#endif

   // check whether we can prune current node
   if ( ngen <= lowerbound )
   {
      cutoff = TRUE;
      SCIPdebugMsg(scip, "Node at depth: %d, cut off due to cliquecover, ngen: %lu, lowerbound: %lu\n", SCIPgetDepth(scip), ngen, lowerbound);

      return SCIP_OKAY;
   }

   size_t* nnodesinclique;
   SCIP_CALL( SCIPallocBufferArray(scip, &nnodesinclique, ngen) );

   size_t& ncliques = probdata->ncliques;

   // if partitioning is close enough to being able to cut current node we search for an improvement
   if ( ncliques - lowerbound <= (size_t) propdata->maximprovementgap )
      SCIP_CALL( improvePartition(scip, probdata, cand, cliquesizes, ngen, lowerbound, ncliques, cutoff) );

   // stop if we can prune now or only want to check for cutoff
   if ( cutoff || propdata->onlycutoff )
   {
      SCIPfreeBufferArray(scip, &nnodesinclique);
      upperbound = probdata->nlocalones + probdata->ncliques + probdata->objoffset;

      return SCIP_OKAY;
   }

   // if partitioning is too weak, we dont look for fixings
   if ( ncliques - lowerbound > (size_t) propdata->maxfixinggap )
   {
      SCIPfreeBufferArray(scip, &nnodesinclique);
      upperbound = probdata->nlocalones + probdata->ncliques + probdata->objoffset;

      return SCIP_OKAY;
   }

   // fix variables to zero if they are contained in more cliques than ncliques - LB + 1 (choosing this node can never improve LB)
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      if ( ncliques <= lowerbound )
      {
         cutoff = TRUE;
         nfixed = 0;
         nfixedzero = 0;
         SCIPdebugMsg(scip, "Node at depth: %d, cut off due to cliquecover, ncliques: %lu, lowerbound: %lu\n", SCIPgetDepth(scip), ncliques, lowerbound);

         break;
      }

      Vertex v = *vit;

      if ( cand[v] < 0 )
         continue;

      assert( SCIPvarGetLbLocal(probdata->vars[v]) < 0.5 );

      // it might happen that nodes in cliques are fixed to zero due to a 1-fixed neighbor
      if ( SCIPvarGetUbLocal(probdata->vars[v]) < 0.5 )
      {
         if ( cliquesizes[cand[v]] == 1 )
            --ncliques;

         --cliquesizes[cand[v]];
         cand[v] = -1;

         // fixing to zero might lead to a cutoff
         if ( ncliques <= lowerbound )
         {
            cutoff = TRUE;
            nfixed = 0;
            nfixedzero = 0;
            SCIPdebugMsg(scip, "Node at depth: %d, cut off due to cliquecover, ncliques: %lu, lowerbound: %lu\n", SCIPgetDepth(scip), ncliques, lowerbound);

            break;
         }
         continue;
      }

      if ( probdata->localdegrees[v] == 0 )
      {
         assert( SCIPvarGetUbLocal(probdata->vars[v]) > 0.5 );
         assert( SCIPvarGetLbLocal(probdata->vars[v]) < 0.5 );

#ifndef NDEBUG
         // all neighbors should be fixed to zero
         AdjacencyIterator ait, aend;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            Vertex w = *ait;
            assert( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 );
         }
#endif
         // node is isolated
         SCIP_Bool infeasible;
         SCIP_Bool fixed;
         if ( propdata->conflictanalysis )
            SCIP_CALL( SCIPinferBinvarProp(scip, probdata->vars[v], TRUE, prop, -1, &infeasible, &fixed) );
         else
            SCIP_CALL( SCIPfixVar(scip, probdata->vars[v], 1.0, &infeasible, &fixed) );

         assert( ! infeasible );
         assert( fixed );

         assert( cliquesizes[cand[v]] == 1 );
         --cliquesizes[cand[v]];
         cand[v] = -1;
         --lowerbound;
         --ncliques;

         ++nfixed;
         SCIPdebugMsg(scip, "Node at depth: %d, fix node %lu to 1.0, as it is isolated.\n", SCIPgetDepth(scip), v);

         continue;
      }

      assert( ncliques > lowerbound );

      // if ncliques = LB + 1, we fix cliques of size 1 to 1
      if ( ncliques == lowerbound + 1 && cliquesizes[cand[v]] == 1 )
      {
         assert( SCIPvarGetUbLocal(probdata->vars[v]) > 0.5 );
         assert( SCIPvarGetLbLocal(probdata->vars[v]) < 0.5 );

         SCIP_Bool infeasible;
         SCIP_Bool fixed;
         if ( propdata->conflictanalysis )
            SCIP_CALL( SCIPinferBinvarProp(scip, probdata->vars[v], TRUE, prop, 0, &infeasible, &fixed) );
         else
            SCIP_CALL( SCIPfixVar(scip, probdata->vars[v], 1.0, &infeasible, &fixed) );

         assert( ! infeasible );
         assert( fixed );

         --cliquesizes[cand[v]];
         cand[v] = -1;
         --lowerbound;
         --ncliques;

         ++nfixed;
         SCIPdebugMsg(scip, "Node at depth: %d, fix node %lu to 1.0, as ncliques = LB + 1, ncliques: %lu\n", SCIPgetDepth(scip), v, ncliques);

         // also fix neighbors
         AdjacencyIterator ait, aend;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            Vertex w = *ait;
            assert( SCIPvarGetLbLocal(probdata->vars[w]) < 0.5 );

            if ( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 )
               continue;

            if ( propdata->conflictanalysis )
               SCIP_CALL( SCIPinferBinvarProp(scip, probdata->vars[w], FALSE, prop, 1, &infeasible, &fixed) );
            else
               SCIP_CALL( SCIPfixVar(scip, probdata->vars[w], 0.0, &infeasible, &fixed) );

            assert( ! infeasible );
            assert( fixed );

            if ( cliquesizes[cand[w]] == 1 )
               --ncliques;

            --cliquesizes[cand[w]];
            cand[w] = -1;

            ++nfixed;
            ++nfixedzero;
         }
#ifndef NDEBUG
         BACScheckCliquesizes(scip, cand, cliquesizes, ngen, n);
#endif
         continue;
      }

      // break if remaining graph is empty or improvement has been found due to the fixings
      if ( lowerbound == 0 )
         break;

      // v needs to be added in at least ngen - LB + 1 cliques.
      assert( probdata->localdegrees[v] >= 0 );
      if ( (size_t) probdata->localdegrees[v] < ncliques - lowerbound + 1 )
         continue;

      // count in how many cliques v could be added
      size_t ncoveredcliques = 1;
      BMSclearMemoryArray(nnodesinclique, ngen);

      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         if ( cand[w] == cand[v] )
            continue;

         if ( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[w]) > 0.5 )
            continue;

         assert( cand[w] >= 0 );
         assert( cand[w] < (int) ngen );

         ++nnodesinclique[cand[w]];

         assert( cliquesizes[cand[w]] > 0 );

         // if clique is covered by neighborhood of v
         if ( nnodesinclique[cand[w]] == cliquesizes[cand[w]] )
            ++ncoveredcliques;

         if ( ncoveredcliques >= ncliques - lowerbound + 1 )
         {  // fix v to zero, neighborhood of v covers to many cliques
            assert( ncoveredcliques > 1 );

            SCIP_Bool infeasible;
            SCIP_Bool fixed;
            if ( propdata->conflictanalysis )
               SCIP_CALL( SCIPinferBinvarProp(scip, probdata->vars[v], FALSE, prop, 0, &infeasible, &fixed) );
            else
               SCIP_CALL( SCIPfixVar(scip, probdata->vars[v], 0.0, &infeasible, &fixed) );

            assert( ! infeasible );
            assert( fixed );

            ++nfixed;
            ++nfixedzero;
            SCIPdebugMsg(scip, "Node at depth: %d, fix node %lu to 0.0, as it can be included in %lu of %lu cliques, while LB = %lu\n", SCIPgetDepth(scip), v, ncoveredcliques, ncliques, lowerbound);

            if ( cliquesizes[cand[v]] == 1 )
               --ncliques;

            --cliquesizes[cand[v]];
            cand[v] = -1;

            break;
         }
      }
#ifndef NDEBUG
      BACScheckCliquesizes(scip, cand, cliquesizes, ngen, n);
#endif
   }

   // only try to improve partition again if fixings to zero have been found and gap is not too large
   if ( nfixedzero > 0 && ncliques - lowerbound < (size_t) propdata->maximprovementgap )
      SCIP_CALL( improvePartition(scip, probdata, cand, cliquesizes, ngen, lowerbound, ncliques, cutoff) );

   if ( cutoff )
      nfixed = 0;

   SCIPfreeBufferArray(scip, &nnodesinclique);
   upperbound = probdata->nlocalones + probdata->ncliques + probdata->objoffset; /*lint !e838*/

   return SCIP_OKAY;
}

/** conflict analysis initiation method of cliquefixing propagator
 *
 * only vars fixed to zero are added to conflict
 * one-fixings are not relevant for cutoff:
 *    - considering only the zero fixings, fixed to 1 vertices are isolated
 *    - in a clique partition they would give one extra clique of size one each, increasing ngen by the number of 1-fixings
 *    - lowerbound would increase by the number of 1-fixings either
 *    - if ngen == lowerbound, then also ngen + |1-fixings| = lowerbound + |1-fixings|
 */
SCIP_RETCODE conflictCliquefixing(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( SCIPisConflictAnalysisApplicable(scip) );

   for (size_t i = 0; i < probdata->n; ++i)
   {
      SCIP_VAR* var = probdata->vars[i];

      // add all vars fixed to zero to conflict set
      if ( SCIPvarGetUbLocal(var) < 0.5 )
         SCIP_CALL( SCIPaddConflictBinvar(scip, var) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */


//! destructor of propagator to free user data (called when SCIP is exiting)
static
SCIP_DECL_PROPFREE(propFreeCliquefixing)
{
   assert( prop != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   SCIPfreeBlockMemory(scip, &propdata);

   return SCIP_OKAY;
}


//! presolving method
static
SCIP_DECL_PROPPRESOL(propPresolCliquefixing)
{
   assert( scip != nullptr );
   assert( result != nullptr );
   assert( nfixedvars != nullptr );

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   // set result
   *result = SCIP_DIDNOTRUN;

   // only works un unweighted graphs
   if ( ! probdata->unweighted )
      return SCIP_OKAY;

   // do not run if disabled
   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   if ( ! propdata->partitioncoloring && ! propdata->partitionlinear && ! propdata->partitiongreedy )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "running cliquefixing presolving.\n");

   SCIP_Bool cutoff = FALSE;
   int nfixedcolor = 0;
   int nfixedlinear = 0;
   int nfixedgreedy = 0;
   SCIP_Real upperbound = SCIPinfinity(scip);

   if ( propdata->partitioncoloring )
   {
      SCIP_CALL( propagateCliquefixing(scip, prop, probdata, PARTITIONING_COLOR, upperbound, cutoff, nfixedcolor) );
   }
   if ( ! cutoff && propdata->partitionlinear )
   {
      SCIP_CALL( propagateCliquefixing(scip, prop, probdata, PARTITIONING_LINEAR, upperbound, cutoff, nfixedlinear) );
   }
   if ( ! cutoff && propdata->partitiongreedy )
   {
      SCIP_CALL( propagateCliquefixing(scip, prop, probdata, PARTITIONING_GREEDY, upperbound, cutoff, nfixedgreedy) );
   }

   // cutoff in presolving implies optimality of already found solution
   if ( cutoff )
   {
      *result = SCIP_CUTOFF;
   }
   else if ( nfixedcolor + nfixedgreedy + nfixedlinear > 0 )
   {
      *result = SCIP_SUCCESS;
      *nfixedvars += nfixedcolor + nfixedgreedy + nfixedlinear;
   }
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY; /*lint !e438*/
}


//! execution method of propagator
static
SCIP_DECL_PROPEXEC(propExecCliquefixing)
{  /*lint --e{715}*/
   assert( scip != nullptr );
   assert( result != nullptr );
   assert( prop != nullptr );

   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   // set result
   *result = SCIP_DIDNOTRUN;

   // only works un unweighted graphs
   if ( ! probdata->unweighted )
      return SCIP_OKAY;

   // if propagation should only run in probing
   if ( propdata->probingonly )
   {
      if ( ! SCIPinProbing(scip) || SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING )
         return SCIP_OKAY;
   }

   if ( SCIPisLT(scip, probdata->localdensity, propdata->mindensity) )
      return SCIP_OKAY;

   if ( ! propdata->partitioncoloring && ! propdata->partitionlinear && ! propdata->partitiongreedy )
      return SCIP_OKAY;

   if ( SCIPinProbing(scip) )
      SCIPdebugMsg(scip, "running cliquefixing propagation in probing.\n");
   else
      SCIPdebugMsg(scip, "running cliquefixing propagation.\n");

   SCIP_Bool cutoff = FALSE;
   int nfixedcolor = 0;
   int nfixedlinear = 0;
   int nfixedgreedy = 0;
   SCIP_Real upperbound = SCIPinfinity(scip);

   if ( propdata->partitioncoloring )
   {
      SCIP_CALL( propagateCliquefixing(scip, prop, probdata, PARTITIONING_COLOR, upperbound, cutoff, nfixedcolor) );
   }
   if ( ! cutoff && propdata->partitionlinear )
   {
      SCIP_CALL( propagateCliquefixing(scip, prop, probdata, PARTITIONING_LINEAR, upperbound, cutoff, nfixedlinear) );
   }
   if ( ! cutoff && propdata->partitiongreedy )
   {
     SCIP_CALL( propagateCliquefixing(scip, prop, probdata, PARTITIONING_GREEDY, upperbound, cutoff, nfixedgreedy) );
   }

   if ( cutoff )
   {
      *result = SCIP_CUTOFF;

      // perform conflict analysis
      if ( propdata->conflictanalysis && SCIPisConflictAnalysisApplicable(scip) )
      {
         SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );
         SCIP_CALL( conflictCliquefixing(scip, probdata ) );
         SCIP_CALL( SCIPanalyzeConflict(scip, 0, nullptr) );
      }
   }
   else if ( nfixedcolor + nfixedgreedy + nfixedlinear > 0 )
   {
      *result = SCIP_REDUCEDDOM;
      SCIP_CALL( SCIPupdateLocalDualbound(scip, upperbound) );
   }
   else
   {
      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( SCIPupdateLocalDualbound(scip, upperbound) );
   }

   return SCIP_OKAY;
}

/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropCliquefixing)
{
   assert( scip != nullptr );
   assert( result != nullptr );
   assert( prop != nullptr );
   assert( infervar != nullptr );
   assert( bdchgidx != nullptr );

   *result = SCIP_DIDNOTFIND;

   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   Graph G = *probdata->G;

   assert( SCIPisConflictAnalysisApplicable(scip) );
   assert( inferinfo >= -1 );
   assert( inferinfo <= 1 );

   Vertex u = (Vertex) SCIPvarGetData(infervar);
   assert( u < probdata->n );

   if ( inferinfo == 0 )
   {
      // node was fixed to zero due to high local degree or based fixed based on clique partition
      // add all vars fixed to zero for non-adjacent vertices

      SCIP_Bool* neighbors;
      SCIP_CALL( SCIPallocClearBufferArray(scip, &neighbors, probdata->n) );

      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(u, G); ait != aend; ++ait)
      {
         Vertex v = *ait;
         assert( v < probdata->n );

         neighbors[v] = TRUE;
      }

      for (size_t i = 0; i < probdata->n; ++i)
      {
         SCIP_VAR* var = probdata->vars[i];

         // skip u and neighbors of u
         if ( i == u || neighbors[i] )
            continue;

         // add all vars fixed to zero to conflict set
         if ( SCIPvarGetUbAtIndex(var, bdchgidx, FALSE) < 0.5 )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, var) );
            *result = SCIP_SUCCESS;
         }
      }
      SCIPfreeBufferArray(scip, &neighbors);
   }
   else if ( inferinfo == -1 )
   {
      // node was isolated and fixed to one
      // add all neighbors to conflict
      assert( boundtype == SCIP_BOUNDTYPE_LOWER );

      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(u, G); ait != aend; ++ait)
      {
         Vertex v = *ait;
         assert( v < probdata->n );

         SCIP_VAR* var = probdata->vars[v];
         assert( SCIPgetVarUbAtIndex(scip, var, bdchgidx, FALSE) < 0.5 );

         SCIP_CALL( SCIPaddConflictBinvar(scip, var) );
         *result = SCIP_SUCCESS;
      }
   }
   else if ( inferinfo == 1 )
   {
      // node was fixed to zero because of a 1-fixed neighbor
      // add neighbor fixed to one to conflict
      assert( boundtype == SCIP_BOUNDTYPE_UPPER );

      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(u, G); ait != aend; ++ait)
      {
         Vertex v = *ait;
         if ( SCIPgetVarLbAtIndex(scip, probdata->vars[v], bdchgidx, TRUE) > 0.5 )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, probdata->vars[v]) );
            *result = SCIP_SUCCESS;
            break;
         }
      }
   }

   return SCIP_OKAY;
}

/*
 * propagator specific interface methods
 */

//! creates the cliquefixing propagator and includes it in SCIP
SCIP_RETCODE BACSincludePropCliquefixing(
   SCIP*                 scip                //! SCIP data structure
   )
{
   SCIP_PROPDATA* propdata = nullptr;
   SCIP_PROP* prop = nullptr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );

   // include propagator
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecCliquefixing, propdata) );
   assert( prop != nullptr );

   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolCliquefixing, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS,
         PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeCliquefixing) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropCliquefixing) );

   // add parameters
   SCIP_CALL( SCIPaddBoolParam(scip,
      "propagating/" PROP_NAME "/probingonly",
      "whether the propagator should only run in probing",
      &propdata->probingonly, FALSE, DEFAULT_PROBINGONLY, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "propagating/" PROP_NAME "/conflictanalysis",
      "whether conflict analysis should be applied",
      &propdata->conflictanalysis, FALSE, DEFAULT_CONFLICTANALYSIS, nullptr, nullptr) );

   SCIP_CALL( SCIPaddRealParam(scip,
      "propagating/" PROP_NAME "/mindensity",
      "minimum density the graph should have such that propagator runs",
      &propdata->mindensity, FALSE, DEFAULT_MINDENSITY, 0.0, 1.0, nullptr, nullptr) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "propagating/" PROP_NAME "/maximprovementgap",
      "maximum gap between partition and lowerbound such that for an improvement is searched",
      &propdata->maximprovementgap, FALSE, DEFAULT_MAXIMPROVEMENTGAP, 0, INT_MAX, nullptr, nullptr) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "propagating/" PROP_NAME "/maxfixinggap",
      "maximum gap between partition and lowerbound such that for fixings is searched",
      &propdata->maxfixinggap, FALSE, DEFAULT_MAXFIXINGGAP, 0, INT_MAX, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "propagating/" PROP_NAME "/partitiongreedy",
      "whether greedy partitioning should be used",
      &propdata->partitiongreedy, FALSE, DEFAULT_PARTITIONGREEDY, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "propagating/" PROP_NAME "/partitionlinear",
      "whether linear partitioning should be used",
      &propdata->partitionlinear, FALSE, DEFAULT_PARTITIONLINEAR, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "propagating/" PROP_NAME "/partitioncoloring",
      "whether coloring partitioning should be used",
      &propdata->partitioncoloring, FALSE, DEFAULT_PARTITIONCOLORING, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "propagating/" PROP_NAME "/onlycutoff",
      "whether only for cutoff should be searched",
      &propdata->onlycutoff, FALSE, DEFAULT_ONLYCUTOFF, nullptr, nullptr) );

   return SCIP_OKAY;
}
