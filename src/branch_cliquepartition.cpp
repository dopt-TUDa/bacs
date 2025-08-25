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
/*    Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                    */
/*                                                                           */
/*    Both are licensed under the Apache License, Version 2.0.               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_cliquepartition.cpp
 * @brief  branching on a clique partition (a.k.a. Balas-Yu) rule
 * @author Annika Jaeger
 * @author Erik Jansen
 * @author Jonas Alker
 * @author Marc Pfetsch
 */

#include <assert.h>

#include "branch_cliquepartition.h"
#include "struct_probdata.h"
#include "cliquepartition.h"

#define BRANCHRULE_NAME            "cliquepartition"
#define BRANCHRULE_DESC            "clique partition branching rules, follows Balas and Yus, Sewell in induced subgraph of unfixed nodes"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define DEFAULT_MINDEPTH           0
#define DEFAULT_MINDENSITY         0.0
#define DEFAULT_MAXABSGAP          INT_MAX
#define DEFAULT_FRACNODECHOICE     0

/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   int                   mindepth;           //< mindepth parameter
   SCIP_Real             mindensity;         //< minimum density such that branching rule is applied
   int                   maxabsgap;          //< maximum absolute solution gap parameter
   int                   fracnodechoice;     //< method of choosinf fracnode
                                             /* 0=in smalles clique with fractional node
                                              * 1=in largest clique with fractional node
                                              * 2=node with larges unfixed neigborhood
                                              */
};


/*
 * Local methods
 */


/** try to reduce number of branches, decrease size of not marked cliques by moving nodes to marked cliques */
static
SCIP_RETCODE reduceNotMarkedCliques(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int*                  cand,               /**< clique partition, representant for every vertex */
   size_t*               cliquesizes,        /**< sizes of cliques */
   size_t                ngen,               /**< maximal clique index (initally generated cliques) */
   size_t&               ncliques,           /**< number of cliques */
   const SCIP_Bool*      marked,             /**< array of marked cliques */
   size_t                ignoreclique        /**< index of clique to be completely ignored */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( cand != nullptr );
   assert( cliquesizes != nullptr );
   assert( marked != nullptr );

   const Graph* G = probdata->G;
   size_t nmoved = 0;

   // init counter
   size_t* nnodesinclique;
   SCIP_CALL( SCIPallocBufferArray(scip, &nnodesinclique, probdata->n) );

   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;

      // skip non covered nodes
      if ( cand[v] == -1 )
         continue;

      // only move not marked vertices
      if ( marked[cand[v]] )
         continue;

      // ignore clique
      if ( cand[v] == (int) ignoreclique )
         continue;

      // check whether v can be added to existing marked clique
      BMSclearMemoryArray(nnodesinclique, ngen);

      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         if ( cand[w] == -1 )
            continue;
         assert( cand[w] >= 0 );
         assert( cand[w] < (int) ngen );

         ++(nnodesinclique[cand[w]]);

         // if cand[w] is covered by neighborhood of v we can add v to this clique
         if ( cand[w] != (int) ignoreclique && marked[cand[w]] && nnodesinclique[cand[w]] == cliquesizes[cand[w]] )
         {
            --cliquesizes[cand[v]];
            if ( cliquesizes[cand[v]] == 0 )
               --ncliques;

            cand[v] = cand[w];
            ++(cliquesizes[cand[w]]);
            ++nmoved;
            break;
         }
      }
#ifndef NDEBUG
      BACScheckCliquesizes(scip, cand, cliquesizes, ngen, probdata->n);
#endif
   }
   SCIPfreeBufferArray(scip, &nnodesinclique);
   if ( nmoved > 0 )
      SCIPdebugMsg(scip, "Depth %d: Moved %lu nodes to marked cliques.\n", SCIPgetDepth(scip), nmoved);
   return SCIP_OKAY;
}

/** get a valid clique partition either by copying a previously computed one or by creating a new one */
static
SCIP_RETCODE getCliquePartition(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int*&                 cand,               /**< clique partition, representant for every vertex */
   size_t*&              cliquesizes,        /**< sizes of cliques */
   size_t&               ngen,               /**< maximal clique index (initally generated cliques) */
   size_t&               ncliques,           /**< number of cliques */
   std::vector<bool>&    fracclique,         /**< indicator whether fractional node included in clique */
   SCIP_Bool             islpbranching       /**< is lp branching? */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );

   const Graph* G = probdata->G;

   // copy partition from probdata
   cand = probdata->cliqueindex;
   cliquesizes = probdata->cliquesizes;
   ngen = probdata->maxcliqueindex;
   ncliques = probdata->ncliques;

   // make sure partition is valid
   if ( SCIPgetCurrentNode(scip) == probdata->cliquenode )
   {
      // init counter
      size_t* nnodesinclique;
      SCIP_CALL( SCIPallocBufferArray(scip, &nnodesinclique, probdata->n) );

      VertexIterator vit, vend;
      for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
      {
         Vertex v = *vit;

         // var is fixed and not included in partition
         if ( probdata->localdegrees[v] == -1 && cand[v] == -1 )
            continue;

         // var is unfixed and included in partition
         if ( probdata->localdegrees[v] >= 0 && cand[v] >= 0 )
            continue;

         // partition includes fixed vertex
         if ( probdata->localdegrees[v] == -1 )
         {
            assert( cand[v] >= 0 );

            --cliquesizes[cand[v]];
            if ( cliquesizes[cand[v]] == 0 )
               --ncliques;
            cand[v] = -1;
            continue;
         }

         // partition does not include unfixed vertex
         if ( cand[v] == -1 )
         {
            assert( probdata->localdegrees[v] >= 0 );

            // check whether v can be added to existing nonempty clique
            BMSclearMemoryArray(nnodesinclique, ngen);

            AdjacencyIterator ait, aend;
            for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
            {
               Vertex w = *ait;

               if ( cand[w] == -1 )
                  continue;
               assert( cand[w] >= 0 );
               assert( cand[w] < (int) ngen );

               ++(nnodesinclique[cand[w]]);

               // if cand[w] is covered by neighborhood of v we can add v to this clique
               if ( nnodesinclique[cand[w]] == cliquesizes[cand[w]] )
               {
                  assert( cand[v] == -1 );
                  cand[v] = cand[w];
                  ++(cliquesizes[cand[w]]);
                  break;
               }
            }
            // if v could not be added to existing clique, create new one here
            if ( cand[v] == -1 )
            {
               cand[v] = (int) ngen;
               cliquesizes[ngen] = 1;
               ++ngen;
               ++ncliques;
            }
         }
      }
#ifndef NDEBUG
      BACScheckCliquesizes(scip, cand, cliquesizes, ngen, probdata->n);
#endif
      SCIPfreeBufferArray(scip, &nnodesinclique);
      SCIPdebugMsg(scip, "Depth %d: Clique partition branching using previously stored partition of size %lu.\n", SCIPgetDepth(scip), ncliques);
   }
   else
   {
      // compute new partition
      SCIP_CALL( computeCliquePartitionGreedy(scip, probdata, cand, cliquesizes, ngen) );
      ncliques = ngen;

      SCIPdebugMsg(scip, "Depth %d: Clique partition branching computed partition of size %lu.\n", SCIPgetDepth(scip), ncliques);
   }

   // update probdata
   probdata->maxcliqueindex = ngen;
   probdata->ncliques = ncliques;
   probdata->cliquenode = SCIPgetCurrentNode(scip);

   // compute fractionalcliques
   if ( islpbranching )
   {
      VertexIterator vit, vend;
      for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
      {
         Vertex v = *vit;
         SCIP_VAR* var = probdata->vars[v];

         // skip fixed nodes
         if ( cand[v] < 0 )
            continue;

         assert( cand[v] < (int) ngen );
         assert( SCIPvarGetLbLocal(var) < 0.5 && SCIPvarGetUbLocal(var) > 0.5 );

         if ( ! SCIPisIntegral(scip, SCIPgetSolVal(scip, nullptr, var)) )
            fracclique[cand[v]] = true; /*lint !e732*/
      }
   }

   return SCIP_OKAY;
}


/** perform branching */
static
SCIP_RETCODE branchPartition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< data of branching rule */
   SCIP_RESULT*          result,             /**< result pointer */
   SCIP_Bool             islpbranching       /**< perform LP branching */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( branchruledata != nullptr );
   assert( result != nullptr );
   assert( probdata->n > 0 );

   *result = SCIP_DIDNOTFIND;

   const Graph* G = probdata->G;

   // number unfixed variables
   assert( probdata->n >= probdata->nlocalones + probdata->nlocalzeros );
   std::size_t nn = probdata->n - probdata->nlocalones - probdata->nlocalzeros;

   if ( nn <= 1 )
      return SCIP_OKAY;

   // compute local lower bound
   size_t lowerbound;
   SCIP_CALL( BACScomputeLocalLowerbound(scip, probdata, lowerbound) );

   // cutoff if remaing graph has not enough nodes
   if ( nn <= lowerbound )
   {
      SCIPdebugMsg(scip, "Cutoff: Free graph too small: free nodes: %lu , lowerbound: %lu.\n", nn, lowerbound);
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   // the branching rule is not helpful if the reduced bound is not positive
   if ( lowerbound == 0 )
      return SCIP_OKAY;

   int* cand;
   size_t* cliquesizes;
   size_t ngen = 0;
   size_t ncliques = 0;
   std::vector<bool> fractionalcliques(probdata->n, false);

   // get a valid clique partition
   SCIP_CALL( getCliquePartition(scip, probdata, cand, cliquesizes, ngen, ncliques, fractionalcliques, islpbranching) );

   assert( cand != nullptr );
   assert( cliquesizes != nullptr );
#ifndef NDEBUG
   BACScheckCliquesizes(scip, cand, cliquesizes, ngen, probdata->n);
# endif

   // if the clique partition is smaller than or equal sized as the lower bound, we prune by node
   if ( ncliques <= lowerbound )
   {
      SCIPdebugMsg(scip, "Cutoff due to clique partition: ncliques: %lu, lowerbound: %lu.\n", ncliques, lowerbound);
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   SCIP_Bool* marked = nullptr;
   SCIP_CALL( SCIPallocClearBufferArray(scip, &marked, ngen) );

   // if ngen = LB + 1, we mark cliques of size 1 to 1
   if ( ncliques == lowerbound + 1 )
   {
      for (std::size_t i = 0; i < ngen; i++)
      {
         if ( cliquesizes[i] == 1 )
            marked[i] = TRUE;
      }
   }

   // fix marked cliques and isolated nodes
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;

      // skip already fixed nodes
      if ( cand[v] < 0 )
         continue;

      assert(  cand[v] < (int) ngen );

      // skip non-marked, non-isolated nodes
      if ( ! marked[cand[v]] && probdata->localdegrees[v] > 0 )
         continue;

      // v is unfixed and (isolated or marked)
      assert( SCIPvarGetLbLocal(probdata->vars[v]) < 0.5 && SCIPvarGetUbLocal(probdata->vars[v]) > 0.5 );
      assert( marked[cand[v]] || probdata->localdegrees[v] == 0 );

      SCIP_Bool infeasible;
      SCIP_Bool fixed;
      SCIP_CALL( SCIPfixVar(scip, probdata->vars[v], 1.0, &infeasible, &fixed) );

      assert( ! infeasible );
      assert( fixed );

      *result = SCIP_REDUCEDDOM;

      --cliquesizes[cand[v]];
      --lowerbound;

      assert( cliquesizes[cand[v]] == 0 );
      --ncliques;
      cand[v] = -1;

      // fix adjacent vertices to 0
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         if ( cand[w] < 0 )
            continue;

         assert( SCIPvarGetLbLocal(probdata->vars[w]) < 0.5 );

         SCIP_CALL( SCIPfixVar(scip, probdata->vars[w], 0.0, &infeasible, &fixed) );

         assert( ! infeasible );
         assert( fixed );

         --cliquesizes[cand[w]];
         if ( cliquesizes[cand[w]] == 0 )
            --ncliques;
         cand[w] = -1;
      }
   }
   // check if fixings lead to cutoff
   if ( ncliques <= lowerbound )
   {
      SCIPdebugMsg(scip, "Cutoff after fixing variables: ncliques: %lu , lowerbound: %lu.\n", ncliques, lowerbound);
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }
   // if fixings were found, we do not branch yet
   if ( *result == SCIP_REDUCEDDOM )
   {
      SCIPdebugMsg(scip, "Clique partition found variable fixings: ncliques: %lu , lowerbound: %lu.\n", ncliques, lowerbound);
      SCIPfreeBufferArray(scip, &marked);
      return SCIP_OKAY;
   }

   int* csizes = nullptr;
   int* cliqueorder = nullptr;
   SCIP_CALL( SCIPallocBufferArray(scip, &cliqueorder, ngen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &csizes, ngen) );

   for (size_t i = 0; i < ngen; ++i)
   {
      csizes[i] = (int) cliquesizes[i];
      cliqueorder[i] = (int) i;
   }
   // sort cliques in descending order by size, such that first ncliques entries refer to nonempty cliques
   SCIPsortDownIntInt(csizes, cliqueorder, (int) ngen);
   SCIPfreeBufferArray(scip, &csizes);

   size_t ssqi = INT_MAX;
   size_t* branchclique = nullptr;

   // if branching on LP solution, we search for fractional clique according to selection rule
   if ( islpbranching )
   {
      if ( branchruledata->fracnodechoice == 0 )
      {
         // compute index of vertex on which primal solution is non-integral and containing clique is minimal
         ssqi = ngen + 1;
         for (int i = (int) ncliques - 1; i >= 0; --i)
         {
            size_t cind = (size_t) cliqueorder[i];

            if ( ! fractionalcliques[cind] )
               continue;

            assert( cliquesizes[cind] > 0 );
            assert ( i == 0 || cliquesizes[cliqueorder[i]] <= cliquesizes[cliqueorder[i - 1]] );

            ssqi = cind;
            break;
         }

         assert( ssqi < ngen );
      }
      else if ( branchruledata->fracnodechoice == 1 )
      {
         // compute index of vertex on which primal solution is non-integral and containing clique is maximal
         ssqi = ngen + 1;
         for (std::size_t i = 0; i < ncliques; ++i)
         {
            int cind = cliqueorder[i];

            assert( cind >= 0 );
            if ( ! fractionalcliques[cind] ) /*lint !e732*/
               continue;

            assert ( cliquesizes[cind] > 0 );
            assert ( i == ncliques || cliquesizes[cliqueorder[i]] >= cliquesizes[cliqueorder[i + 1]] );

            ssqi = (size_t) cind;
            break;
         }
         assert( ssqi < ngen );
      }
      else if ( branchruledata->fracnodechoice == 2 )
      {
         // compute index of vertex on which primal solution is non-integral and neighborhood is maximal
         int maxdeg = -1;
         Vertex fracnode = probdata->n;
         for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
         {
            Vertex v = *vit;
            SCIP_Real value = SCIPgetSolVal(scip, nullptr, probdata->vars[v]);

            if ( cand[v] < 0 )
               continue;

            if ( ! SCIPisFeasIntegral(scip, value) )
            {
               assert( probdata->localdegrees[v] > 0 );

               if ( probdata->localdegrees[v] > maxdeg )
               {
                  maxdeg = probdata->localdegrees[v];
                  fracnode = v;
               }
            }
         }
         ssqi = (size_t) cand[fracnode];

         assert( ssqi < ngen );
         assert( fracnode < probdata->n );
      }
      else
      {
         abort();
      }

      // store clique containing fracnode
      SCIP_CALL( SCIPallocBufferArray(scip, &branchclique, cliquesizes[ssqi]) );
      size_t k = 0;
      for ( size_t i = 0; i < probdata->n; ++i )
      {
         if ( k == cliquesizes[ssqi] )
            break;

         if ( cand[i] == (int) ssqi )
         {
            branchclique[k] = i;
            ++k;
         }
      }
   }

   // we mark the first lowerbound many cliques, ignoring ssqi
   size_t nmarked = 0;           // we dont need to reset marked, if it was used, we would have returned earlier
   for (size_t i = 0; i < ncliques; ++i)
   {
      size_t cind = (size_t) cliqueorder[i];

      // stop, when enough cliques have been marked
      if ( nmarked == lowerbound )
         break;

      // ignore ssqi
      if ( cind == ssqi)
         continue;

      marked[cind] = TRUE;
      ++nmarked;
   }

   // in LP branching ssqi should not be marked
   assert( ! islpbranching || ! marked[ssqi] );

   // try to reduce number of branches, decrease size of not marked cliques by moving nodes to marked cliques while ignoring ssqi
   SCIP_CALL( reduceNotMarkedCliques(scip, probdata, cand, cliquesizes, ngen, ncliques, marked, ssqi) );

   // check if swaps lead to cutoff
   if ( ncliques <= lowerbound )
   {
      if ( islpbranching )
         SCIPfreeBufferArray(scip, &branchclique);

      SCIPfreeBufferArray(scip, &cliqueorder);
      SCIPfreeBufferArray(scip, &marked);

      SCIPdebugMsg(scip, "Cutoff after moving vertices to other cliques: ncliques: %lu , lowerbound: %lu.\n", ncliques, lowerbound);
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   Vertex** cliquelist;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cliquelist, ngen) );
   std::vector<size_t> sizes(ngen, 0);
   for (size_t i = 0; i < ngen; ++i)
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(cliquelist[i]), cliquesizes[i]) );

   for (size_t v = 0; v < probdata->n; ++v)
   {
      if ( cand[v] == -1 )
         continue;

      assert( cand[v] >= 0 && cand[v] < (int) ngen );

      cliquelist[cand[v]][sizes[cand[v]]] = v; /*lint !e732*/
      ++sizes[cand[v]]; /*lint !e732*/
   }

#ifndef NDEBUG
   for (size_t i = 0; i < ngen; ++i)
   {
      assert( cliquesizes[i] == sizes[i] );
      for (size_t j = 0; j < cliquesizes[i]; ++j)
         assert( cand[cliquelist[i][j]] == (int) i );
   }
#endif

   // count how many cliques present per node
   size_t dualbound = probdata->nlocalones + lowerbound + (size_t) probdata->objoffset;

   // loop through unmarked nodes
   size_t nbranch = 0;
   for ( size_t i = lowerbound; i < ncliques; ++i )
   {
      size_t cind = (size_t) cliqueorder[i];
      assert( cind < ngen );

      // ignore marked cliques
      if ( marked[cind] )
         continue;

      // ignore fractional clique
      if ( cind == ssqi )
         continue;

      for ( size_t index = 0; index < cliquesizes[cind]; ++index)
      {
         Vertex v = cliquelist[cind][index];

         assert( cand[v] == (int) cind );

         // increase counter of cliques in childnode
         if ( index == 0 )
            ++dualbound;

         assert ( SCIPvarGetLbLocal(probdata->vars[v]) < 0.5 && SCIPvarGetUbLocal(probdata->vars[v]) > 0.5 );

         // create branch
         SCIP_NODE* childnode;
         SCIP_CALL( SCIPcreateChild(scip, &childnode, 0.0, SCIPgetLocalTransEstimate(scip)) );
         SCIP_CALL( SCIPupdateNodeDualbound(scip, childnode, (double) dualbound) );
         assert( dualbound <= probdata->nlocalones + ncliques + probdata->objoffset);

         // fix first node to 1
         SCIP_CALL( SCIPchgVarLbNode(scip, childnode, probdata->vars[v], 1.0) );

         ++nbranch;
         *result = SCIP_BRANCHED;

         // fix remaining nodes to 0
         for ( size_t i2 = i; i2 < ncliques; ++i2 )
         {
            size_t cind2 = (size_t) cliqueorder[i2];
            assert( cind2 < ngen );

            // clique should not be marked
            assert ( ! marked[cind2] );

            // ignore fractional clique
            if ( cind2 == ssqi )
               continue;

            for ( size_t index2 = 0; index2 < cliquesizes[cind2]; ++index2)
            {
               // increment index2 if its lower than index
               if ( cind == cind2 && index2 <= index )
                  continue;

               Vertex w = cliquelist[cind2][index2];
               assert( cand[w] == (int) cind2 );

               assert ( SCIPvarGetLbLocal(probdata->vars[w]) < 0.5 && SCIPvarGetUbLocal(probdata->vars[w]) > 0.5 );

               // fix to 0
               SCIP_CALL( SCIPchgVarUbNode(scip, childnode, probdata->vars[w], 0.0) );
            }
         }

         // if no LP branching we are done here
         if ( ! islpbranching )
            continue;

         // fix frac node and its clique to 0
         assert( branchclique != nullptr );
         for (size_t k = 0; k < cliquesizes[ssqi]; ++k)
         {
            SCIP_CALL( SCIPchgVarUbNode(scip, childnode, probdata->vars[branchclique[k]], 0.0) );
         }
      }
   }


   for (size_t i = 0; i < ngen; ++i)
      SCIPfreeBlockMemoryArray(scip, &(cliquelist[i]), cliquesizes[i]);
   SCIPfreeBlockMemoryArray(scip, &cliquelist, ngen);

   if ( islpbranching )
   {
      // increase counter of cliques in childnode
      ++dualbound;

      // create additional branch fixing lastclique
      assert( branchclique != nullptr );
      for (std::size_t i = 0; i < cliquesizes[ssqi]; ++i)
      {
         SCIP_NODE* childnode;
         SCIP_CALL( SCIPcreateChild(scip, &childnode, 0.0, SCIPgetLocalTransEstimate(scip)) );
         SCIP_CALL( SCIPupdateNodeDualbound(scip, childnode, (double) dualbound) );
         assert( dualbound <= probdata->nlocalones + ncliques + probdata->objoffset);

         SCIP_CALL( SCIPchgVarLbNode(scip, childnode, probdata->vars[branchclique[i]], 1.0) );

         for (size_t j = i + 1; j < cliquesizes[ssqi]; ++j)
         {
            SCIP_CALL( SCIPchgVarUbNode(scip, childnode, probdata->vars[branchclique[j]], 0.0) );
         }
         ++nbranch;
         *result = SCIP_BRANCHED;
      }

      SCIPfreeBufferArray(scip, &branchclique);
   }

   SCIPdebugMsg(scip, "Clique partition branching created %lu branches (branched on %lu cliques).\n", nbranch, ncliques - lowerbound);
   SCIPfreeBufferArray(scip, &cliqueorder);
   SCIPfreeBufferArray(scip, &marked);

   return SCIP_OKAY;
}


/*
 * Callback methods of branching rule
 */

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpCliquePartition)
{
   assert( scip != nullptr );
   assert( branchrule != nullptr );
   assert( result != nullptr );

   SCIP_BRANCHRULEDATA* branchruledata = SCIPbranchruleGetData(branchrule);
   assert( branchruledata != nullptr );

   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   *result = SCIP_DIDNOTRUN;

   // do not run too early in the tree
   if ( SCIPgetDepth(scip) < branchruledata->mindepth )
      return SCIP_OKAY;

   // do not run on sparse graphs
   if ( SCIPisLE(scip, probdata->localdensity, branchruledata->mindensity) )
      return SCIP_OKAY;

   // only run when gap is not too large
   if ( SCIPisGT(scip, SCIPgetUpperbound(scip) - SCIPgetLowerbound(scip), branchruledata->maxabsgap) )
      return SCIP_OKAY;

   // only run on unweighted graphs
   if ( ! probdata->unweighted )
      return SCIP_OKAY;

   // run branching rule
   SCIPdebugMsg(scip, "Running clique partition branching rule ...\n");
   SCIP_CALL( branchPartition(scip, probdata, branchruledata, result, TRUE) );

   // update dual bound at current node
   SCIP_Real upperbound = probdata->nlocalones + probdata->ncliques + (size_t) probdata->objoffset;
   SCIP_CALL( SCIPupdateLocalDualbound(scip, upperbound) );

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsCliquePartition)
{
   assert( scip != nullptr );
   assert( branchrule != nullptr );
   assert( result != nullptr );

   SCIP_BRANCHRULEDATA* branchruledata = SCIPbranchruleGetData(branchrule);
   assert( branchruledata != nullptr );

   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   *result = SCIP_DIDNOTRUN;

   // do not run too early in the tree
   if ( SCIPgetDepth(scip) < branchruledata->mindepth )
      return SCIP_OKAY;

   // do not run on sparse graphs
   if ( SCIPisLE(scip, probdata->localdensity, branchruledata->mindensity) )
      return SCIP_OKAY;

   // only run when gap is not too large
   if ( SCIPisGT(scip, SCIPgetUpperbound(scip) - SCIPgetLowerbound(scip), branchruledata->maxabsgap) )
      return SCIP_OKAY;

   // only run on unweighted graphs
   if ( ! probdata->unweighted )
      return SCIP_OKAY;

   // run branching rule
   SCIPdebugMsg(scip, "Running clique partition branching rule ...\n");
   SCIP_CALL( branchPartition(scip, probdata, branchruledata, result, FALSE) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeCliquePartition)
{
   assert( scip != nullptr );
   assert( branchrule != nullptr );

   SCIP_BRANCHRULEDATA* branchruledata = SCIPbranchruleGetData(branchrule);

   SCIPfreeBlockMemoryNull(scip, &branchruledata);

   return SCIP_OKAY;
}


/*
 * branching rule specific interface methods
 */

/** creates the clique partition branching rule and includes it in SCIP */
SCIP_RETCODE BACSincludeBranchruleCliquePartition(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != nullptr );

   SCIP_BRANCHRULEDATA* branchruledata = nullptr;
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );

   SCIP_BRANCHRULE* branchrule = nullptr;

   /* create deg branching rule data */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert( branchrule != nullptr );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpCliquePartition) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsCliquePartition) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeCliquePartition) );

   SCIP_CALL( SCIPaddIntParam(scip, "branching/" BRANCHRULE_NAME "/mindepth",
         "minimum depth of branching rule",
         &branchruledata->mindepth, FALSE, DEFAULT_MINDEPTH, 0, INT_MAX, nullptr, nullptr) );

   SCIP_CALL( SCIPaddRealParam(scip, "branching/" BRANCHRULE_NAME "/mindensity",
         "minimum density such that branching rule is applied",
         &branchruledata->mindensity, FALSE, DEFAULT_MINDENSITY, 0, INT_MAX, nullptr, nullptr) );

   SCIP_CALL( SCIPaddIntParam(scip, "branching/" BRANCHRULE_NAME "/maxabsgap",
         "maximum absolute solution gap",
         &branchruledata->maxabsgap, FALSE, DEFAULT_MAXABSGAP, 0, INT_MAX, nullptr, nullptr) );

   SCIP_CALL( SCIPaddIntParam(scip, "branching/" BRANCHRULE_NAME "/fracnodechoice",
         "method of choosing fracnode",
         &branchruledata->fracnodechoice, FALSE, DEFAULT_FRACNODECHOICE, 0, INT_MAX, nullptr, nullptr) );

   return SCIP_OKAY;
}
