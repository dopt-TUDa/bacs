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

/**@file   cliquepartition.cpp
 * @brief  functions to compute cliquepartition on local graph
 * @author Jonas Alker
 */

#include "cliquepartition.h"
#include "struct_probdata.h"


#ifndef NDEBUG
void BACScheckCliquesizes(
   SCIP*                 scip,               /**< SCIP pointer */
   int*                  cand,               /**< clique partition, representant for every vertex */
   size_t*               cliquesizes,        /**< sizes of cliques */
   size_t                ngen,               /**< number of initially generated cliques */
   size_t                n                   /**< number of vertices */
   )
{
   assert( cand!= nullptr );
   assert( cliquesizes != nullptr );
   assert( ngen < n );

   size_t* cliquecounter;
   (void) SCIPallocClearBufferArray(scip, &cliquecounter, (int) ngen);
   for (size_t i = 0; i < n; ++i)
   {
      if ( cand[i] < 0 )
         continue;

      assert( cand[i] < (int) ngen );

      ++(cliquecounter[cand[i]]);
   }

   for (size_t i = 0; i < ngen; ++i)
   {
      assert( cliquecounter[i] == cliquesizes[i] );
   }
   SCIPfreeBufferArray(scip, &cliquecounter);
}
#endif


/** compute lowerbound on current locally valid graph respecting the current primalbound, the objective offset and locally fixed variables */
SCIP_RETCODE BACScomputeLocalLowerbound(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   size_t&               lowerbound          /**< lowerbound regarding free nodes */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( probdata->unweighted );

   // compute local lower bound by subtracting number of fixed variables from best primal bound
   SCIP_Real objval = SCIPgetPrimalbound(scip);
   assert( SCIPisFeasIntegral(scip, objval) );
   assert( SCIPisFeasGE(scip, objval, 0.0) );

   // compute local lower bound
   if ( SCIPisGT(scip, (double) probdata->nlocalones + probdata->objoffset, objval) )
      lowerbound = 0;   // we already found an improvemt with fixings if we are feasible
   else
      lowerbound = (size_t) SCIPfeasRound(scip, objval - probdata->objoffset - probdata->nlocalones);

   return SCIP_OKAY;
}


/** try to improve given partition by removing cliques from the partition and moving contained vertices to other cliques
 *
 *  if successful, the method runs again until no further improvement can be found or lowerbound cliques have been reached
 */
SCIP_RETCODE improvePartition(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int*                  cand,               /**< clique partition, representant for every vertex */
   size_t*               cliquesizes,        /**< sizes of cliques */
   size_t                ngen,               /**< number of initially generated cliques */
   size_t                lowerbound,         /**< lowerbound regarding free nodes */
   size_t&               ncliques,           /**< current number of cliques */
   SCIP_Bool&            cutoff              /**< whether a cutoff was detected */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( probdata->unweighted );
   assert( cliquesizes != nullptr );

   const Graph* G = probdata->G;
   size_t n = probdata->n;

   size_t* nswapablenodes;       //!< counter how many nodes per clique can be moved to another clique
   size_t* targetclique;         //!< index of clique a node can we moved to
   size_t* nnodesinclique;       //!< counter for a node v how many nodes per clique are covered by v's neighborhood
   SCIP_CALL( SCIPallocBufferArray(scip, &nswapablenodes, ngen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &targetclique, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nnodesinclique, ngen) );

   for ( size_t ncurrcliques = ncliques; ncurrcliques > lowerbound; --ncurrcliques )
   {
      // check if ngen can be reduced by one:
      SCIP_Bool reducedcliques = FALSE;

      BMSclearMemoryArray(nswapablenodes, ngen);

      VertexIterator vit, vend;
      for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
      {
         Vertex v = *vit;

         if ( cand[v] < 0 )
            continue;

         assert( SCIPvarGetUbLocal(probdata->vars[v]) > 0.5 && SCIPvarGetLbLocal(probdata->vars[v]) < 0.5 );

         // check whether v can be added in other clique than itself:
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

            // if clique is covered by neighborhood of v
            if ( nnodesinclique[cand[w]] >= cliquesizes[cand[w]] )
            {
               ++(nswapablenodes[cand[v]]);
               targetclique[v] = (size_t) cand[w];

               assert( targetclique[v] < ngen );
               assert( cand[v] >= 0 );
               assert( cand[v] != cand[w] );
               assert( cliquesizes[cand[w]] > 0 );

               // check if clique cand[v] can be removed
               if ( nswapablenodes[cand[v]] == cliquesizes[cand[v]] )
               {
                  int emptyclique = cand[v];
                  reducedcliques = TRUE;

                  if ( ncurrcliques == lowerbound + 1 )
                  {
                     cutoff = TRUE;
                     //printf("Node at depth: %d, cut off due to improvement, ngen: %lu, lowerbound: %lu\n", SCIPgetDepth(scip), ngen, lowerbound);
                     SCIPdebugMsg(scip, "Node at depth: %d, cut off due to cliquecover, ngen: %lu, lowerbound: %lu\n", SCIPgetDepth(scip), ngen, lowerbound);

                     SCIPfreeBufferArray(scip, &nnodesinclique);
                     SCIPfreeBufferArray(scip, &targetclique);
                     SCIPfreeBufferArray(scip, &nswapablenodes);

                     return SCIP_OKAY;
                  }

                  // reallocate all nodes in clique cand[v]:
                  if ( cliquesizes[cand[v]] == 1 )
                  {
                     assert( nswapablenodes[cand[v]] == 1 );

                     ++(cliquesizes[targetclique[v]]);
                     --(cliquesizes[cand[v]]);
                     cand[v] = (int) targetclique[v];
#ifndef NDEBUG
                     BACScheckCliquesizes(scip, cand, cliquesizes, ngen, n);
#endif
                  }
                  else
                  {
                     assert( nswapablenodes[cand[v]] == cliquesizes[cand[v]] );

                     for (size_t j = 0; j < n; ++j)
                     {
                        assert( nswapablenodes[emptyclique] >= cliquesizes[cand[v]] );

                        if ( cand[j] != emptyclique )
                           continue;

                        assert( cand[j] >= 0 );
                        assert( targetclique[j] < ngen );
                        assert( (int) targetclique[j] != cand[j] );

                        ++(cliquesizes[targetclique[j]]);
                        --(cliquesizes[cand[j]]);
                        cand[j] = (int) targetclique[j];
#ifndef NDEBUG
                        BACScheckCliquesizes(scip, cand, cliquesizes, ngen, n);
#endif
                        if ( cliquesizes[emptyclique] == 0 )
                           break;
                     }
                  }
               }
               // we have found a targetclique
               break;
            }
         }

         if ( reducedcliques )
         {
            --ncliques;
            break;
         }
      }

      if ( ! reducedcliques )
         break;
   }
   SCIPfreeBufferArray(scip, &nnodesinclique);
   SCIPfreeBufferArray(scip, &targetclique);
   SCIPfreeBufferArray(scip, &nswapablenodes);

   return SCIP_OKAY;
}


/** compute clique partition in O(n + m)
 *
 *  this partitioning algorithm is inspired by a simple coloring heuristic
 *
 *  we start with an empty partitioning
 *
 *  for all vertices v in G we check whether v can be added to an existing clique in the partitioning by scanning the neighborhood of v
 *  and checking whether a clique in the current partitioning is covered by the neighborhood of v, then we add v oder create a new clique
 *
 *  we scan every vertex exactly once and all adjacent edges at most once
 *
 *  @todo running time is not linear currently, as we reset nnodesinclique counter, we iterate over adjacency of v again to obtain truely linear running time
 */
SCIP_RETCODE computeCliquePartitionColoring(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int*                  cand,               /**< clique partition, representant for every vertex */
   size_t*               cliquesizes,        /**< sizes of cliques */
   size_t&               ngen                /**< number of initially generated cliques */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( cand != nullptr );
   assert( cliquesizes != nullptr );

   ngen = 0;
   const Graph* G = probdata->G;

   size_t* nnodesinclique;
   SCIP_CALL( SCIPallocBufferArray(scip, &nnodesinclique, probdata->n) );

   Vertex* vertices;
   int* degs;
   SCIP_CALL( SCIPallocBufferArray(scip, &vertices, probdata->n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &degs, probdata->n) );

   // initialize
   for (size_t i = 0; i < probdata->n; ++i)
   {
      vertices[i] = i;
      degs[i] = probdata->localdegrees[i];
      cliquesizes[i] = 0;
      cand[i] = -1;
   }

   // sort vertices
   SCIPsortIntPtr(degs, (void**) vertices, (int) probdata->n);
   SCIPfreeBufferArray(scip, &degs);

   for ( size_t i = 0; i < probdata->n; ++i )
   {
      Vertex v = vertices[i];

      // skip nodes already covered
      if ( SCIPvarGetUbLocal(probdata->vars[v]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[v]) > 0.5 )
         continue;

      assert( cand[v] == -1 );

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
      }
   }
   SCIPfreeBufferArray(scip, &vertices);
   SCIPfreeBufferArray(scip, &nnodesinclique);

   return SCIP_OKAY;
}


/** compute clique partition in O(m)
 *
 *  start with a trivial partition of the graph G where all nodes are in the same subset
 *
 *  for all vertices v in G we split the subset containing v such that all non-adjacent vertices of v
 *  are in a different subset that v itself
 *
 *  after iterating over all vertices every vertex is connected with every other vertex in the same subset - a clique
 *
 *  @todo running time is not linear currently, as we reset nnodesinclique counter, we iterate over adjacency of v again to obtain truely linear running time
 */
SCIP_RETCODE computeCliquePartitionLinear(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int*                  cand,               /**< clique partition, representant for every vertex */
   size_t*               cliquesizes,        /**< sizes of cliques */
   size_t&               ngen                /**< number of initially generated cliques */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( cand != nullptr );
   assert( cliquesizes != nullptr );

   ngen = 0;
   const Graph* G = probdata->G;

   size_t* nnodesinclique;
   SCIP_CALL( SCIPallocBufferArray(scip, &nnodesinclique, probdata->n) );

   Vertex* vertices;
   int* degs;
   SCIP_CALL( SCIPallocBufferArray(scip, &vertices, probdata->n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &degs, probdata->n) );
   size_t nfixed = 0;

   // initialize
   for (size_t i = 0; i < probdata->n; ++i)
   {
      vertices[i] = i;
      degs[i] = probdata->localdegrees[i];

      cliquesizes[i] = 0;
      if ( SCIPvarGetUbLocal(probdata->vars[i]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[i]) > 0.5 )
      {
         cand[i] = -1;
         ++nfixed;
      }
      else
      {
         cand[i] = 0;
         ++(cliquesizes[0]);
      }
   }
   ngen = 1;

   // sort vertices
   SCIPsortIntPtr(degs, (void**) vertices, (int) probdata->n);
   SCIPfreeBufferArray(scip, &degs);

   int nextemptyclique = 1;
   for ( size_t i = nfixed; i < probdata->n; ++i )
   {
      Vertex v = vertices[i];

      assert( cand[v] >= 0 );
      assert( cand[v] < (int) ngen + 1);

      // check whether v can be moved to existing nonempty clique
      BMSclearMemoryArray(nnodesinclique, ngen + 1); // if nextemptyclique != ngen all cliques have labels in [0, ngen]

      SCIP_Bool foundclique = FALSE;
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         if ( cand[w] == -1 )
            continue;
         assert( cand[w] >= 0 );
         assert( cand[w] < (int) ngen + 1);

         ++(nnodesinclique[cand[w]]);

         if ( nnodesinclique[cand[w]] == cliquesizes[cand[w]] )
         {
            --(cliquesizes[cand[v]]);
            cand[v] = cand[w];
            ++(cliquesizes[cand[w]]);
            foundclique = TRUE;
            break;
         }
      }
      if ( foundclique )
         continue;

      // move v and all of its neighbors from the same clique in a new clique
      int oldclique = cand[v];
      int newclique = nextemptyclique;

      assert( cliquesizes[newclique] == 0 );
      assert( cliquesizes[oldclique] > 0 );

      cand[v] = newclique;
      --(cliquesizes[oldclique]);
      ++(cliquesizes[newclique]);

      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         // skip fixed nodes and nodes that were not in the same clique as v
         if ( cand[w] != oldclique )
            continue;

         assert( cand[w] >= 0 );
         assert( cand[w] == oldclique );
         assert( cliquesizes[oldclique] > 0 );
         assert( cliquesizes[newclique] > 0 );

         cand[w] = newclique;
         --(cliquesizes[oldclique]);
         ++(cliquesizes[newclique]);
      }
      assert( cliquesizes[newclique] > 0 );

      // check whether a new clique was created
      if( cliquesizes[oldclique] == 0 )
      {
         nextemptyclique = oldclique;
      }
      else
      {
         // new clique has been created by splitting oldclique, all labels in [0, ngen - 1] should be used
         // next clique gets label ngen
         ++ngen;
         nextemptyclique = (int) ngen;
      }
   }

   // possibly move all nodes with cand == ngen to nextemptyclique
   if ( nextemptyclique < (int) ngen )
   {
      for ( size_t i = nfixed; i < probdata->n; ++i )
      {
         Vertex v = vertices[i];

         if ( cand[v] == (int) ngen )
         {
            cand[v] = nextemptyclique;
            --(cliquesizes[ngen]);
            ++(cliquesizes[nextemptyclique]);

            if ( cliquesizes[ngen] == 0 )
               break;
         }
      }
   }
   SCIPfreeBufferArray(scip, &vertices);
   SCIPfreeBufferArray(scip, &nnodesinclique);

   return SCIP_OKAY;
}


//! compute greedy clique partition in O(n*deg^2)
SCIP_RETCODE computeCliquePartitionGreedy(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int*                  cand,               /**< clique partition, representant for every vertex */
   size_t*               cliquesizes,        /**< sizes of cliques */
   size_t&               ngen                /**< number of initially generated cliques */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( cand != nullptr );
   assert( cliquesizes != nullptr );

   ngen = 0;
   const Graph* G = probdata->G;

   // initialize
   for (size_t i = 0; i < probdata->n; ++i)
   {
      cand[i] = -1;
      cliquesizes[i] = 0;
   }

   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;

      // skip nodes already covered
      if ( SCIPvarGetUbLocal(probdata->vars[v]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[v]) > 0.5 )
         continue;

      if ( cand[v] >= 0 )
         continue;

      cand[v] = (int) ngen;
      size_t csize = 1;

      // sorted adjacency list
      size_t nneighbors = 0;
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         // skip fixed nodes
         if ( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[w]) > 0.5 )
            continue;

         // count seen (unfixed) neighbors
         if ( (int) nneighbors >= probdata->localdegrees[v] )
            break;
         ++nneighbors;

         // skip nodes already covered
         if ( cand[w] >= 0 )
            continue;

         // check whether w can be added to clique
         if ( csize == 1 )
         {
            ++csize;
            cand[w] = (int) ngen;
            continue;
         }

         unsigned int nneigh = 0;
         AdjacencyIterator bit, bend;
         for (boost::tie(bit, bend) = boost::adjacent_vertices(w, *G); bit != bend; ++bit)
         {
            Vertex x = *bit;
            if ( cand[x] == (int) ngen )
            {
               ++nneigh;
               if ( nneigh == csize )
                  break;
            }
         }

         // if the neighbors of w cover all nodes in the current clique, we can extend it
         if ( nneigh == csize )
         {
            ++csize;
            cand[w] = (int) ngen;
         }
      }
      cliquesizes[ngen] = csize;
      ++ngen;
   }

   return SCIP_OKAY;
}
