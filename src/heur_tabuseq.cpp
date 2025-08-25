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

/**@file   heur_tabuseq.cpp
 * @ingroup DEFPLUGINS_HEUR
 * @brief  sequential improvement tabu search primal heuristic
 * @author Marc Pfetsch
 *
 * We take an existing solution and try to find a stable set that is one larger. During the loop, we maintain a solution
 * that is not stable, i.e., it has a certain number of violated edges. If this number reaches 0, we have found a stable
 * set and increase the size. Each move consists of picking a non-tabu node x in the solution and non-tabu node y that
 * is not in the solution. Then the best pair produces a new solution by exchanging their roles w.r.t. the solution.
 *
 * The principle move structure was introduced in@n
 * STABULUS: A technique for finding stable sets in large graphs with tabu search@n
 * C. Friden, A. Hertz and D. de Werra@n
 * Computing 42, 35–44 (1989), DOI: 10.1007/BF02243141
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <set>

#include "heur_tabuseq.h"
#include "graph.h"
#include "struct_probdata.h"
#include <scip/scip_sol.h>
#include <scip/scip_randnumgen.h>
#include <scip/sol.h>

#define HEUR_NAME            "tabuseq"
#define HEUR_DESC            "sequential improvement tabu search for stable set solutions"
#define HEUR_DISPCHAR        'T'
#define HEUR_PRIORITY        5000
#define HEUR_FREQ            1
#define HEUR_FREQOFS         0
#define HEUR_MAXDEPTH        -1
#define HEUR_TIMING          SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP     FALSE /**< does the heuristic use a secondary SCIP instance? */

// defines for heuristic
#define TABUSEQ_maxIter          30000     //!< maximal number of iterations
#define TABUSEQ_fixed            100000    //!< fixed nodes get huge iter number, assert( TABUSEQ_fixed > TABUSEQ_maxIter )
#ifdef SCIP_DEBUG
#define TABUSEQ_freq             10000     //!< frequency for output
#endif


/** primal heuristic data */
struct SCIP_HeurData
{
   int                   lastsolindex;       //!< value of last solution tried to improve
   SCIP_RANDNUMGEN*      randnumgen;         //!< randon number generator
};


//! determine the number of conflicts of the solution @a S
size_t numberConcflicts(
   const Graph*          G,                  /**< graph */
   const SCIP_Bool*      S                   /**< solution */
   )
{
   size_t nconflicts = 0;
   EdgeIterator eit, eend;
   for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
   {
      Vertex s = boost::source(*eit, *G);
      Vertex t = boost::target(*eit, *G);

      if ( S[s] && S[t] )
         ++nconflicts;
   }
   return nconflicts;
}

#ifndef NDEBUG
//! check whether adjacency data is correct
SCIP_Bool checkAData(
   const Graph*          G,                  /**< graph */
   const SCIP_Bool*      S,                  /**< solution */
   const size_t*         A                   /**< number of neighbors in solution */
   )
{
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      size_t nadj = 0;
      AdjacencyIterator ait, aend;
      for (boost::tie(ait,aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         if ( S[*ait] )
            ++nadj;
      }
      if ( nadj != A[v] )
         return FALSE;
   }
   return TRUE;
}

//! check whether solution data is correct
SCIP_Bool checkSData(
   const Graph*          G,                  /**< graph */
   const SCIP_Bool*      S,                  /**< solution bitset */
   const Vertex*         S_arr,              /**< solution array */
   size_t                size                /**< solution size */
   )
{
   for (size_t i = 0; i < size; ++i)
   {
      Vertex v = S_arr[i];
      if ( ! S[v] )
         return FALSE;
   }

   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      if ( ! S[v] )
         continue;
      SCIP_Bool found = FALSE;
      for (size_t i = 0; i < size; ++i)
      {
         if ( v == S_arr[i] )
         {
            found = TRUE;
            break;
         }
      }
      if ( ! found )
         return FALSE;
   }
   return TRUE;
}

//! check size
size_t checkSize(
   size_t                n,                  /**< number of nodes */
   const SCIP_Bool*      S                   /**< solution */
   )
{
   size_t size = 0;
   for (size_t i = 0; i < n; ++i)
   {
      if ( S[i] )
         ++size;
   }
   return size;
}
#endif

/** calculation of tabu lengths
 *
 *  The choice is explained in
 *  "An adaptive multistart tabu search approach to solve the maximum clique problem" - Qinghua Wu, Jin-Kao Hao
 */
static
void calcTabuLength(
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   size_t                nconflicts,         /**< current number of conflicts in S */
   size_t&               Tin,                /**< number of tabu iterations for leaving node */
   size_t&               Tout                /**< number of tabu iterations for entering node */
   )
{
   assert( heurdata != nullptr );
   assert( heurdata->randnumgen != nullptr );

   Tin = 10 + (size_t) (nconflicts + (size_t) SCIPrandomGetInt(heurdata->randnumgen, 1, 10 + (int) nconflicts));
   Tout = 10 + (size_t) (0.6 * nconflicts + (size_t) SCIPrandomGetInt(heurdata->randnumgen, 1, 10 + (int) nconflicts));
}


/** random variable deciding whether a random swap should be executed
 *
 *  The choice is explained in
 *  "An adaptive multistart tabu search approach to solve the maximum clique problem" - Qinghua Wu, Jin-Kao Hao
 */
static
SCIP_Bool shouldDoRandomMove(
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   size_t                nnodes,             /**< size of graph */
   size_t                size,               /**< size of current S */
   size_t                nconflicts          /**< current number of conflicts in S */
   )
{
   int l = (int) (size * (size - 1) / 2 - nconflicts);

   SCIP_Real p = MIN( (l + 2.0) / (SCIP_Real) nnodes, 0.1);
   if ( SCIPrandomGetReal(heurdata->randnumgen, 0.0, 1.0) < p )
      return TRUE;
   return FALSE;
}

//! run tabuseq seach algorithm to improve upon any solution
static
SCIP_RETCODE tabuseqSearch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   SCIP_HEUR*            heur,               /**< current heuristic */
   SCIP_SOL*             bestsol,            /**< initial solution to start tabuseq search with */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( heur != nullptr );
   assert( bestsol != nullptr );
   assert( result != nullptr );

   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);

   size_t n = probdata->n;
   const Graph* G = probdata->G;
   assert( G != nullptr );
   size_t bestobj = (size_t) SCIPgetSolOrigObj(scip, bestsol);

   // allocate data
   size_t* T;        // iteration at which tabuseq node will expire to be tabuseq
   size_t* A;           // number of adjacent nodes in solution
   SCIP_Bool* S;     // bitset of solution
   Vertex* S_arr;    // array of solution
   Vertex* marked;   // mark neighbors of nodein
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &T, n) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &A, n) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &S, n) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &S_arr, n) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &marked, n) );

   // initialize data
   // - node globally fixed to 0: set T value to TABUSEQ_fixed (will never be added to solution)
   // - node globally fixed to 1: add to solutions and fix all neighbors to zero
   size_t solsize = 0;
   size_t maxsize = n;  // number nodes that are not fixed to 0
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      assert( SCIPisEQ(scip, boost::get(vertex_weight_t(), *G, v), 1.0) ); // only works for unweighted instances

      SCIP_VAR* var = probdata->vars[v];
      if ( SCIPvarGetUbGlobal(var) < 0.5 && T[v] < TABUSEQ_fixed )
      {
         T[v] = TABUSEQ_fixed;
         --maxsize;
      }
      else if ( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         T[v] = TABUSEQ_fixed;
         S[v] = TRUE;
         S_arr[solsize] = v;
         ++solsize;

         // mark neighbors of variables that are fixed to 1
         AdjacencyIterator ait, aend;
         for (boost::tie(ait,aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            if ( T[*ait] < TABUSEQ_fixed )
            {
               T[*ait] = TABUSEQ_fixed;
               --maxsize;
            }
         }
      }
   }

   // initialize to given solution
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      SCIP_VAR* var = probdata->vars[v];

      if ( SCIPgetSolVal(scip, bestsol, var) > 0.5 )
      {
         // skip nodes already in the solution
         if ( S[v] )
            continue;

         // skip nodes that are globally fixed to 0.0
         if ( T[v] > TABUSEQ_maxIter )
            continue;

         assert( ! S[v] );
         assert( T[v] == 0 );

         S[v] = TRUE;
         S_arr[solsize] = v;
         ++solsize;
      }
   }
   assert( checkSData(G, S, S_arr, solsize) );

   // initialize A data
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;

      AdjacencyIterator ait, aend;
      for (boost::tie(ait,aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         if ( S[*ait] )
            ++A[v];
      }
   }
   assert( checkAData(G, S, A) );

   // loop for iteratively increasing size of solution
   for (size_t size = bestobj + 1; size < maxsize; ++size)
   {
      // fill S until solsize == size
      while ( solsize < size )
      {
         // search for node v not in solution with minimal A[v]
         size_t s = Graph::null_vertex();
         size_t bestA = INT_MAX;

         for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
         {
            Vertex v = *vit;

            if ( S[v] || T[v] > TABUSEQ_maxIter)
               continue;

            if ( A[v] < bestA )
            {
               s = v;
               bestA = A[v];
               if ( bestA == 0 )
                  break;
            }
         }
         // break if sol of size bestobj + 1 could not be built with local fixings
         if ( bestA == INT_MAX )
            break;

         assert( s < n );

         // include node into current solution
         assert( ! S[s] );
         S[s] = TRUE;
         S_arr[solsize] = s;
         ++solsize;
         assert( solsize == checkSize(n, S) );

         // update A
         AdjacencyIterator ait, aend;
         for (boost::tie(ait,aend) = boost::adjacent_vertices(s, *G); ait != aend; ++ait)
            ++A[*ait];

         assert( checkAData(G, S, A) );
      }
      assert( solsize == size );
      assert( checkSData(G, S, S_arr, solsize) );

      // initialize tabu list
      for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
      {
         Vertex v = *vit;
         if ( T[v] < TABUSEQ_fixed )
            T[v] = 0;
      }

      // initialize the number of conflicts (= number of violated edges)
      size_t nconflicts = numberConcflicts(G, S);
      size_t bestnconflicts = nconflicts;

      // main tabu loop with fixed number of iterations
      for (size_t iter = 1; iter <= TABUSEQ_maxIter; ++iter)
      {
         // exit loop if we have found a feasible solution
         if ( nconflicts == 0 )
            break;

         // search for best pair of nodes to be exchanged in solution
         int bestdelta = INT_MAX;
         Vertex bestnodein = Graph::null_vertex();
         Vertex bestnodeout = Graph::null_vertex();
         size_t besti = INT_MAX;
         size_t maxAinS = 0;

         // determine maximum number of conflicts in S
         for (size_t i = 0; i < size; ++i)
         {
            if ( T[S_arr[i]] > iter )
               continue;

            if ( A[S_arr[i]] > maxAinS )
               maxAinS = A[S_arr[i]];
         }

         // initialize marked
         for (size_t i = 0; i < n; ++i)
            marked[i] = Graph::null_vertex();

         // find swap, we choose bestnodein with maximum number of conflicts
         for (size_t i = 0; i < size; ++i)
         {
            Vertex nodein = S_arr[i];
            assert( S[nodein] );

            // choose nodein as bestnodein if not fixed
            if ( T[nodein] == TABUSEQ_fixed )
               continue;

            // only accept maximal number of conflicts
            if ( A[nodein] < maxAinS )
               continue;

            assert( nodein != Graph::null_vertex() );
            assert( nodein < n );
            assert( A[nodein] == maxAinS || T[nodein] > iter );

            SCIP_Bool aspiration = FALSE;

            // first search in neighborhood
            AdjacencyIterator ait, aend;
            for (boost::tie(ait, aend) = boost::adjacent_vertices(nodein, *G); ait != aend && ! aspiration; ++ait)
            {
               Vertex nodeout = *ait;

               // consider nodes not in the solution and not tabu
               if ( S[nodeout] || T[nodeout] > iter )
                  continue;

               marked[nodeout] = nodein;

               // compute change in nconflicts, for adjacent nodes A[nodeout] is off by one:
               int delta = (int) A[nodeout] - 1 - (int) A[nodein];

               // we accept if we number of conflicts decreases or nodein is not tabu and currently best swap found
               if ( delta < 0 || ( T[nodein] < iter && delta < bestdelta ) )
               {
                  bestdelta = delta;
                  bestnodeout = nodeout;
                  bestnodein = nodein;
                  besti = i;

                  if ( delta < 0 )
                  {
                     aspiration = TRUE;
                     break;
                  }
               }
            }

            VertexIterator vit2, vend2;
            for (boost::tie(vit2, vend2) = boost::vertices(*G); vit2 != vend2 && ! aspiration; ++vit2)
            {
               Vertex nodeout = *vit2;

               // neighbors of nodein have already been covered
               if ( marked[nodeout] == nodein )
                  continue;

               // consider nodes not in the solution and not tabu
               if ( S[nodeout] || T[nodeout] > iter )
                  continue;

               // compute change in nconflicts:
               int delta = (int) A[nodeout] - (int) A[nodein];

               // we accept if we number of conflicts decreases or nodein is not tabu and currently best swap found
               if ( delta < 0 || ( T[nodein] < iter && delta < bestdelta ) )
               {
                  bestdelta = delta;
                  bestnodeout = nodeout;
                  bestnodein = nodein;
                  besti = i;

                  if ( delta < 0 )
                  {
                     aspiration = TRUE;
                     break;
                  }
               }
            }
         }

         // do a random swap if no swap has been found yet or with a certain probabilty when no progress was made
         if ( bestdelta == INT_MAX || ( bestdelta >= 0 && shouldDoRandomMove(heurdata, n, size, nconflicts)) )
         {
            // choose nodein randomly by respecting global fixings
            size_t starti = (size_t) SCIPrandomGetInt(heurdata->randnumgen, 0, (int) size - 1);
            besti = starti;
            bestnodein = S_arr[besti];
            while ( T[S_arr[besti]] >= TABUSEQ_fixed )
            {
               ++besti;
               if ( besti == size )
                  besti = 0;
               if ( besti == starti )
               {
                  bestnodein = Graph::null_vertex();
                  break;
               }
               bestnodein = S_arr[besti]; /*lint !e838*/
            }

            assert( S[bestnodein] || bestnodein == Graph::null_vertex() );

            // choose nodeout sufficiently bad by respecting global fixings
            size_t startnodeout = (size_t) SCIPrandomGetInt(heurdata->randnumgen, 0, (int) n - 1);
            bestnodeout = startnodeout;
            while ( S[bestnodeout] || T[bestnodeout] >= TABUSEQ_fixed || A[bestnodeout] < (size_t) (size * probdata->m / (n * n)))
            {
               ++bestnodeout;
               if ( bestnodeout == n )
                  bestnodeout = 0;
               if ( bestnodeout == startnodeout )
               {
                  bestnodeout = Graph::null_vertex();
                  break;
               }
            }

            if ( bestnodeout != Graph::null_vertex() )
            {
               assert( ! S[bestnodeout] );
               assert( T[bestnodeout] < TABUSEQ_fixed );
               assert( bestnodeout < n );
            }
         }

         // no exchange was found - may happen for very small graphs where everything is tabu
         if ( bestnodein == Graph::null_vertex() || bestnodeout == Graph::null_vertex() )
            break;

         assert( bestnodein != Graph::null_vertex() );
         assert( bestnodeout != Graph::null_vertex() );
         assert( S[bestnodein] );
         assert( ! S[bestnodeout] );
         assert( S_arr[besti] == bestnodein );
         assert( besti < size );

         // change solution
         S[bestnodein] = FALSE;
         S[bestnodeout] = TRUE;
         S_arr[besti] = bestnodeout;

         // update tabu list
         size_t Tin;
         size_t Tout;
         calcTabuLength(heurdata, nconflicts, Tin, Tout);

         T[bestnodein] = iter + Tin;
         T[bestnodeout] = iter + Tout;

         // update A and check whether the nodes are neighbors
         size_t neighbors = 0;
         AdjacencyIterator ait, aend;
         for (boost::tie(ait,aend) = boost::adjacent_vertices(bestnodeout, *G); ait != aend; ++ait)
         {
            ++A[*ait];
            if ( *ait == bestnodein )
               ++neighbors;
         }
         for (boost::tie(ait,aend) = boost::adjacent_vertices(bestnodein, *G); ait != aend; ++ait)
            --A[*ait];

         // new number of conflicts
         nconflicts = nconflicts + A[bestnodeout] - A[bestnodein] + neighbors;

         assert( nconflicts == numberConcflicts(G, S) );
         assert( checkAData(G, S, A) );
         assert( checkSData(G, S, S_arr, size) );

#ifdef SCIP_DEBUG
         // possibly output
         if ( iter == 1 || iter % TABUSEQ_freq == 0 || nconflicts < bestnconflicts )
         {
            SCIPinfoMessage(scip, nullptr, "size: %lu\titer: %7lu\tconf: %4lu\tbest: %4lu\tnode in: %4lu\t\tnode out: %4lu\tdelta: %d\n",
               size, iter, nconflicts, bestnconflicts, bestnodein, bestnodeout, bestdelta);
         }
#endif

         if ( nconflicts < bestnconflicts )
            bestnconflicts = nconflicts;
      }

      // if we found a feasible solution
      if ( nconflicts == 0 )
      {
         // pass solution to SCIP
         SCIP_SOL* sol;
         SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
         for (size_t i = 0; i < size; ++i)
         {
            SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->vars[S_arr[i]], 1.0) );
         }
         assert( SCIPisEQ(scip, (SCIP_Real) size + SCIPgetOrigObjoffset(scip), SCIPgetSolOrigObj(scip, sol)) );

         SCIP_Bool stored;
         SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, TRUE, FALSE, &stored) );

         if ( stored )
         {
            *result = SCIP_FOUNDSOL;
            SCIPdebugMsg(scip, "Tabuseq search found solution of size %lu.\n", size);
         }
         else
            SCIPdebugMsg(scip, "Tabuseq search did not find new solution.\n");

#ifdef SCIP_DEBUG
         SCIPinfoMessage(scip, nullptr, "Found solution of size %lu. Increase size ...\n", size);
#endif
      }
      else
         break;
   }
   SCIPfreeBlockMemoryArray(scip, &marked, n);
   SCIPfreeBlockMemoryArray(scip, &S_arr, n);
   SCIPfreeBlockMemoryArray(scip, &S, n);
   SCIPfreeBlockMemoryArray(scip, &A, n);
   SCIPfreeBlockMemoryArray(scip, &T, n);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTabuseq)
{
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( result != nullptr );
   static_assert(TABUSEQ_fixed > TABUSEQ_maxIter, "TABUSEQ_fixed is not large enough.");

   *result = SCIP_DIDNOTRUN;

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   size_t n = probdata->n;
   if ( n == 0 )
      return SCIP_OKAY;

   // the heuristic currently only works for unweighted graphs
   if ( ! probdata->unweighted )
      return SCIP_OKAY;

   SCIP_SOL* bestsol = SCIPgetBestSol(scip);
   if ( bestsol == nullptr )
      return SCIP_OKAY;

   SCIP_HeurData* heurdata = SCIPheurGetData(heur);
   int index = SCIPsolGetIndex(bestsol);

   if ( index == heurdata->lastsolindex )
      return SCIP_OKAY;
   heurdata->lastsolindex = index;

   SCIPdebugMsg(scip, "Running tabuseq search heuristic ...\n");
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( tabuseqSearch(scip, probdata, heur, bestsol, result) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeTabuseq)
{
   assert( scip != nullptr );
   assert( heur != nullptr );

   SCIP_HEURDATA* heurdata;
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != nullptr );

   SCIPfreeRandom(scip, &(heurdata->randnumgen));
   SCIPfreeBlockMemoryNull(scip, &heurdata);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates tabuseq heuristic and includes it in SCIP */
SCIP_RETCODE BACSincludeHeurTabuseq(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata = nullptr;
   SCIP_HEUR* heur = nullptr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   heurdata->lastsolindex = -1;
   SCIP_CALL( SCIPcreateRandom(scip, &(heurdata->randnumgen), 0, TRUE) );

   SCIP_CALL(SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecTabuseq, heurdata));

   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeTabuseq) );

   assert( heur != nullptr );

   return SCIP_OKAY;
}
