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

/**@file   heur_dynamicdeg.cpp
 * @brief  primal heuristic based on dynamic degrees
 * @author Jonas Alker
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "heur_dynamicdeg.h"
#include "graph.h"
#include "struct_probdata.h"
#include "scip/scip_sol.h"
#include "probdata_bacs.h"

#define HEUR_NAME            "dynamicdeg"
#define HEUR_DESC            "greedy search for feasible primal solutions based on dynamic node degrees"
#define HEUR_DISPCHAR        'D'
#define HEUR_PRIORITY        10000
#define HEUR_FREQ            1
#define HEUR_FREQOFS         0
#define HEUR_MAXDEPTH        -1
#define HEUR_TIMING          ( SCIP_HEURTIMING_BEFOREPRESOL | SCIP_HEURTIMING_BEFORENODE )
#define HEUR_USESSUBSCIP     FALSE /**< does the heuristic use a secondary SCIP instance? */


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecDynamicdeg)
{
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( result != nullptr );

   SCIPdebugMsg(scip, "Running dynamicdeg heuristic ...\n");

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   const Graph* G = probdata->G;

#ifndef NDEBUG
   SCIP_CALL( BACScheckLocalDegrees(scip, probdata) );
#endif

   int* degs;
   SCIP_Real* quot;
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &degs, probdata->localdegrees, probdata->n) ); /*lint !e530 */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &quot, probdata->n) );

   for (size_t i = 0; i < probdata->n; ++i)
   {
      if ( degs[i] != 0 )
         quot[i] = SCIPvarGetObj(probdata->vars[i]) / degs[i];
      else
         quot[i] = SCIP_REAL_MAX;
   }

   // here the changes in comparison to greedydeg begin:
   SCIP_SOL* sol;
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   // set containing all already fixed variables,
   std::vector<bool> blockset(probdata->n, false);
   unsigned int blockcount = 0;

   // stored fixed variables and adjust degrees similar to computeLocalDegrees()
   for (size_t i = 0; i < probdata->n; i++)
   {
      SCIP_VAR* var = probdata->vars[i];

      if ( SCIPvarGetUbLocal(var) < 0.5 && ! blockset[i] )
      {
         blockset[i] = true;
         blockcount++;
      }
      else if ( SCIPvarGetLbLocal(var) > 0.5 )
      {
         assert( ! blockset[i] );

         SCIP_CALL( SCIPsetSolVal(scip, sol, var, 1.0) );
         blockset[i] = true;
         blockcount++;

         AdjacencyIterator ait, aend;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(i, *G); ait != aend; ++ait)
         {
            Vertex w = *ait;

            assert( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 );
            if ( ! blockset[w] )
            {
               ++blockcount;
               blockset[w] = true;
            }
         }
      }
   }

   // try remaining variables in greedy order
   for (size_t iter = 0 ; iter < probdata->n ; ++iter)
   {
      assert( blockcount <= probdata->n );

      if ( blockcount == probdata->n )
         break;

      size_t bestnode = 0;
      SCIP_Real bestquot = - SCIPinfinity(scip);

      for (size_t i = 0 ; i < probdata->n ; ++i)
      {
         SCIP_VAR* var = probdata->vars[i];
         size_t nodeind = (size_t) SCIPvarGetData(var);

         if ( blockset[nodeind] )
            continue;

         if ( quot[i] > bestquot )
         {
            bestnode = nodeind;
            bestquot = quot[i];
         }
      }

      assert( SCIPisFeasGE(scip, bestquot, 0.0) );

      SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->vars[bestnode], 1.0) );
      blockset[bestnode] = true;
      ++blockcount;

      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(bestnode, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;
         if ( blockset[w] )
            continue;

         ++blockcount;
         blockset[w] = true;

         // adjust quot for neighbors of fixed to zero variables
         AdjacencyIterator bit, bend;
         for (boost::tie(bit, bend) = boost::adjacent_vertices(w, *G); bit != bend; ++bit)
         {
            Vertex x = *bit;
            if ( blockset[x] )
               continue;

            assert( degs[x] > 0 );
            --degs[x];

            if ( degs[x] != 0 )
               quot[x] = SCIPvarGetObj(probdata->vars[x]) / degs[x];
            else
               quot[x] = SCIP_REAL_MAX;
         }
      }
   }
   assert( blockcount == probdata->n );

   SCIPfreeBlockMemoryArray(scip, &quot, probdata->n);
   SCIPfreeBlockMemoryArray(scip, &degs, probdata->n);

   SCIP_Bool stored;
   SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, TRUE, FALSE, &stored) );

   if ( stored )
   {
      *result = SCIP_FOUNDSOL;
      SCIPdebugMsg(scip, "Dynamic Degree heuristic found solution.\n");
   }
   else
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMsg(scip, "Dynamic Degree heuristic did not find new solution.\n");
   }

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the dynamicdeg primal heuristic and includes it in SCIP */
SCIP_RETCODE BACSincludeHeurDynamicdeg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   constexpr SCIP_HEURDATA* heurdata = nullptr;
   SCIP_HEUR* heur = nullptr;

   /* create dynamicdeg primal heuristic data */
   SCIP_CALL(SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecDynamicdeg, heurdata));
   assert( heur != nullptr );

   return SCIP_OKAY;
}
