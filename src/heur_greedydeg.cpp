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

/**@file   heur_greedydeg.cpp
 * @brief  greedy degree primal heuristic with smart degree computation
 * @author Erik Jansen
 * @author Jonas Alker
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "heur_greedydeg.h"
#include "graph.h"
#include "struct_probdata.h"
#include "scip/scip_sol.h"
#include "probdata_bacs.h"

#define HEUR_NAME            "greedydeg"
#define HEUR_DESC            "greedy search for feasible primal solutions based on node degree (smart degree computation)"
#define HEUR_DISPCHAR        'G'
#define HEUR_PRIORITY        10000
#define HEUR_FREQ            1
#define HEUR_FREQOFS         0
#define HEUR_MAXDEPTH        -1
#define HEUR_TIMING          ( SCIP_HEURTIMING_BEFOREPRESOL | SCIP_HEURTIMING_BEFORENODE )
#define HEUR_USESSUBSCIP     FALSE /**< does the heuristic use a secondary SCIP instance? */


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecGreedyDeg)
{
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( result != nullptr );

   SCIPdebugMsg(scip, "Running greedy degree search heuristic ...\n");

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   const Graph* G = probdata->G;

#ifndef NDEBUG
   SCIP_CALL( BACScheckLocalDegrees(scip, probdata) );
#endif

   int* degs;
   SCIP_Real* quot;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &degs, probdata->localdegrees, probdata->n) ); /*lint !e530 */
   SCIP_CALL( SCIPallocBufferArray(scip, &quot, probdata->n) );

   for (size_t i = 0; i < probdata->n; ++i)
   {
      if ( degs[i] != 0 )
         quot[i] = SCIPvarGetObj(probdata->vars[i]) / degs[i];
      else
         quot[i] = SCIP_REAL_MAX;
   }

   SCIP_VAR** vars;
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &vars, probdata->vars, probdata->n) ); /*lint !e530 */

   SCIPsortDownRealPtr(quot, (void**) vars, (int) probdata->n);

   SCIPfreeBufferArray(scip, &quot);
   SCIPfreeBufferArray(scip, &degs);

   SCIP_SOL* sol;
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   // set containing all already fixed variables
   std::vector<bool> blockset(probdata->n, false);
   unsigned int blockcount = 0;

   // stored fixed variables
   for (size_t i = 0; i < probdata->n; i++)
   {
      assert( blockcount <= probdata->n );

      SCIP_VAR* var = vars[i];
      size_t nodeind = (size_t) SCIPvarGetData(var);

      if ( SCIPvarGetUbLocal(var) < 0.5 && ! blockset[nodeind] )
      {
         blockset[nodeind] = true;
         blockcount++;
      }
      else if ( SCIPvarGetLbLocal(var) > 0.5 )
      {
         assert( ! blockset[nodeind] );

         SCIP_CALL( SCIPsetSolVal(scip, sol, var, 1.0) );
         blockset[nodeind] = true;
         blockcount++;

         AdjacencyIterator ait, aend;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(nodeind, *G); ait != aend; ++ait)
         {
            Vertex w = *ait;
            if ( ! blockset[w] )
            {
               ++blockcount;
               blockset[w] = true;
            }
         }
      }
   }

   // try remaining variables in greedy order
   for (size_t i = 0 ; i < probdata->n ; i++)
   {
      assert( blockcount <= probdata->n );

      if ( blockcount == probdata->n )
         break;

      SCIP_VAR* var = vars[i];
      size_t nodeind = (size_t) SCIPvarGetData(var);

      if ( blockset[nodeind] )
         continue;

      SCIP_CALL( SCIPsetSolVal(scip, sol, var, 1.0) );
      blockset[nodeind] = true;
      blockcount++;

      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(nodeind, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;
         if ( ! blockset[w] )
         {
            ++blockcount;
            blockset[w] = true;
         }
      }
   }
   assert( blockcount == probdata->n );

   SCIPfreeBlockMemoryArray(scip, &vars, probdata->n);

   SCIP_Bool stored;
   SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, TRUE, FALSE, &stored) );

   if ( stored )
   {
      *result = SCIP_FOUNDSOL;
      SCIPdebugMsg(scip, "Greedy Degree heuristic found solution.\n");
   }
   else
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMsg(scip, "Greedy Degree heuristic did not find new solution.\n");
   }

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the greedy degree primal heuristic and includes it in SCIP */
SCIP_RETCODE BACSincludeHeurGreedyDeg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   constexpr SCIP_HEURDATA* heurdata = nullptr;
   SCIP_HEUR* heur = nullptr;

   /* create greedydeg primal heuristic data */
   SCIP_CALL(SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecGreedyDeg, heurdata));
   assert( heur != nullptr );

   return SCIP_OKAY;
}
