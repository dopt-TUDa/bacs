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

/**@file   heur_greedyrounding.cpp
 * @ingroup DEFPLUGINS_HEUR
 * @brief  greedy rounding primal heuristic
 * @author Erik Jansen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "heur_greedyrounding.h"
#include "graph.h"
#include "struct_probdata.h"
#include "scip/scip_sol.h"
#include "probdata_bacs.h"

#define HEUR_NAME            "greedyrounding"
#define HEUR_DESC            "greedy search for feasible primal solutions based on rounded LP solution"
#define HEUR_DISPCHAR        'G'
#define HEUR_PRIORITY        100000
#define HEUR_FREQ            1
#define HEUR_FREQOFS         0
#define HEUR_MAXDEPTH        -1
#define HEUR_TIMING          SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP     FALSE /**< does the heuristic use a secondary SCIP instance? */


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecGreedyRounding)
{
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( result != nullptr );
   assert( SCIPhasCurrentNodeLP(scip) );

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand or if relaxation solution is available */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is bigger than the cutoff bound */
   if ( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Running greedy rounding search heuristic ...\n");

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   const Graph* G = probdata->G;
   SCIP_VAR** vars;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, probdata->vars, probdata->n) ); /*lint !e530 */

   SCIP_Real* primsol;
   SCIP_CALL( SCIPallocBufferArray(scip, &primsol, probdata->n) );

   for (size_t i = 0; i < probdata->n; ++i)
   {
      SCIP_VAR* var = probdata->vars[i];
      primsol[i] = SCIPgetSolVal(scip, nullptr, var);
      if ( SCIPvarGetUbLocal(var) > 0.5 && SCIPvarGetLbLocal(var) < 0.5 )
      {
         if ( primsol[i] > 0.5 )
            primsol[i] = 1.0;
         else
            primsol[i] = 0.0;
      }
   }

#ifndef NDEBUG
   SCIP_CALL( BACScheckLocalDegrees(scip, probdata) );
#endif

   int* degs = nullptr;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &degs, probdata->localdegrees, probdata->n) );
   assert( degs != nullptr );

   SCIPsortDownRealIntPtr(primsol, degs, (void**) vars, (int) probdata->n);

#ifndef NDEBUG
   SCIP_Real roundval = 0.0;
   SCIP_Real oldroundval;

   // assert the order of sorted vars in descending order
   for (size_t i = 0; i < probdata->n; i++)
   {
      oldroundval = roundval;
      roundval = primsol[i];

      assert( ( SCIPisFeasEQ(scip, 0, roundval) || SCIPisFeasEQ(scip, 1, roundval) ) );

      if ( i != 0 )
         assert( SCIPisFeasLE(scip, roundval, oldroundval) );
   }
#endif

   size_t vit = 0; // array iterator
   int l = 0;   // length of array of same size

   while ( vit < probdata->n )
   {
      if ( (vit != probdata->n - 1) && SCIPisFeasEQ(scip, primsol[vit], primsol[vit + 1]) )
         l++;
      else
      {
         int baseidx = (int) vit - l;
         SCIPsortIntPtr(&degs[baseidx], (void**) &vars[baseidx], l+1);

#ifndef NDEBUG
         // assert the order of sorted vars in descending order
         for (int i = 1; i < l; i++)
         {
            assert( degs[baseidx + i - 1] <= degs[baseidx + i] );
         }
#endif

         l = 0;
      }
      ++vit;
   }

   SCIPfreeBufferArray(scip, &degs);
   SCIPfreeBufferArray(scip, &primsol);

   SCIP_SOL* sol;
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   // set containing all already fixed variables
   std::vector<bool> blockset(probdata->n,false);
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
         ++blockcount;
      }
      else if ( SCIPvarGetLbLocal(var) > 0.5 )
      {
         assert( ! blockset[nodeind] );

         SCIP_CALL( SCIPsetSolVal(scip, sol, var, 1.0) );
         blockset[nodeind] = true;
         ++blockcount;

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

   // try left variables in greedy order
   for (size_t i = 0; i < probdata->n; i++)
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
      ++blockcount;

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

#ifdef SCIP_DEBUG
   SCIP_Bool feasible;
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, TRUE, TRUE, FALSE, &feasible) );
   assert( feasible);
#endif

   SCIPfreeBufferArray(scip, &vars);

   SCIP_Bool stored;
   SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, TRUE, FALSE, &stored) );

   if ( stored )
   {
      *result = SCIP_FOUNDSOL;
      SCIPdebugMsg(scip, "Greedy Rounding heuristic found solution.\n");
   }
   else
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMsg(scip, "Greedy Rounding heuristic did not find new solution.\n");
   }

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the greedy rounding primal heuristic and includes it in SCIP */
SCIP_RETCODE BACSincludeHeurGreedyRounding(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   constexpr SCIP_HEURDATA* heurdata = nullptr;
   SCIP_HEUR* heur = nullptr;

   /* create greedyrounding primal heuristic data */
   SCIP_CALL(SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecGreedyRounding, heurdata));
   assert( heur != nullptr );

   return SCIP_OKAY;
}
