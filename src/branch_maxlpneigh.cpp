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

/**@file   branch_maxlpneigh.cpp
 * @brief  maximal lp value in neighborhood branching rule
 * @author Annika Jaeger
 * @author Erik Jansen
 * @author Jonas Alker
 * @author Marc Pfetsch
 */

#include <assert.h>

#include "branch_maxlpneigh.h"
#include "graph.h"
#include "struct_probdata.h"


#define BRANCHRULE_NAME            "maxlpneigh"
#define BRANCHRULE_DESC            "maximal lp -value in neighborhood branching rule in induced subgraph of unfixed nodes"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0


/*
 * Callback methods of branching rule
 */

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpMaxlpneigh)
{
   assert( scip != nullptr );
   assert( branchrule != nullptr );
   assert( result != nullptr );

   SCIPdebugMsg(scip, "Running maxlp neighborhood branching rule ...\n");

   *result = SCIP_DIDNOTFIND;

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );
   assert( probdata->n > 0 );

   const Graph* G = probdata->G;

   SCIP_Real* neighlpval;
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &neighlpval, probdata->n) );

   EdgeIterator eit, eend;
   for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
   {
      Vertex s = boost::source(*eit, *G);
      Vertex t = boost::target(*eit, *G);

      assert( s < probdata->n );
      assert( t < probdata->n );

      if ( SCIPvarGetUbLocal(probdata->vars[s]) > 0.5 && SCIPvarGetLbLocal(probdata->vars[s]) < 0.5
         && SCIPvarGetUbLocal(probdata->vars[t]) > 0.5 && SCIPvarGetLbLocal(probdata->vars[t]) < 0.5 )
      {
         neighlpval[t] += SCIPgetSolVal(scip, nullptr, probdata->vars[s]);
         neighlpval[s] += SCIPgetSolVal(scip, nullptr, probdata->vars[t]);
      }
   }

   // find vertex with maximal lp value sum in neighborhood
   Vertex v = probdata->n;
   SCIP_Real maxlpneigh = 0.0;
   for (size_t i = 0; i < probdata->n; ++i)
   {
      if ( SCIPisGT(scip, neighlpval[i], maxlpneigh) && ! SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, nullptr, probdata->vars[i])) )
      {
         maxlpneigh = neighlpval[i];
         v = i;
      }
   }

   if ( SCIPisGT(scip, maxlpneigh, 0.0) )
   {
      assert( v < probdata->n );
      assert( ! SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, nullptr, probdata->vars[v])) );

      // actually perform branching
      SCIP_CALL( SCIPbranchVar(scip, probdata->vars[v], nullptr, nullptr, nullptr) );
      *result = SCIP_BRANCHED;
   }

   SCIPfreeBlockMemoryArray(scip, &neighlpval, probdata->n);

   return SCIP_OKAY;
}


/*
 * branching rule specific interface methods
 */

/** creates the maxlpneigh branching rule and includes it in SCIP */
SCIP_RETCODE BACSincludeBranchruleMaxlpneigh(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != nullptr );

   constexpr SCIP_BRANCHRULEDATA* branchruledata = nullptr;
   SCIP_BRANCHRULE* branchrule = nullptr;

   /* create maxlpneigh branching rule data */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert( branchrule != nullptr );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpMaxlpneigh) );

   return SCIP_OKAY;
}
