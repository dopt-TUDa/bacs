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

/**@file   branch_degree.cpp
 * @brief  branching rule based on local degree
 * @author Annika Jaeger
 * @author Erik Jansen
 * @author Jonas Alker
 * @author Marc Pfetsch
 */

#include <assert.h>

#include "branch_degree.h"
#include "struct_probdata.h"


#define BRANCHRULE_NAME            "degree"
#define BRANCHRULE_DESC            "local degree branching rule, allows branching on variables or neigborhoods"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define DEFAULT_MINDEGBRANCHING   FALSE      //< branch on node with smallest local degree instead of largest?
#define DEFAULT_NEIGHBRANCHING    FALSE      //< branch on neighborhood instead of variable, creating only 1-branches (but many more)


/*
 * Data structures
 */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Bool             mindegbranching;    //< branch on node with smallest local degree instead of largest?
   SCIP_Bool             neighbranching;     //< branch on neighborhood instead of variable, creating only 1-branches (but many more)
};


/*
 * Local methods
 */

/** perform maximal degree branching */
static
SCIP_RETCODE branchDegree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULEDATA*  branchruledata,     /**< data of branching rule */
   SCIP_RESULT*          result,             /**< result pointer */
   SCIP_Bool             islpbranching       /**< perform LP branching */
   )
{
   assert( scip != nullptr );
   assert( branchruledata != nullptr );
   assert( result != nullptr );

   *result = SCIP_DIDNOTFIND;

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );
   assert( probdata->n > 0 );

   const Graph G = *probdata->G;
   int* degs = probdata->localdegrees;

   // find vertex to branch on based on degree
   Vertex v = probdata->n;
   int degree = -1;
   if ( branchruledata->mindegbranching )
      degree = INT_MAX;

   for (size_t i = 0; i < probdata->n; ++i)
   {
      // skip fixed nodes
      if ( degs[i] == -1 )
         continue;

      assert( SCIPvarGetLbLocal(probdata->vars[i]) < 0.5 && SCIPvarGetUbLocal(probdata->vars[i]) > 0.5 );

      // when doing lp branching, only branch on non-integral values
      if ( islpbranching && SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, nullptr, probdata->vars[i])) )
         continue;

      if ( degs[i] == 0 )
      {
#ifndef NDEBUG
         // assert all neighbors are fixed to zero
         AdjacencyIterator ait, aend;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(i, G); ait != aend; ++ait)
         {
            assert( probdata->localdegrees[*ait] < 0 );
            assert( SCIPvarGetUbLocal(probdata->vars[*ait]) < 0.5 );
         }
#endif
         SCIP_Bool infeasible;
         SCIP_Bool fixed;
         SCIP_CALL( SCIPfixVar(scip, probdata->vars[i], 1.0, &infeasible, &fixed) );

         assert( ! infeasible );
         assert( fixed );

         *result = SCIP_REDUCEDDOM;
         continue;
      }

      // check for larger degree
      if ( ! branchruledata->mindegbranching && degs[i] > degree )
      {
         degree = degs[i];
         v = i;
      }

      // or check for smaller degree
      if ( branchruledata->mindegbranching && degs[i] < degree )
      {
         degree = degs[i];
         v = i;
      }
   }
   if ( *result == SCIP_REDUCEDDOM )
      return SCIP_OKAY;

   assert( v < probdata->n );
   assert( degree > 0 && degree < (int) probdata->n);
   assert( ! islpbranching || ! SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, nullptr, probdata->vars[v])) );

   // actually perform branching
   if ( ! branchruledata->neighbranching )
   {
      // simply branch on variable v
      SCIP_CALL( SCIPbranchVar(scip, probdata->vars[v], nullptr, nullptr, nullptr) );

      SCIPdebugMsg(scip, "degree branching rule branched on variable %zu\n", v);
      *result = SCIP_BRANCHED;
   }
   else
   {
      // in an optimal solution either v has value 1 or at least one neighbor of v
      size_t nbranch = 0;

      // create a 1-child for node v
      SCIP_NODE* childnode;
      SCIP_CALL( SCIPcreateChild(scip, &childnode, 0.0, SCIPgetLocalTransEstimate(scip)) );

      SCIP_CALL( SCIPchgVarLbNode(scip, childnode, probdata->vars[v], 1.0) );
      ++nbranch;

      // create a 1-child for every neighbor of node v
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices((Vertex) v, G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         assert( SCIPvarGetLbLocal(probdata->vars[w]) < 0.5 ); // neighbors of v cannot be already fixed to one

         // dont create a branch for already fixed neighbors
         if ( degs[w] == -1 )
            continue;

         assert( SCIPvarGetUbLocal(probdata->vars[w]) > 0.5 );

         SCIP_NODE* neighbornode;
         SCIP_CALL( SCIPcreateChild(scip, &neighbornode, 0.0, SCIPgetLocalTransEstimate(scip)) );

         SCIP_CALL( SCIPchgVarUbNode(scip, neighbornode, probdata->vars[v], 0.0) );
         SCIP_CALL( SCIPchgVarLbNode(scip, neighbornode, probdata->vars[w], 1.0) );

#ifndef NDEBUG
         // neighbors of w cannot be already fixed to one
         AdjacencyIterator cit, cend;
         for (boost::tie(cit, cend) = boost::adjacent_vertices((Vertex) w, G); cit != cend; ++cit)
         {
            Vertex z = *cit;
            assert( SCIPvarGetLbLocal(probdata->vars[z]) < 0.5 );
         }
#endif
         // set all later vertices in neighborhood to zero to create disjoint childnodes
         AdjacencyIterator bit, bend;
         bend = aend;
         for (bit = ait + 1; bit != bend; ++bit)
         {
            Vertex x = *bit;

            if ( degs[x] == -1 )
               continue;

            assert( SCIPvarGetLbLocal(probdata->vars[x]) < 0.5 && SCIPvarGetUbLocal(probdata->vars[x]) > 0.5 );

            SCIP_CALL( SCIPchgVarUbNode(scip, neighbornode, probdata->vars[x], 0.0) );
         }
         ++nbranch;
      }

      assert( nbranch == (size_t) degree + 1 );
      SCIPdebugMsg(scip, "degree branching rule created %zu branches\n", nbranch);
      *result = SCIP_BRANCHED;
   }
   return SCIP_OKAY;
}


/*
 * Callback methods of branching rule
 */

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpDegree)
{
   assert( scip != nullptr );
   assert( branchrule != nullptr );
   assert( result != nullptr );

   SCIP_BRANCHRULEDATA* branchruledata = SCIPbranchruleGetData(branchrule);
   assert( branchruledata != nullptr );

   SCIPdebugMsg(scip, "Running degree branching rule ...\n");

   SCIP_CALL( branchDegree(scip, branchruledata, result, TRUE) );

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsDegree)
{
   assert( scip != nullptr );
   assert( branchrule != nullptr );
   assert( result != nullptr );

   SCIP_BRANCHRULEDATA* branchruledata = SCIPbranchruleGetData(branchrule);
   assert( branchruledata != nullptr );

   SCIPdebugMsg(scip, "Running pseudo degree branching rule ...\n");

   SCIP_CALL( branchDegree(scip, branchruledata, result, FALSE) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeDegree)
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

/** creates the deg branching rule and includes it in SCIP */
SCIP_RETCODE BACSincludeBranchruleDegree(
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
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpDegree) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsDegree) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeDegree) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/" BRANCHRULE_NAME "/mindegbranching",
         "branch on node with smallest local degree instead of largest?",
         &branchruledata->mindegbranching, FALSE, DEFAULT_MINDEGBRANCHING, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/" BRANCHRULE_NAME "/neighbranching",
         "branch on neighborhood instead of variable, creating only 1-branches (but many more)",
         &branchruledata->neighbranching, FALSE, DEFAULT_NEIGHBRANCHING, nullptr, nullptr) );

   return SCIP_OKAY;
}
