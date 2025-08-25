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

/**@file   prop_dominance.cpp
 * @ingroup PROPAGATORS
 * @brief  propagator to fix nodes to zero whose neigborhood contains neighborhood of a node that has been branched to 0
 * @author Annika Jaeger
 *
 * This propagator is only correct if it runs directly after branching before any other propagation methods were executed.
 * To ensure this, we use PROP_TIMING = SCIP_PROPTIMING_BEFORELP and a very high PROP_PRIORITY.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "prop_dominance.h"
#include "graph.h"
#include "struct_probdata.h"

#include <scip/prop_symmetry.h>


#define PROP_NAME                  "dominance"
#define PROP_DESC                  "propagator to fix nodes to zero whose neigborhood contains neighborhood of a node that has been branched to 0"
#define PROP_TIMING                 SCIP_PROPTIMING_BEFORELP
#define PROP_PRIORITY             100000000  //!< priority of the propagator (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PROP_FREQ                         2  //!< propagation frequency
#define PROP_DELAY                    FALSE  //!< should propagation method be delayed, if other propagators found reductions?

#define DEFAULT_ONLYFORNOSYMMETRY      TRUE  //!< whether the propagator should only be run in the tree if there are no symmetries

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_Bool             probingonly;        //!< whether the propagator should only run in probing
   SCIP_Bool             onlyfornosymmetry;  //!< whether the propagator should only be run in the tree if there are no symmetries
   SCIP_Bool             enabled;            //!< whether the propagator is enabled
   SCIP_Node*            lastnode;           //!< last node the propagator has been used at
};


/*
 * Local methods
 */

//! perform propagation
static
SCIP_RETCODE propagateDominance(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_VAR*             svar,               /**< variable that was branched to zero */
   SCIP_Bool&            cutoff,             /**< whether we detected a cutoff */
   size_t&               nfixed              /**< number of variables fixed */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( svar != nullptr );

   cutoff = FALSE;
   nfixed = 0;

   const Graph* G = probdata->G;

   Vertex s = (Vertex) SCIPvarGetData(svar);
   assert( s != Graph::null_vertex() );

   // compute size of neighborhood of s
   size_t nneigh = 0;
   AdjacencyIterator ait, aend;
   for (boost::tie(ait, aend) = boost::adjacent_vertices(s, *G); ait != aend; ++ait)
   {
      if ( SCIPvarGetUbLocal(SCIPvarGetTransVar(probdata->vars[*ait])) < 0.5 )
         continue;

      // If v has a neighbor fixed to 1.0, then s should have been fixed to 0.0,
      // however, this does not seem to be the case here (we branched on s), so we exit.
      if ( SCIPvarGetLbLocal(SCIPvarGetTransVar(probdata->vars[*ait])) > 0.5 )
         return SCIP_OKAY;

      ++nneigh;
   }

   // initialize vector counting how often nodes appear in neighborhoods of neighbors of s
   std::vector<size_t> domnodes(probdata->n, 0);

   // loop over neighbors of s
   for (boost::tie(ait, aend) = boost::adjacent_vertices(s, *G); ait != aend; ++ait)
   {
      Vertex neigh = *ait;
      SCIP_VAR* neighvar = SCIPvarGetTransVar(probdata->vars[neigh]);

      // skip fixed node
      if ( SCIPvarGetUbLocal(neighvar) < 0.5 || SCIPvarGetLbLocal(neighvar) > 0.5 )
         continue;

      // loop over neighbors of neigh
      AdjacencyIterator ait2, aend2;
      for (boost::tie(ait2, aend2) = boost::adjacent_vertices(neigh, *G); ait2 != aend2; ++ait2)
      {
         Vertex t = *ait2;
         SCIP_VAR* tvar = SCIPvarGetTransVar(probdata->vars[t]);

         // skip fixed node
         if ( SCIPvarGetUbLocal(tvar) < 0.5 || SCIPvarGetLbLocal(tvar) > 0.5 )
            continue;
         assert( s != t );

         ++domnodes[t];
         assert( domnodes[t] <= nneigh );

         if ( domnodes[t] == nneigh )
         {
            SCIPdebugMsg(scip, " -> can fix <%s> to 0, because its neighborhood dominates the neighborhood of <%s> which is branched to 0.\n",
               SCIPvarGetName(tvar), SCIPvarGetName(svar));

            // neighborhood of s is subset of neighborhood of t and s is branched to 0. This can only be optimal if t is also 0 -> fix t to 0.0
            SCIP_Bool tightened = FALSE;
            SCIP_CALL( SCIPinferBinvarProp(scip, tvar, FALSE, prop, (int) s, &cutoff, &tightened) );
            ++nfixed;
            assert( ! cutoff );
            assert( tightened );
         }
      }
   }

   return SCIP_OKAY;
}



/*
 * Callback methods of propagator
 */

//! solving process initialization method of propagator (called when branch and bound process is about to begin)
static
SCIP_DECL_PROPINITSOL(propInitsolDominance)
{
   assert( prop != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   // if propagation should be performed in the tree, we currently turn off symmetry handling to avoid conflicts, unless we turn it off below anyways
   if ( SCIPpropGetFreq(prop) >= 0 && ! propdata->onlyfornosymmetry )
   {
      if ( ! propdata->probingonly )
      {
         SCIP_RETCODE retcode = SCIPsetIntParam(scip, "propagating/symmetry/freq", -1);
         if ( retcode == SCIP_OKAY )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "Turning off symmetry handling to avoid conflicts with dominance propagation.\n");
         else
            assert( retcode == SCIP_PARAMETERUNKNOWN );
      }
   }

   return SCIP_OKAY;
}


//! destructor of propagator to free user data (called when SCIP is exiting)
static
SCIP_DECL_PROPFREE(propFreeDominance)
{
   assert( prop != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   SCIPfreeBlockMemory(scip, &propdata);

   return SCIP_OKAY;
}


//! execution method of propagator
static
SCIP_DECL_PROPEXEC(propExecDominance)
{  /*lint --e{715}*/
   assert( scip != nullptr );
   assert( result != nullptr );
   assert( prop != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   *result = SCIP_DIDNOTRUN;
   // if propagation should only run in probing
   if ( SCIPinProbing(scip) )
      return SCIP_OKAY;

   // do not run if disabled
   if ( ! propdata->enabled )
      return SCIP_OKAY;

   // possibly only run in the tree if there are no symmetries
   if ( propdata->onlyfornosymmetry && SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING && SCIPgetSymmetryNGenerators(scip) > 0 )
      return SCIP_OKAY;

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   // get current node
   SCIP_NODE* node = SCIPgetCurrentNode(scip);

   // do not run if propagator was already applied at this node
   if ( node == propdata->lastnode )
      return SCIP_OKAY;

   // get number of domain changes at current node
   int nconsprop;
   int nprop;
   SCIPnodeGetNDomchg(node, nullptr, &nconsprop, &nprop);

   // do not run if constraint handlers or propagators have already found reductions (may happen when node gets repropagated)
   if ( nconsprop > 0 || nprop > 0 )
      return SCIP_OKAY;

   // get branching variable and its value
   int branchvarssize = 1;
   SCIP_VAR* branchvars;
   SCIP_Real branchbounds;
   SCIP_BOUNDTYPE boundtypes;
   int nbranchvars;

   SCIPnodeGetParentBranchings(node, &branchvars, &branchbounds, &boundtypes, &nbranchvars, branchvarssize);

   // do nothing if no node has been branched on
   if ( nbranchvars == 0 )
      return SCIP_OKAY;

   // do nothing if it was branched on one more than one variable
   if ( nbranchvars > branchvarssize )
      return SCIP_OKAY;

   // skip if node is not branched to 0
   if ( boundtypes != SCIP_BOUNDTYPE_UPPER || ! SCIPisZero(scip, branchbounds) )
      return SCIP_OKAY;

   propdata->lastnode = node;
   SCIP_Bool cutoff;
   size_t nfixed;
   SCIPdebugMsg(scip, "running dominance propagation.\n");
   SCIP_CALL( propagateDominance(scip, prop, probdata, branchvars, cutoff, nfixed) );

   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( nfixed > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropDominance)
{
   assert( scip != nullptr );
   assert( prop != nullptr );
   assert( infervar != nullptr );
   assert( bdchgidx != nullptr );
   assert( result != nullptr );

   SCIPdebugMsg(scip, "Propagation resolution method of <%s>.\n", PROP_NAME);
   *result = SCIP_DIDNOTFIND;

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );
   assert( 0 <= inferinfo && inferinfo < (int) probdata->n );
   assert( boundtype == SCIP_BOUNDTYPE_UPPER );
   const Graph* G = probdata->G;

   // neighborhood of s is contained in the one of t
   Vertex s = (Vertex) (size_t) inferinfo;
   Vertex t = (Vertex) SCIPvarGetData(infervar);
   SCIP_VAR* svar = SCIPvarGetTransVar(probdata->vars[s]);
   SCIP_VAR* tvar = infervar;
   SCIPdebugMsg(scip, "neighborhood of <%s> was contained in <%s> when <%s> was branched to zero.\n", SCIPvarGetName(svar), SCIPvarGetName(tvar), SCIPvarGetName(svar));

   // check preconditions for original propagation
   assert( SCIPvarGetTransVar(probdata->vars[t]) == tvar );
   assert( ! boost::edge(s, t, *G).second );
   // s was branched to zero
   assert( SCIPgetVarUbAtIndex(scip, svar, bdchgidx, FALSE) < 0.5 );
   assert( SCIPgetVarLbAtIndex(scip, svar, bdchgidx, FALSE) < 0.5 );
   // The following assert might fail, because conflict analysis can determine that the variables could have been fixed
   // to 0 higher up in the tree.
   // assert( SCIPgetVarUbAtIndex(scip, tvar, bdchgidx, FALSE) > 0.5 );
   assert( SCIPgetVarLbAtIndex(scip, tvar, bdchgidx, FALSE) < 0.5 );

   SCIP_CALL( SCIPaddConflictUb(scip, svar, bdchgidx) );
   *result = SCIP_SUCCESS;

   // set up iterator for neighborhood of t
   AdjacencyIterator tit, tend;
   boost::tie(tit, tend) = boost::adjacent_vertices(t, *G);

   // loop over neighborhood of s
   AdjacencyIterator sit, send;
   for (boost::tie(sit, send) = boost::adjacent_vertices(s, *G); sit != send; ++sit)
   {
      assert ( *sit != t );

      if ( tit != tend )
      {
         // ignore common neighbors
         if ( *sit == *tit )
            continue;

         // s has neighbor which t does not have - thus this variable needs to be fixed to 0 for a valid reduction
         if ( *sit < *tit )
         {
            SCIP_VAR* var = probdata->vars[*sit];
            assert( SCIPvarGetUbAtIndex(var, bdchgidx, FALSE) < 0.5 );
            SCIPdebugMsg(scip, "fixings in neighborhood of <%s>: <%s> = 0.\n", SCIPvarGetName(svar), SCIPvarGetName(var));
            SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
            continue;
         }
      }

      // ignore neighbors of t which s does not have
      while ( tit != tend && *sit > *tit )
         ++tit;
   }

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

//! creates the dominance propagator and includes it in SCIP
SCIP_RETCODE BACSincludePropDominance(
   SCIP*                 scip                //! SCIP data structure
   )
{
   SCIP_PROPDATA* propdata = nullptr;
   SCIP_PROP* prop = nullptr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );
   propdata->enabled = TRUE;
   propdata->lastnode = nullptr;

   // include propagator
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecDominance, propdata) );
   assert( prop != nullptr );

   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeDominance) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolDominance) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropDominance) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/onlyfornosymmetry",
         "whether the propagator should only be run in the tree if there are no symmetries",
         &propdata->onlyfornosymmetry, FALSE, DEFAULT_ONLYFORNOSYMMETRY, nullptr, nullptr) );

   return SCIP_OKAY;
}
