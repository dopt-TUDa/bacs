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

/**@file   prop_comp.cpp
 * @brief  propagator for handling connected components
 * @author Marc Pfetsch
 *
 * This propagator treats connected components. Modified from the corresponding SCIP constraint handler.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "prop_comp.h"
#include "probdata_bacs.h"

#include "blockmemshell/memory.h"
#include "scip/debug.h"
#include "scip/pub_cons.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_sol.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/prop_symmetry.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_dialog.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_pricer.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include <string.h>

#include "transfer_stats.h"
#include "sepa_neigh.h"
#include "prop_neighborhoods.h"
#include "heur_greedylp.h"
#include "heur_greedydeg.h"
#include "heur_greedyrounding.h"
#include "heur_tabu.h"
#include "heur_dynamicdeg.h"

#define PROP_NAME                  "comp"
#define PROP_DESC                  "connected components propagator"
#define PROP_TIMING                 SCIP_PROPTIMING_BEFORELP
#define PROP_PRIORITY               -9999999 //!< propagation priority
#define PROP_FREQ                         1 //!< propagation frequency
#define PROP_DELAY                     FALSE //!< should propagation method be delayed, if other propagators found reductions?
#define PROP_PRESOL_PRIORITY         -999999 //!< priority of the propagator (>= 0: before, < 0: after constraint handlers)
#define PROP_PRESOL_MAXROUNDS              1 //!< maximal number of propving rounds the propver participates in (-1: no limit)
#define PROP_PRESOLTIMING           SCIP_PRESOLTIMING_EXHAUSTIVE //!< timing of the presolving method (fast, medium, or exhaustive)

// default values for parameters
#define DEFAULT_MAXDEPTH              5      //!< maximum depth of a node to run components detection (-1: disable component detection during solving)
#define DEFAULT_NODELIMIT       10000LL      //!< maximum number of nodes to be solved in subproblems during solving
#define DEFAULT_MAXSIZE             100      //!< maximum size of components considered in the tree (-1: no limit)
#define DEFAULT_ONLYSMALL         FALSE      //!< only treat small components
#define DEFAULT_MINSIZE_SYM          70      //!< minimal size of component to include symmetry handling for

#define DEFAULT_ONLYFORNOSYMMETRY      TRUE  //!< whether the propagator should only be run in the tree if there are no symmetries

// other defines
#define SUBSCIP_OUTPUT_FREQ        1000      //!< frequency in which output is produced for subscips


/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_Longint          nodelimit;          //!< maximum number of nodes to be solved in subproblems during solving
   int                   maxdepth;           //!< maximum depth of a node to run components detection (-1: disable component detection during solving)
   int                   maxsize;            //!< maximal size of components to consider (-1: no limit)
   SCIP_Bool             enabled;            //!< whether the propagator is enabled
   SCIP_Bool             onlysmall;          //!< only treat small components
   int                   minsizesym;         //!< minimal size of component to include symmetry handling for
   size_t*               compsizes;          //!< sizes of found components in the tree for statistic output
   int                   ncompsizes;         //!< length of compsizes
   SCIP_Bool             onlyfornosymmetry;  //!< whether the propagator should only be run in the tree if there are no symmetries
};


/*
 * Callback methods of event handler for displaying information of subscip
 */

#define EVENTHDLR_NAME         "subscip-display"
#define EVENTHDLR_DESC         "event handler for displaying information of subscip"

//! event handler data
struct SCIP_EventhdlrData
{
   SCIP*                 origscip;           //!< pointer to original SCIP
   int                   compidx;            //!< index of component
};


//! destructor of event handler to free user data (called when SCIP is exiting)
static
SCIP_DECL_EVENTFREE(eventFreeSubscip)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert( scip != nullptr );
   assert( eventhdlr != nullptr );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0 );

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert( eventhdlrdata != nullptr );

   SCIPfreeMemory(scip, &eventhdlrdata);
   SCIPeventhdlrSetData(eventhdlr, nullptr);

   return SCIP_OKAY;
}


//! initialization method of event handler (called after problem was transformed)
static
SCIP_DECL_EVENTINIT(eventInitSubscip)
{  /*lint --e{715}*/
   assert( scip != nullptr );
   assert( eventhdlr != nullptr );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0 );

   // cat node change events
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, nullptr, nullptr) );

   return SCIP_OKAY;
}


//! execution method of event handler
static
SCIP_DECL_EVENTEXEC(eventExecSubscip)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   SCIPdebugMsg(scip, "Exec method of eventhdlr <%s>.\n", EVENTHDLR_NAME);
   assert( eventhdlr != nullptr );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0 );
   assert( event != nullptr );
   assert( scip != nullptr );

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert( eventhdlrdata != nullptr );

   SCIP_Longint nnodes = SCIPgetNNodes(scip);
   if ( nnodes % SUBSCIP_OUTPUT_FREQ == 0 )
   {
      SCIP* origscip = eventhdlrdata->origscip;
      SCIPverbMessage(origscip, SCIP_VERBLEVEL_HIGH, nullptr,
         "   (%.1fs) Solving component %d: #nodes = %" SCIP_LONGINT_FORMAT ", left = %d, dual = %g, primal = %g, gap = %.2f %%\n",
         SCIPgetSolvingTime(origscip), eventhdlrdata->compidx, nnodes, SCIPgetNNodesLeft(scip),
         SCIPgetDualbound(scip), SCIPgetPrimalbound(scip), 100.0 * SCIPgetGap(scip));
   }

   return SCIP_OKAY;
}


//! creates event handler for subscip output
static
SCIP_RETCODE SCIPincludeEventHdlrSubscip(
   SCIP*                 subscip,            //!< subSCIP
   SCIP*                 origscip,           //!< original SCIP
   int                   compidx             //!< component index
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata = nullptr;
   SCIP_EVENTHDLR* eventhdlr = nullptr;

   SCIP_CALL( SCIPallocMemory(subscip, &eventhdlrdata) );
   eventhdlrdata->origscip = origscip;
   eventhdlrdata->compidx = compidx;

   // include event handler into SCIP
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecSubscip, eventhdlrdata) );
   assert( eventhdlr != nullptr );

   //! set non fundamental callbacks via setter functions
   SCIP_CALL( SCIPsetEventhdlrFree(subscip, eventhdlr, eventFreeSubscip) );
   SCIP_CALL( SCIPsetEventhdlrInit(subscip, eventhdlr, eventInitSubscip) );

   return SCIP_OKAY;
}


/*
 * Local methods
 */

//! create a sub-SCIP
static
SCIP_RETCODE createSubscip(
   SCIP*                 scip,               //!< main SCIP data structure
   SCIP**                subscip,            //!< pointer to store created sub-SCIP
   int                   minsizesym,         //!< minimal size of component to include symmetry handling for
   int                   size                //!< size of component
   )
{
   SCIP_Bool success;
   SCIP_Bool copytables = FALSE;

   // create a new SCIP instance
   SCIP_CALL( SCIPcreate(subscip) );

#ifdef SCIP_MORE_DEBUG // we print statistics later, so we need to copy statistics tables
   copytables = TRUE;
#endif

   // copy plugins, we omit pricers (because we do not run if there are active pricers) and dialogs
   /* Boolean parameters: readers, pricers, conshdlrs, conflicthdlrs, presolvers, relaxators, separators,
    * cutselectors, propagators, heuristics, eventhdlrs, nodeselectors, branchrules, (iisfinders), displays, dialogs, tables,
    * exprhdlrs, nlpis, passmessagehdlr, valid
    */
#if SCIP_VERSION >= 1000
   SCIP_CALL( SCIPcopyPlugins(scip, *subscip, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, copytables,
         TRUE, TRUE, TRUE, TRUE, &success) );
#else
   SCIP_CALL( SCIPcopyPlugins(scip, *subscip, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, copytables,
         TRUE, TRUE, TRUE, TRUE, &success) );
#endif

   if ( ! success )
   {
      SCIP_CALL( SCIPfree(subscip) );
      *subscip = nullptr;
      return SCIP_OKAY;
   }

   if ( size >= minsizesym )
   {
      if ( SCIPfindProp(*subscip, "symmetry") == nullptr )
      {
         // because of a bug in SCIP, we add an empty root and display dialog (this can be remove in later SCIP versions)
         SCIP_DIALOG* root;
         SCIP_CALL( SCIPincludeDialog(*subscip, &root, nullptr, nullptr, nullptr, nullptr,
               "main", "empty main menu", TRUE, nullptr) );
         SCIP_CALL( SCIPsetRootDialog(*subscip, root) );
         assert( !SCIPdialogHasEntry(root, "display") );
         SCIP_DIALOG* submenu;
         SCIP_CALL( SCIPincludeDialog(*subscip, &submenu, nullptr, nullptr, nullptr, nullptr, "display", "display", TRUE, nullptr) );
         SCIP_CALL( SCIPaddDialogEntry(*subscip, root, submenu) );
         SCIP_CALL( SCIPreleaseDialog(*subscip, &submenu) );
         SCIP_CALL( SCIPreleaseDialog(*subscip, &root) );

         // now add symmetry propagator
         SCIP_CALL( SCIPincludePropSymmetry(*subscip) );
      }
   }

   SCIP_CALL( BACSincludeSepaNeigh(*subscip) );
   SCIP_CALL( BACSincludePropNeighborhoods(*subscip) );
   SCIP_CALL( BACSincludeHeurGreedyLP(*subscip) );
   SCIP_CALL( BACSincludeHeurGreedyDeg(*subscip) );
   SCIP_CALL( BACSincludeHeurGreedyRounding(*subscip) );
   SCIP_CALL( BACSincludeHeurTabu(*subscip) );
   SCIP_CALL( BACSincludeHeurDynamicdeg(*subscip) );

   // copy parameter settings
   SCIP_CALL( SCIPcopyParamSettings(scip, *subscip) );

   // disable presolving
   SCIP_CALL( SCIPsetPresolving(*subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   // disable output, unless in extended debug mode
#ifndef SCIP_MORE_DEBUG
   SCIP_CALL( SCIPsetIntParam(*subscip, "display/verblevel", 0) );
#endif

   return SCIP_OKAY;
}


//! create a stable set sub-problem
static
SCIP_RETCODE createStableSetSubproblem(
   SCIP*                 scip,               //!< original SCIP
   SCIP*                 subscip,            //!< target SCIP to be filled
   const char*           probname,           //!< name for subproblem
   int                   ncomponents,        //!< number of components
   int*                  components,         //!< array to store the components
   int                   compidx             //!< index of component to treat
   )
{
   assert( scip != nullptr );
   assert( subscip != nullptr );
   assert( components != nullptr );
   assert( compidx < ncomponents );

   // get original problem data
   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);

   // create data of subproblem
   SCIP_PROBDATA* subprobdata;
   SCIP_CALL( BACSinitProblem(subscip, &subprobdata) );

   // construct graph
   Graph* SG = new Graph();

   Vertex* vertexmap;
   SCIP_CALL( SCIPallocBufferArray(scip, &vertexmap, probdata->n) );

   // create nodes
   const Graph G = *probdata->G;
   subprobdata->unweighted = TRUE;
   for (size_t i = 0; i < probdata->n; ++i)
   {
      assert( -1 <= components[i] && components[i] < ncomponents );
      if ( components[i] == compidx )
      {
         Vertex v = boost::add_vertex(*SG);

         // treat weights
         SCIP_Real weight = boost::get(vertex_weight_t(), G, i);
         if ( ! SCIPisEQ(scip, weight, 1.0) )
            subprobdata->unweighted = FALSE;
         boost::put(vertex_weight_t(), *SG, v, weight);

         vertexmap[i] = v;
      }
      else
         vertexmap[i] = Graph::null_vertex();
   }

   // create edges
   EdgeIterator eit, eend;
   for (boost::tie(eit, eend) = boost::edges(G); eit != eend; ++eit)
   {
      Vertex s = boost::source(*eit, G);
      Vertex t = boost::target(*eit, G);

      if ( components[s] != compidx || components[t] != compidx )
         continue;

      Vertex subs = vertexmap[s];
      Vertex subt = vertexmap[t];
      assert( subs != Graph::null_vertex() && subt != Graph::null_vertex() );

      // create edge
      (void) boost::add_edge(subs, subt, *SG);
   }
   SCIPfreeBufferArray(scip, &vertexmap);

   subprobdata->G = SG;
   subprobdata->n = boost::num_vertices(*SG);
   subprobdata->m = boost::num_edges(*SG);

   SCIP_CALL( BACSinitDatastructures(subscip, subprobdata) );
   SCIP_CALL( BACSsetupProblem(subscip, probname, subprobdata) );

   // update subscip depth
   SCIPsetSubscipDepth(subscip, SCIPgetSubscipDepth(scip) + 1);

   return SCIP_OKAY;
}


//! solve a given sub-SCIP up to the given limits
static
SCIP_RETCODE solveSubscip(
   SCIP*                 scip,               //!< main SCIP
   SCIP*                 subscip,            //!< sub-SCIP to solve
   SCIP_Longint          nodelimit           //!< node limit
   )
{
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   SCIP_Bool avoidmemout;

   assert( scip != nullptr );
   assert( subscip != nullptr );

   // set time limit
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if ( ! SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);

   // Substract the memory already used by the main SCIP and the estimated memory usage of external software.
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if ( ! SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip) / 1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip) / 1048576.0;
   }

   // check if mem limit needs to be avoided
   SCIP_CALL( SCIPgetBoolParam(scip, "misc/avoidmemout", &avoidmemout) );

   /* Abort if no time is left or not enough memory (we don't abort in this case if misc_avoidmemout == TRUE)
    * to create a copy of SCIP, including external memory usage. */
   if ( avoidmemout && memorylimit <= 0.0 )
   {
      SCIPdebugMsg(scip, "--> not solved (not enough memory left)\n");
      return SCIP_OKAY;
   }
   else if ( timelimit <= 0.0 )
   {
      SCIPdebugMsg(scip, "--> not solved (not enough time left)\n");
      return SCIP_OKAY;
   }

   /* SCIPcopyLimits will set wrong time limits since it does not take into account time spent already in the sub-SCIP;
    * nevertheless, we call it to set the memory limit and unset all other limits, if set in the main SCIP. */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );

   // set time and memory limit for the subproblem
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );

   // set node limit
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );

   // solve the subproblem
   SCIP_CALL( SCIPsolve(subscip) );

#ifdef SCIP_MORE_DEBUG
   SCIP_CALL( SCIPprintBestSol(subscip, nullptr, FALSE) );
   SCIP_CALL( SCIPprintStatistics(subscip, nullptr) );
#endif

   return SCIP_OKAY;
}


//! solve a connected component during presolving and evaluate the result
static
SCIP_RETCODE solveAndEvalSubscip(
   SCIP*                 scip,               //!< SCIP main data structure
   SCIP_PROPDATA*        propdata,           //!< the propagator data
   SCIP*                 subscip,            //!< sub-SCIP to be solved
   SCIP_SOL*             neworigsol,         //!< solution of original problem to be constructed
   SCIP_VAR**            vars,               //!< array of variables copied to this component
   SCIP_VAR**            subvars,            //!< array of sub-SCIP variables corresponding to the vars array
   int                   nsubvars,           //!< number of variables in sub-SCIP
   int*                  nfixedvars,         //!< pointer to store the number of fixed variables
   int*                  nchgbds,            //!< pointer to store the number of bound tightenings
   SCIP_Bool*            cutoff,             //!< whether a cutoff was detected
   SCIP_Bool*            solved              //!< pointer to store if the problem was solved to optimality
   )
{
   assert( scip != nullptr );
   assert( propdata != nullptr );
   assert( subscip != nullptr );
   assert( vars != nullptr );
   assert( nfixedvars != nullptr );
   assert( nchgbds != nullptr );
   assert( cutoff != nullptr );
   assert( solved != nullptr );

   *solved = FALSE;
   *cutoff = FALSE;

   // try to transfer solution to subscip
   SCIP_SOL* bestsol = SCIPgetBestSol(scip);
   if ( bestsol != nullptr )
   {
      SCIP_SOL* subsol;
      SCIP_CALL( SCIPcreateSol(subscip, &subsol, nullptr) );

      // set solution values
      for (int v = 0; v < nsubvars; ++v)
      {
         assert( subvars[v] != nullptr );
         SCIP_CALL( SCIPsetSolVal(subscip, subsol, subvars[v], SCIPgetSolVal(scip, bestsol, vars[v])) );
      }

      assert( SCIPgetStage(subscip) == SCIP_STAGE_PROBLEM );

      SCIP_Bool feasible;
      SCIP_CALL( SCIPcheckSolOrig(subscip, subsol, &feasible, FALSE, FALSE) );
      if ( feasible )
      {
         SCIPdebugMsg(scip, "Transfered solution of original problem to subscip with primal bound %.9g.\n",
            SCIPgetSolOrigObj(subscip, subsol));
         SCIP_CALL( SCIPaddSolFree(subscip, &subsol, &feasible) );
      }
      else
      {
         SCIP_CALL( SCIPfreeSol(subscip, &subsol) );
      }
   }

   // solve subproblem; set node limit if we run in the tree
   if ( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      SCIP_CALL( solveSubscip(scip, subscip, -1) );
   }
   else
   {
      assert( SCIPgetStage(scip) == SCIP_STAGE_SOLVING );
      SCIP_CALL( solveSubscip(scip, subscip, propdata->nodelimit) );
   }

   // use result of sub-SCIP
   if ( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
   {
      SCIP_SOL* sol = SCIPgetBestSol(subscip);
      assert( sol != nullptr );

      SCIP_Bool feasible;
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPcheckSolOrig(subscip, sol, &feasible, TRUE, TRUE) );
#else
      SCIP_CALL( SCIPcheckSolOrig(subscip, sol, &feasible, FALSE, FALSE) );
#endif

      SCIPdebugMsg(scip, "--> solved to optimality: time = %.2f, value = %g, solution is%s feasible\n", SCIPgetSolvingTime(subscip), SCIPgetPrimalbound(subscip), feasible ? "" : " not");

      if ( feasible )
      {
         // fix variables to the values of the optimal solution of the subproblem
         for (int i = 0; i < nsubvars; ++i)
         {
            SCIP_VAR* var = vars[i];
            assert( var != nullptr );
            assert( SCIPvarIsBinary(var) );
            assert( SCIPvarGetLbLocal(var) < 0.5 && SCIPvarGetUbLocal(var) > 0.5 );

            SCIP_VAR* subvar = subvars[i];
            assert( subvar != nullptr );
            assert( SCIPvarIsBinary(subvar) );

            // get solution value from optimal solution of the sub-SCIP
            SCIP_Real subval = SCIPgetSolVal(subscip, sol, subvar);
            assert( SCIPisFeasIntegral(subscip, subval) );

            SCIP_Bool infeasible;
            SCIP_Bool fixed;
            if ( subval > 0.5 )
            {
               if ( neworigsol != nullptr )
               {
                  SCIP_CALL( SCIPsetSolVal(scip, neworigsol, var, 1.0) );
               }
               SCIP_CALL( SCIPfixVar(scip, var, 1.0, &infeasible, &fixed) );
            }
            else
            {
               if ( neworigsol != nullptr )
               {
                  SCIP_CALL( SCIPsetSolVal(scip, neworigsol, var, 0.0) );
               }
               SCIP_CALL( SCIPfixVar(scip, var, 0.0, &infeasible, &fixed) );
            }
            SCIPvarMarkDeleteGlobalStructures(var);
            assert( ! infeasible );
            assert( fixed );

            ++(*nfixedvars);
         }
         *solved = TRUE;
      }
   }
   else if ( SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE )
   {
      SCIPdebugMsg(scip, "--> subproblem is infeasible\n");
      *cutoff = TRUE;
   }
   else if ( SCIPgetStatus(subscip) == SCIP_STATUS_UNBOUNDED || SCIPgetStatus(subscip) == SCIP_STATUS_INFORUNBD )
   {
      SCIPerrorMessage("Subproblem is unbounded - this should not happen.\n");
      return SCIP_INVALIDDATA;
   }
   else
   {
      SCIPdebugMsg(scip, "--> solving interrupted (status = %d, time = %.2f)\n", SCIPgetStatus(subscip), SCIPgetSolvingTime(subscip));

      SCIP_Bool infeasible;
      SCIP_Bool tightened;
      int ntightened = 0;

      SCIP_SOL* sol = SCIPgetBestSol(subscip);

      // Transfer global bounds from subscip to the original problem: this is valid since any reduction on the
      // subproblem is valid for the original problem, even dual reductions (note that the optimal value on the subscip
      // is a lower bound on the value of the original, because we can set the variables of the other components to 0).
      for (int i = 0; i < nsubvars; ++i)
      {
         assert( subvars[i] != nullptr );

         SCIP_VAR* var = vars[i];
         assert( var != nullptr );
         assert( SCIPvarIsBinary(var) );
         assert( SCIPvarGetLbLocal(var) < 0.5 && SCIPvarGetUbLocal(var) > 0.5 );

         SCIP_VAR* subvar = subvars[i];
         assert( subvar != nullptr );
         assert( SCIPvarIsBinary(subvar) );

         // get solution value from optimal solution of the sub-SCIP
         if ( neworigsol != nullptr && sol != nullptr )
         {
            SCIP_Real subval = SCIPgetSolVal(subscip, sol, subvar);
            assert( SCIPisFeasIntegral(subscip, subval) );
            if ( subval > 0.5 )
            {
               SCIP_CALL( SCIPsetSolVal(scip, neworigsol, var, 1.0) );
            }
            else
            {
               SCIP_CALL( SCIPsetSolVal(scip, neworigsol, var, 0.0) );
            }
         }

         // tighten bounds
         SCIP_CALL( SCIPtightenVarLb(scip, vars[i], SCIPvarGetLbGlobal(subvar), FALSE, &infeasible, &tightened) );
         assert( ! infeasible );
         if ( tightened )
            ++ntightened;

         SCIP_CALL( SCIPtightenVarUb(scip, vars[i], SCIPvarGetUbGlobal(subvar), FALSE, &infeasible, &tightened) );
         assert( ! infeasible );
         if ( tightened )
            ++ntightened;
      }

      *nchgbds += ntightened;

      SCIPdebugMsg(scip, "--> tightened %d bounds of variables due to global bounds in the sub-SCIP\n", ntightened);

      // TODO: Check whether other information can be transferred, e.g., conflict constraints and branching information.
   }

   return SCIP_OKAY;
}


//! determine connected components through BFS, ignoring nodes whose variable is fixed
static
SCIP_RETCODE BFS(
   SCIP*                 scip,               //!< SCIP main data structure
   const Graph&          G,                  //!< graph
   SCIP_VAR**            vars,               //!< corresponding variables
   int*                  components,         //!< array to store the component index for each variable
   int&                  ncomponents         //!< pointer to store the number of components
   )
{
   assert( scip != nullptr );
   assert( components != nullptr );
   assert( vars != nullptr );

   ncomponents = 0;

   // compute components by BFS
   size_t n = boost::num_vertices(G);
   Vertex* Q; // queue
   SCIP_CALL( SCIPallocBufferArray(scip, &Q, (int) n) );

   for (size_t i = 0; i < n; ++i)
      components[i] = -1;

   size_t vertexidx = 0;
   while ( vertexidx < n )
   {
      if ( SCIPvarGetLbLocal(vars[vertexidx]) < 0.5 && SCIPvarGetUbLocal(vars[vertexidx]) > 0.5 )
      {
         // mark vertex
         components[vertexidx] = ncomponents;

         // init queue
         Q[0] = (Vertex) vertexidx;
         size_t startq = 0;
         size_t endq = 0;

         // loop through Q until empty
         while ( startq <= endq )
         {
            Vertex v = Q[startq++];
            assert( startq <= n );

            AdjacencyIterator ait, aend;
            for (boost::tie(ait, aend) = boost::adjacent_vertices(v, G); ait != aend; ++ait)
            {
               Vertex w = *ait;

               // do not connect through fixed vertices
               if ( SCIPvarGetLbLocal(vars[w]) < 0.5 && SCIPvarGetUbLocal(vars[w]) > 0.5 )
               {
                  if ( components[w] < 0 )
                  {
                     components[w] = ncomponents;
                     Q[++endq] = w;
                     assert( endq <= n );
                  }
               }
            }
         }
         ++ncomponents;
      }

      // find next vertex that is not yet covered
      do
      {
         ++vertexidx;
      }
      while ( vertexidx < n && components[vertexidx] >= 0 );
   }
   SCIPfreeBufferArray(scip, &Q);

   return SCIP_OKAY;
}


//! find nontrivial components in the problem
static
SCIP_RETCODE findComponents(
   SCIP*                 scip,               //!< SCIP main data structure
   SCIP_PROBDATA*        probdata,           //!< problem data structure
   int                   maxsize,            //!< maximal size of component (-1 if unbounded)
   int*                  components,         //!< array to store the component index for each variable
   SCIP_VAR**            compvars,           //!< array to store variables of each component, should have enough size for all variables
   int*                  begcompvars,        //!< start points of components in compvars array
   int&                  ncomponents,        //!< pointer to store the number of components
   int&                  ntoolarge,          //!< number of components that are too large
   int&                  nsmall              //!< number of small componets (size <= 3)
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( components != nullptr );
   assert( compvars != nullptr );
   assert( begcompvars != nullptr );
   assert( probdata->n > 0 );

   ncomponents = 0;
   ntoolarge = 0;
   nsmall = 0;

   if ( probdata->n == 0 )
      return SCIP_OKAY;

   // determine connected components by BFS, ignoring fixed nodes
   SCIP_CALL( BFS(scip, *probdata->G, probdata->vars, components, ncomponents) );

   // exit if we only have at most one component
   if ( ncomponents <= 1 )
      return SCIP_OKAY;

   // set up permutation and reverse permutation of components for sorting
   int* perm;
   int* revperm;
   int* ncompvars;
   SCIP_CALL( SCIPallocBufferArray(scip, &ncompvars, ncomponents) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, ncomponents) );
   SCIP_CALL( SCIPallocBufferArray(scip, &revperm, ncomponents) );

   // initialize size of components
   for (int c = 0; c < ncomponents; ++c)
   {
      ncompvars[c] = 0;
      perm[c] = c;
   }

   // determine sizes of the components
   for (size_t i = 0; i < probdata->n; ++i)
   {
      if ( components[i] >= 0 )
      {
         assert( components[i] < ncomponents );
         ++ncompvars[components[i]];
      }
   }

   // sort components according to size
   SCIPsortIntInt(ncompvars, perm, ncomponents);

   // compute reverse mapping
   for (int c = 0; c < ncomponents; ++c)
   {
      assert( 0 <= perm[c] && perm[c] < ncomponents );
      revperm[perm[c]] = c;
   }

   // replace component numbers
   for (size_t i = 0; i < probdata->n; ++i)
   {
      if ( components[i] >= 0 )
         components[i] = revperm[components[i]];
   }

   SCIPfreeBufferArray(scip, &revperm);
   SCIPfreeBufferArray(scip, &perm);

#ifndef NDEBUG
   for (int c = 0; c < ncomponents; ++c)
   {
      int size = 0;
      for (size_t i = 0; i < probdata->n; ++i)
      {
         if ( components[i] == c )
            ++size;
      }
      assert( size == ncompvars[c] );
   }
#endif

   // determine start of variables
   begcompvars[0] = 0;
   for (int c = 0; c < ncomponents; ++c)
   {
      assert( ncompvars[c] > 0 );

      // count number of compents that are too large
      if ( maxsize >= 0 && ncompvars[c] > maxsize )
         ++ntoolarge;

      // count number of small components
      if ( ncompvars[c] <= 3 )
         ++nsmall;

      begcompvars[c + 1] = begcompvars[c] + ncompvars[c];
      ncompvars[c] = 0;
   }
   assert( begcompvars[ncomponents] <= (int) probdata->n );

   // fill in variables
   for (size_t i = 0; i < probdata->n; ++i)
   {
      int c = components[i];
      if ( c >= 0 )
      {
         assert( c < ncomponents );
         int idx = begcompvars[c] + ncompvars[c];
         assert( 0 <= idx && idx < begcompvars[c + 1] );
         compvars[idx] = probdata->vars[i];
         ++ncompvars[c];
      }
   }
   SCIPfreeBufferArray(scip, &ncompvars);

   return SCIP_OKAY;
}


//! handle small components explicitly
static
SCIP_RETCODE handleSmallComponents(
   SCIP*                 scip,               //!< SCIP main data structure
   SCIP_PROBDATA*        probdata,           //!< problem data structure
   SCIP_SOL*             neworigsol,         //!< solution to be constructed
   int                   comp,               //!< component number for output
   SCIP_VAR**            compvars,           //!< array to store variables of given component
   int                   size,               //!< size of component
   int&                  nfixedsingle,       //!< number of variables fixed in single components
   int&                  nfixeddouble,       //!< number of variables fixed in components of size 2
   int&                  nfixedtriple,       //!< number of variables fixed in components of size 3
   int*                  nfixedvars,         //!< number of variables that were fixed
   SCIP_Bool*            cutoff,             //!< whether a cutoff was detected
   SCIP_Bool&            treatedcomp         //!< whether component has been treated
   )
{
   treatedcomp = FALSE;

   // treat components of size 1
   if ( size <= 1 )
   {
      // make sure single variables in a component are fixed
      SCIP_VAR* var = compvars[0];
      assert( SCIPvarGetLbLocal(var) < 0.5 && SCIPvarGetUbLocal(var) > 0.5 );

      // check objective to be sure - dual presolving should have fixed variables with negative objective
      if ( SCIPvarGetObj(var) >= 0.0 )
      {
         SCIPdebugMsg(scip, "Fixing single variable <%s> in component %d to 1.\n", SCIPvarGetName(var), comp);
         SCIP_Bool fixed;
         SCIP_CALL( SCIPfixVar(scip, var, 1.0, cutoff, &fixed) );
         SCIPvarMarkDeleteGlobalStructures(var);
         if ( neworigsol != nullptr )
         {
            SCIP_CALL( SCIPsetSolVal(scip, neworigsol, var, 1.0) );
         }
         assert( ! *cutoff );
         assert( fixed );
         ++nfixedsingle;
         ++(*nfixedvars);
         treatedcomp = TRUE;
      }
   }
   else if ( size == 2 )
   {
      // we have a component of two nodes that are connected by an edge -> fix one var to 1 and the other to 0
      SCIP_VAR* var1 = compvars[0];
      SCIP_VAR* var2 = compvars[1];
      assert( SCIPvarGetLbLocal(var1) < 0.5 && SCIPvarGetUbLocal(var1) > 0.5 );
      assert( SCIPvarGetLbLocal(var2) < 0.5 && SCIPvarGetUbLocal(var2) > 0.5 );

      // check objective to be sure - dual presolving should have fixed variables with negative objective
      if ( SCIPvarGetObj(var1) >= 0.0 && SCIPvarGetObj(var2) >= 0.0 )
      {
         if ( SCIPvarGetObj(var1) < SCIPvarGetObj(var2) )
            SCIPswapPointers((void**) &var1, (void**) &var2);

         SCIPdebugMsg(scip, "Fixing variable <%s> in component %d to 1, <%s> to 0.\n", SCIPvarGetName(var1), comp, SCIPvarGetName(var2));
         SCIP_Bool fixed;
         SCIP_CALL( SCIPfixVar(scip, var1, 1.0, cutoff, &fixed) );
         assert( ! *cutoff );
         if ( fixed )
            ++(*nfixedvars);
         SCIP_CALL( SCIPfixVar(scip, var2, 0.0, cutoff, &fixed) );
         assert( ! *cutoff );
         if ( fixed )
            ++(*nfixedvars);
         SCIPvarMarkDeleteGlobalStructures(var1);
         SCIPvarMarkDeleteGlobalStructures(var2);
         if ( neworigsol != nullptr )
         {
            SCIP_CALL( SCIPsetSolVal(scip, neworigsol, var1, 1.0) );
            SCIP_CALL( SCIPsetSolVal(scip, neworigsol, var2, 0.0) );
         }
         nfixeddouble += 2;
         treatedcomp = TRUE;
      }
   }
   else if ( size == 3 )
   {
      // we have a component of three nodes that are connected
      SCIP_VAR* var1 = compvars[0];
      SCIP_VAR* var2 = compvars[1];
      SCIP_VAR* var3 = compvars[2];
      assert( SCIPvarGetLbLocal(var1) < 0.5 && SCIPvarGetUbLocal(var1) > 0.5 );
      assert( SCIPvarGetLbLocal(var2) < 0.5 && SCIPvarGetUbLocal(var2) > 0.5 );
      assert( SCIPvarGetLbLocal(var3) < 0.5 && SCIPvarGetUbLocal(var3) > 0.5 );
      SCIP_Real obj1 = SCIPvarGetObj(var1);
      SCIP_Real obj2 = SCIPvarGetObj(var2);
      SCIP_Real obj3 = SCIPvarGetObj(var3);

      // check objective to be sure - dual presolving should have fixed variables with negative objective
      if ( obj1 >= 0.0 && obj2 >= 0.0 && obj3 >= 0.0 )
      {
         // get data about graph
         Vertex v1 = (Vertex) SCIPvarGetData(var1);
         Vertex v2 = (Vertex) SCIPvarGetData(var2);
         Vertex v3 = (Vertex) SCIPvarGetData(var3);
         SCIP_Bool edgev1v2 = FALSE;
         SCIP_Bool edgev1v3 = FALSE;
         SCIP_Bool edgev2v3 = FALSE;

         // check whether variable var1 is incident to var 2 and var3
         AdjacencyIterator ait, aend;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v1, *probdata->G); ait != aend; ++ait)
         {
            if ( *ait == v2 )
            {
               edgev1v2 = TRUE;
               if ( edgev1v3 )
                  break;
            }

            if ( *ait == v3 )
            {
               edgev1v3 = TRUE;
               if ( edgev1v2 )
                  break;
            }
         }

         // check whether variable var2 is incident to var3
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v2, *probdata->G); ait != aend; ++ait)
         {
            if ( *ait == v3 )
            {
               edgev2v3 = TRUE;
               break;
            }
         }
         assert( edgev1v2 + edgev1v3 + edgev2v3 >= 2 );

         // determine new values
         SCIP_Real newval1 = 0.0;
         SCIP_Real newval2 = 0.0;
         SCIP_Real newval3 = 0.0;
         if ( ! edgev1v2 )
         {
            if ( obj1 + obj2 >= obj3 )
            {
               newval1 = 1.0;
               newval2 = 1.0;
            }
            else
               newval3 = 1.0;
         }
         else if ( ! edgev1v3 )
         {
            if ( obj1 + obj3 >= obj2 )
            {
               newval1 = 1.0;
               newval3 = 1.0;
            }
            else
               newval2 = 1.0;
         }
         else if ( ! edgev2v3 )
         {
            if ( obj2 + obj3 >= obj1 )
            {
               newval2 = 1.0;
               newval3 = 1.0;
            }
            else
               newval1 = 1.0;
         }
         else
         {
            // in this case, the three nodes form a clique and we pick a variable with largest weight
            if ( obj1 >= obj2 && obj1 >= obj3 )
               newval1 = 1.0;
            else if ( obj2 >= obj1 && obj2 >= obj3 )
               newval2 = 1.0;
            else
            {
               assert( obj3 >= obj1 && obj3 >= obj2 );
               newval3 = 1.0;
            }
         }

         SCIPdebugMsg(scip, "Component %d: Fixing var <%s> = %f, <%s> = %f, <%s> = %f.\n", comp, SCIPvarGetName(var1), newval1, SCIPvarGetName(var2), newval2, SCIPvarGetName(var3), newval3);

         SCIP_Bool fixed;
         SCIP_CALL( SCIPfixVar(scip, var1, newval1, cutoff, &fixed) );
         assert( ! *cutoff );
         assert( fixed );
         ++(*nfixedvars);
         SCIP_CALL( SCIPfixVar(scip, var2, newval2, cutoff, &fixed) );
         assert( ! *cutoff );
         assert( fixed );
         ++(*nfixedvars);
         SCIP_CALL( SCIPfixVar(scip, var3, newval3, cutoff, &fixed) );
         assert( ! *cutoff );
         assert( fixed );
         ++(*nfixedvars);

         SCIPvarMarkDeleteGlobalStructures(var1);
         SCIPvarMarkDeleteGlobalStructures(var2);
         SCIPvarMarkDeleteGlobalStructures(var3);
         if ( neworigsol != nullptr )
         {
            SCIP_CALL( SCIPsetSolVal(scip, neworigsol, var1, newval1) );
            SCIP_CALL( SCIPsetSolVal(scip, neworigsol, var2, newval2) );
            SCIP_CALL( SCIPsetSolVal(scip, neworigsol, var3, newval3) );
         }
         nfixedtriple += 3;
         treatedcomp = TRUE;
      }
   }

   return SCIP_OKAY;
}


//! run fixing loop
static
SCIP_RETCODE handleComponents(
   SCIP*                 scip,               //!< SCIP main data structure
   SCIP_PROPDATA*        propdata,           //!< propagator data
   SCIP_PROBDATA*        probdata,           //!< problem data structure
   SCIP_Bool             inpresol,           //!< whether we are called from the prsolving routine
   int                   maxsize,            //!< maximal size of components to solve (-1 if unbounded)
   int*                  nfixedvars,         //!< number of variables that were fixed
   int*                  nchgbds,            //!< number of bounds that were changed
   SCIP_Bool*            cutoff              //!< whether a cutoff was detected
   )
{
   assert( scip != nullptr );
   assert( propdata != nullptr );
   assert( probdata != nullptr );
   assert( nfixedvars != nullptr );
   assert( nchgbds != nullptr );
   assert( cutoff != nullptr );

   *cutoff = FALSE;

   // define needed variables
   SCIP_VAR** compvars;
   int* components;
   int* begcompvars;

   // allocate memory for sorted components
   SCIP_CALL( SCIPallocBufferArray(scip, &components, probdata->n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &compvars, probdata->n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &begcompvars, probdata->n + 1) );

   // find components
   int ncomponents;
   int ntoolarge;
   int nsmall;
   SCIP_CALL( findComponents(scip, probdata, maxsize, components, compvars, begcompvars, ncomponents, ntoolarge, nsmall) );

   // exit if at most one nontrival and not too large component is left
   if ( ncomponents - ntoolarge <= 1 )
   {
      SCIPfreeBufferArray(scip, &begcompvars);
      SCIPfreeBufferArray(scip, &compvars);
      SCIPfreeBufferArray(scip, &components);
      return SCIP_OKAY;
   }

   // give output in presolving
   if ( inpresol )
   {
      if ( propdata->onlysmall )
      {
         if ( nsmall > 0 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "   (%.1fs) Treat %d small components ...\n",
               SCIPgetSolvingTime(scip), nsmall);
         }
         else
         {
            SCIPfreeBufferArray(scip, &begcompvars);
            SCIPfreeBufferArray(scip, &compvars);
            SCIPfreeBufferArray(scip, &components);
            return SCIP_OKAY;
         }
      }
      else
      {
         if ( maxsize >= 0 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "   (%.1fs) Found %d connected components of size in [%d,%d] in the graph ...\n",
               SCIPgetSolvingTime(scip), ncomponents - ntoolarge, 2, maxsize);
         }
         else
         {
            assert( ntoolarge == 0 );
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "   (%.1fs) Found %d connected components in the graph ...\n",
               SCIPgetSolvingTime(scip), ncomponents);
         }
      }
   }

   // create solution of the original problem to be set-up by the components
   SCIP_SOL* neworigsol = nullptr;
   if ( SCIPgetBestSol(scip) != nullptr )
   {
      SCIP_CALL( SCIPcreateSolCopy(scip, &neworigsol, SCIPgetBestSol(scip)) );
   }

   // loop over all components
   int nsolved = 0;
   int compidx = 0;
   int comp;
   int nfixedsingle = 0;
   int nfixeddouble = 0;
   int nfixedtriple = 0;
   for (comp = 0; comp < ncomponents && ! SCIPisStopped(scip); ++comp)
   {
      assert( 0 <= begcompvars[comp] && begcompvars[comp] < (int) probdata->n );

      int size = begcompvars[comp + 1] - begcompvars[comp];
      assert( size > 0 );

      // exit for components that are too large
      if ( maxsize >= 0 && size > maxsize )
         break;

      // treat components of size 1 to 3
      if ( size <= 3 )
      {
         SCIP_Bool treatedcomp;
         SCIP_CALL( handleSmallComponents(scip, probdata, neworigsol, comp, &compvars[begcompvars[comp]], size, nfixedsingle, nfixeddouble, nfixedtriple, nfixedvars, cutoff, treatedcomp) );
         assert( ! *cutoff );
         assert( treatedcomp );
         ++nsolved; // treat component as solved

         if ( ! inpresol )
         {
            assert( size < propdata->ncompsizes );
            ++propdata->compsizes[size];
         }

         continue;
      }

      // exit if we only want to treat small components
      if ( propdata->onlysmall )
         break;

      // if there is only one component left, let's solve this in the main SCIP
      if ( nsolved == ncomponents - 1 )
         break;

      // create subscip
      SCIP* subscip;
      SCIP_CALL( createSubscip(scip, &subscip, propdata->minsizesym, size) );

      if ( subscip == nullptr )
         break;

      // add event handler if we are in presolving
      if ( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
      {
         SCIP_CALL( SCIPincludeEventHdlrSubscip(subscip, scip, compidx) );
      }

      // create subproblem
      char name[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp_%d", SCIPgetProbName(scip), comp);
      SCIP_CALL( createStableSetSubproblem(scip, subscip, name, ncomponents, components, comp) );

      SCIP_PROBDATA* subprobdata;
      subprobdata = SCIPgetProbData(subscip);
      assert( subprobdata != nullptr );

      // define for later use
      SCIP_VAR** subvars = subprobdata->vars;
      SCIP_VAR** vars = &compvars[begcompvars[comp]];
      int nsubvars = size;

      // set up debug solution
#if ( defined(WITH_DEBUG_SOLUTION) && SCIP_VERSION_MAJOR >= 10 )
      if ( SCIPdebugSolIsEnabled(scip) )
      {
         SCIP_SOL* debugsol;
         SCIP_Real val;

         SCIPdebugSolEnable(scip);
         SCIP_CALL( SCIPdebugGetSol(scip, &debugsol) );

         // set solution values in the debug solution if it is available
         if ( debugsol != nullptr )
         {
            // clear solution of sub-SCIP, since the variables do not match to the main SCIP
            SCIP_CALL( SCIPdebugClearSol(subscip) );
            SCIPdebugSolEnable(subscip);

            SCIP_Real obj = 0.0;
            for (int i = 0; i < size; ++i)
            {
               assert( subvars[i] != nullptr );
               SCIP_CALL( SCIPdebugGetSolVal(scip, vars[i], &val) );
               obj += SCIPvarGetObj(vars[i]) * val;
               SCIP_CALL( SCIPdebugAddSolVal(subscip, subvars[i], val) );
            }
            SCIPdebugMsg(scip, "Value of debug solution on component = %g.\n", obj);
         }
      }
#endif

      // solve the subproblem and evaluate the result, i.e., apply fixings of variables
      SCIP_Bool solved;
      int noldchgbds = *nchgbds;
      int noldfixedvars = *nfixedvars;

      SCIP_CALL( solveAndEvalSubscip(scip, propdata, subscip, neworigsol, vars, subvars, nsubvars, nfixedvars, nchgbds, cutoff, &solved) );

      if ( solved )
      {
         ++nsolved;

         if ( ! inpresol )
         {
            assert( size < propdata->ncompsizes );
            ++propdata->compsizes[size];
         }
      }

      if ( inpresol )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr,
            "   (%.1fs) Component %d: size %d, fixed %d variables, tightened %d bounds, #nodes = %" SCIP_LONGINT_FORMAT ".\n",
            SCIPgetSolvingTime(scip), compidx++, size, *nfixedvars - noldfixedvars, *nchgbds - noldchgbds, SCIPgetNNodes(subscip));
      }

      // transfer statistics from subscip to scip
      SCIP_CALL( SCIPtransferStatistics(subscip, scip) );

      SCIP_CALL( SCIPfree(&subscip) );

      if ( *cutoff )
         break;  // if the component is infeasible, this holds for the complete problem as well
      else if ( ! solved )
         break;  // stop if we could not solve the last problem
   }

   if ( inpresol && nfixedsingle + nfixeddouble + nfixedtriple > 0 )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "   (%.1fs) Fixed %d, %d, %d variables in components of size 1, 2, 3.\n", SCIPgetSolvingTime(scip), nfixedsingle, nfixeddouble, nfixedtriple);

   if ( nsolved == ncomponents - 1 && SCIPgetDepth(scip) <= 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "   (%.1fs) Continue with solving remaining component of size %d ...\n",
         SCIPgetSolvingTime(scip), begcompvars[nsolved + 1] - begcompvars[nsolved]);
   }

   if ( neworigsol != nullptr )
   {
      SCIP_Bool stored;
      SCIP_CALL( SCIPtrySolFree(scip, &neworigsol, FALSE, FALSE, FALSE, FALSE, FALSE, &stored) );
   }

   SCIPfreeBufferArray(scip, &begcompvars);
   SCIPfreeBufferArray(scip, &compvars);
   SCIPfreeBufferArray(scip, &components);

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */

//! copy method for propagator plugins (called when SCIP copies plugins)
static
SCIP_DECL_PROPCOPY(propCopyComp)
{  /*lint --e{715}*/
   assert( scip != nullptr );
   assert( prop != nullptr );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   // call inclusion method of propagator
   SCIP_CALL( BACSincludePropComp(scip) );

   return SCIP_OKAY;
}


//! destructor of propagator to free user data (called when SCIP is exiting)
static
SCIP_DECL_PROPFREE(propFreeComp)
{  /*lint --e{715}*/
   assert( scip != nullptr );
   assert( prop != nullptr );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   /* free propagator data */
   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);

   assert( propdata != nullptr );

   if ( SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_MINIMAL )
   {
      unsigned int maxindcomp = 0;
      size_t numbcomp = 0;

      for (int i = 0; i < propdata->ncompsizes; ++i)
      {
         if ( propdata->compsizes[i] > 0 )
            maxindcomp = (unsigned int) i;

         numbcomp += propdata->compsizes[i];
      }

      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\nTotal number of components fixed in the tree:\t\t%8zu\n", numbcomp);

      if ( numbcomp > 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Component sizes:\n");
         for (unsigned int i = 1; i <= maxindcomp; ++i)
         {
            int w = (int) ceil(log10(MAX(i, propdata->compsizes[i])+0.1));
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "%*u ", w, i);
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\n");
         for (unsigned int i = 1; i <= maxindcomp; ++i)
         {
            int w = (int) ceil(log10(MAX(i, propdata->compsizes[i])+0.1));
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "%*lu ", w, propdata->compsizes[i]);
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\n");
      }
   }

   SCIPfreeBlockMemoryArray(scip, &propdata->compsizes, propdata->ncompsizes);
   SCIPfreeBlockMemory(scip, &propdata);
   SCIPpropSetData(prop, nullptr);

   return SCIP_OKAY;
}


//! solving process initialization method of propagator (called when branch and bound process is about to begin)
static
SCIP_DECL_PROPINITSOL(propInitsolComp)
{
   assert( prop != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   // if propagation should be performed in the tree, we currently turn off symmetry handling to avoid conflicts, unless we turn it off below anyways
   if ( SCIPpropGetFreq(prop) >= 0 && ! propdata->onlyfornosymmetry )
   {
      SCIP_RETCODE retcode = SCIPsetIntParam(scip, "propagating/symmetry/freq", -1);
      if ( retcode == SCIP_OKAY )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "Turning off symmetry handling to avoid conflicts with component propagation.\n");
      else
         assert( retcode == SCIP_PARAMETERUNKNOWN );
   }

   return SCIP_OKAY;
}


//! presolving initialization method of propagator (called when presolving is about to begin)
static
SCIP_DECL_PROPINITPRE(propInitPre)
{  /*lint --e{715}*/
   assert( scip != nullptr );
   assert( prop != nullptr );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   if ( propdata == nullptr )
      return SCIP_OKAY;

   // this presolver only works if the original problem just contains one clique constraint
   if ( SCIPgetNOrigConss(scip) != 1 )
      propdata->enabled = FALSE;
   else
   {
      SCIP_CONS** conss = SCIPgetOrigConss(scip);
      assert( conss != nullptr );
      if ( strcmp(SCIPconsGetName(conss[0]), "clique") != 0 )
         propdata->enabled = FALSE;
   }

   if ( ! propdata->enabled )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "Turned off components presolver, because the problem does not contain exaclty one clique constraint.\n");

   return SCIP_OKAY;
}

//! execution method of propagator
static
SCIP_DECL_PROPEXEC(propExecComp)
{  /*lint --e{715}*/
   assert( scip != nullptr );
   assert( prop != nullptr );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );
   assert( SCIPgetNActivePricers(scip) == 0 );

   *result = SCIP_DIDNOTRUN;

   // do not run, if not all variables are explicitly known
   if ( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   // do not run in probing or in repropagation
   if ( SCIPinProbing(scip) || SCIPinRepropagation(scip) )
      return SCIP_OKAY;

   // avoid recursive call
   if ( SCIPgetSubscipDepth(scip) > 0 )
      return SCIP_OKAY;

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   if ( probdata == nullptr )
      return SCIP_OKAY;

   // we do not want to run, if there are no variables left
   if ( probdata->n == 0 )
      return SCIP_OKAY;

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   // do not run if disabled
   if ( ! propdata->enabled )
      return SCIP_OKAY;

   // do not try to detect independent components if the depth is too high
   if ( SCIPgetDepth(scip) > propdata->maxdepth )
      return SCIP_OKAY;

   // check for a reached timelimit
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   // the comp propagator performs dual reductions
   if ( ! SCIPallowStrongDualReds(scip) || ! SCIPallowWeakDualReds(scip) )
      return SCIP_OKAY;

   // possibly only run in the tree if there are no symmetries
   if ( propdata->onlyfornosymmetry && SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING && SCIPgetSymmetryNGenerators(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   // search and possibly handle components
   SCIP_Bool cutoff;
   int nfixedvars = 0;
   int nchgbds = 0;
   SCIP_CALL( handleComponents(scip, propdata, probdata, FALSE, propdata->maxsize, &nfixedvars, &nchgbds, &cutoff) );

   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( nfixedvars + nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


//! presolving method
static
SCIP_DECL_PROPPRESOL(propPresolComp)
{  /*lint --e{715}*/
   assert( scip != nullptr );
   assert( prop != nullptr );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );
   assert( result != nullptr );

   *result = SCIP_DIDNOTRUN;

   if ( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING || SCIPinProbing(scip) )
      return SCIP_OKAY;

   // do not run, if not all variables are explicitly known
   if ( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   // only call the components presolving, if presolving would be stopped otherwise
   if ( ! SCIPisPresolveFinished(scip) )
      return SCIP_OKAY;

   // the comp propagator performs dual reductions
   if ( ! SCIPallowStrongDualReds(scip) || ! SCIPallowWeakDualReds(scip) )
      return SCIP_OKAY;

   // avoid recursive call
   if ( SCIPgetSubscipDepth(scip) > 0 )
      return SCIP_OKAY;

   // check for a reached timelimit
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   // get problem data
   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   if ( probdata == nullptr )
      return SCIP_OKAY;

   // search and possibly handle components
   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   // do not run if disabled
   if ( ! propdata->enabled )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIP_Bool cutoff;
   int noldfixedvars = *nfixedvars;
   int noldchgbds = *nchgbds;
   // do not use maxsize but rely on node limit
   SCIP_CALL( handleComponents(scip, propdata, probdata, TRUE, -1, nfixedvars, nchgbds, &cutoff) );

   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( (*nfixedvars - noldfixedvars) + (*nchgbds - noldchgbds) > 0 )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


//! creates the comp propagator and includes it in SCIP
SCIP_RETCODE BACSincludePropComp(
   SCIP*                 scip                //!< SCIP data structure
   )
{
   SCIP_PROPDATA* propdata = nullptr;
   SCIP_PROP* prop;

   assert( scip != nullptr );

   // create comp data and allocate memory
   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );
   propdata->enabled = TRUE;

   // include propagator
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecComp, propdata) );
   assert( prop != nullptr );

   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolComp) );
   SCIP_CALL( SCIPsetPropInitpre(scip, prop, propInitPre) );
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyComp) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeComp) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolComp, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS,
         PROP_PRESOLTIMING) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxsize",
         "maximal size of components to consider in the tree (-1: no limit)",
         &propdata->maxsize, FALSE, DEFAULT_MAXSIZE, -1, INT_MAX, nullptr, nullptr) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/maxdepth",
         "maximum depth of a node to run components detection (-1: disable component detection during solving)",
         &propdata->maxdepth, FALSE, DEFAULT_MAXDEPTH, -1, INT_MAX, nullptr, nullptr) );

   SCIP_CALL( SCIPaddLongintParam(scip,
         "propagating/" PROP_NAME "/nodelimit",
         "limit on number of nodes in tree for solving subproblems during solving",
         &propdata->nodelimit, FALSE, DEFAULT_NODELIMIT, -1LL, SCIP_LONGINT_MAX, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/onlysmall",
         "only treat small components",
         &propdata->onlysmall, FALSE, DEFAULT_ONLYSMALL, nullptr, nullptr) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/minsizesym",
         "minimal size of component to include symmetry handling for",
         &propdata->minsizesym, FALSE, DEFAULT_MINSIZE_SYM, -1, INT_MAX, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/onlyfornosymmetry",
         "whether the propagator should only be run in the tree if there are no symmetries",
         &propdata->onlyfornosymmetry, FALSE, DEFAULT_ONLYFORNOSYMMETRY, nullptr, nullptr) );

   if ( propdata->maxsize >= 0 )
   {
      propdata->ncompsizes = propdata->maxsize + 1;
   }
   else
   {
      SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
      propdata->ncompsizes = (int) probdata->n;
   }

   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &propdata->compsizes, propdata->ncompsizes) );

   return SCIP_OKAY;
}
