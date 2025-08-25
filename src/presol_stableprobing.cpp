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

/**@file   presol_stableprobing.cpp
 * @brief  presolver to perform probing
 * @author Marc Pfetsch
 *
 * This code is based on scip/presol_probing.c.
 */

#include "presol_stableprobing.h"
#include "graph.h"
#include "struct_probdata.h"

#define PRESOL_NAME            "stableprobing"
#define PRESOL_DESC            "probing presolver specialized to stable set problems"
#define PRESOL_PRIORITY         -500000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS              0 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */


/*
 * Default parameter settings
 */

#define DEFAULT_MAXRUNS              1  /**< maximal number of runs, probing participates in (-1: no limit) */
#define DEFAULT_PROPROUNDS          -1  /**< maximal number of propagation rounds in probing subproblems */
#define DEFAULT_MAXTOTALUSELESS   1000  /**< maximal number of succ. probings without fixings, bound changes,
                                         *   and implications, until probing is aborted (0: don't abort) */

//! presolver data
struct SCIP_PresolData
{
   int                   maxruns;            /**< maximal number of runs, probing participates in (-1: no limit) */
   int                   proprounds;         /**< maximal number of propagation rounds in probing subproblems */
   int                   maxtotaluseless;    /**< maximal number of succ. probings without fixings, bound changes,
                                              *   and implications, until probing is aborted (0: don't abort) */
};



/** applies and evaluates probing of a single variable in the given direction */
static
SCIP_RETCODE applyProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   proprounds,         /**< maximal number of propagation rounds in probing subproblems */
   SCIP_VAR**            vars,               /**< variables */
   size_t                n,                  /**< number of nodes */
   SCIP_VAR*             probvar,            /**< variable to probe on */
   SCIP_Bool             probingdir,         /**< value to fix probing variable to */
   SCIP_Real*            impllbs,            /**< array to store lower bounds after applying implications and cliques */
   SCIP_Real*            implubs,            /**< array to store upper bounds after applying implications and cliques */
   SCIP_Real*            proplbs,            /**< array to store lower bounds after full propagation */
   SCIP_Real*            propubs,            /**< array to store upper bounds after full propagation */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   )
{  /*lint --e{818}*/
   assert( scip != nullptr );
   assert( vars != nullptr );
   assert( probvar != nullptr );
   assert( impllbs != nullptr );
   assert( implubs != nullptr );
   assert( proplbs != nullptr );
   assert( propubs != nullptr );
   assert( cutoff != nullptr );
   assert( SCIPvarGetType(probvar) == SCIP_VARTYPE_BINARY );
   assert( SCIPvarGetLbLocal(probvar) < 0.5 );
   assert( SCIPvarGetUbLocal(probvar) > 0.5 );

   SCIPdebugMsg(scip, "applying probing on variable <%s> == %u (impls=%d/%d)\n",
      SCIPvarGetName(probvar), probingdir, SCIPvarGetNImpls(probvar, FALSE), SCIPvarGetNImpls(probvar, TRUE) );

   // start probing mode
   SCIP_CALL( SCIPstartProbing(scip) );

   // fix variable
   if ( probingdir == FALSE )
   {
      SCIP_CALL( SCIPchgVarUbProbing(scip, probvar, 0.0) );
   }
   else
   {
      SCIP_CALL( SCIPchgVarLbProbing(scip, probvar, 1.0) );
   }

   // apply propagation of implication graph and clique table
   SCIP_CALL( SCIPpropagateProbingImplications(scip, cutoff) );
   if ( !(*cutoff) )
   {
      for (size_t i = 0; i < n; ++i)
      {
         impllbs[i] = SCIPvarGetLbLocal(vars[i]);
         implubs[i] = SCIPvarGetUbLocal(vars[i]);
      }
   }

   // apply propagation
   if ( !(*cutoff) )
   {
      SCIP_CALL( SCIPpropagateProbing(scip, proprounds, cutoff, nullptr) );
   }

   // evaluate propagation
   if ( !(*cutoff) )
   {
      for (size_t i = 0; i < n; ++i)
      {
         proplbs[i] = SCIPvarGetLbLocal(vars[i]);
         propubs[i] = SCIPvarGetUbLocal(vars[i]);
      }
   }

   // exit probing mode
   SCIP_CALL( SCIPendProbing(scip) );

   return SCIP_OKAY;
}



/*
 * Callback methods of presolver
 */

/** execution method of presolver */
SCIP_DECL_PRESOLEXEC(presolExecStableprobing)
{  /*lint --e{715}*/
   assert( scip != nullptr );
   assert( presol != nullptr );
   assert( result != nullptr );
   assert( nfixedvars != nullptr );
   assert( nchgbds != nullptr );

   *result = SCIP_DIDNOTRUN;

   SCIP_PRESOLDATA* presoldata = SCIPpresolGetData(presol);

   // check, if probing should be applied in the current run
   if ( presoldata->maxruns >= 0 && SCIPgetNRuns(scip) > presoldata->maxruns )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "executing probing (used %.1f sec)\n", SCIPpresolGetTime(presol));

   *result = SCIP_DIDNOTFIND;

   // get data
   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   // get sizes and variables
   size_t n = probdata->n;
   SCIP_VAR** vars = probdata->vars;
   assert( vars != nullptr );

   int oldnfixedvars = *nfixedvars;
   int oldnchgbds = *nchgbds;
   int nimplications = 0;
   SCIP_Bool cutoff = FALSE;
   int cnt = 0;

   // get storage for computed bounds
   SCIP_Real* zeroimpllbs;
   SCIP_Real* zeroimplubs;
   SCIP_Real* zeroproplbs;
   SCIP_Real* zeropropubs;
   SCIP_Real* oneimpllbs;
   SCIP_Real* oneimplubs;
   SCIP_Real* oneproplbs;
   SCIP_Real* onepropubs;

   SCIP_CALL( SCIPallocBufferArray(scip, &zeroimpllbs, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zeroimplubs, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zeroproplbs, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zeropropubs, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oneimpllbs, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oneimplubs, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oneproplbs, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &onepropubs, n) );

   // check all variables
   unsigned int maxtotaluseless = (unsigned int) (presoldata->maxtotaluseless > 0 ? presoldata->maxtotaluseless : INT_MAX);
   unsigned int ntotaluseless = 0;
   for (size_t i = 0; i < n && ! cutoff; ++i)
   {
      SCIP_VAR* probvar = vars[i];

      // ignore variables, that were fixed, aggregated, or deleted in prior probings
      if ( ! SCIPvarIsActive(probvar) || SCIPvarIsDeleted(probvar) || SCIPvarGetLbGlobal(probvar) > 0.5 || SCIPvarGetUbGlobal(probvar) < 0.5 )
         continue;

      ++cnt;
      ++ntotaluseless;

      if ( ntotaluseless >= maxtotaluseless )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr,
            "   (%.1fs) stableprobing aborted: %u/%u successive totally useless probings\n", SCIPgetSolvingTime(scip),
            ntotaluseless, maxtotaluseless);
         break;
      }

      // check whether probing should be aborted
      if ( SCIPisStopped(scip) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "   (%.1fs) stableprobing aborted: solving stopped\n", SCIPgetSolvingTime(scip));
         break;
      }

      // display probing status
      if ( cnt % 100 == 0 )
      {
         SCIP_VERBLEVEL verblevel = (cnt % 1000 == 0 ? SCIP_VERBLEVEL_HIGH : SCIP_VERBLEVEL_FULL);
         SCIPverbMessage(scip, verblevel, nullptr,
            "   (%.1fs) stableprobing: %d/%lu - %d fixings, %d implications, %d bound changes\n",
            SCIPgetSolvingTime(scip), cnt, n, *nfixedvars, nimplications, *nchgbds);
      }

      // apply probing for fixing the variable to 0
      SCIP_Bool localcutoff = FALSE;
      SCIP_CALL( applyProbing(scip, presoldata->proprounds, vars, n, probvar, FALSE,
            zeroimpllbs, zeroimplubs, zeroproplbs, zeropropubs, &localcutoff) );

      if ( localcutoff )
      {
         SCIP_Bool fixed;

         // the variable can be fixed to 1
         SCIPdebugMsg(scip, "fixing stableprobing variable <%s> to 1.0\n", SCIPvarGetName(probvar));
         SCIP_CALL( SCIPfixVar(scip, probvar, 1.0, &cutoff, &fixed) );
         assert( fixed );
         ++(*nfixedvars);
         ntotaluseless = 0;
         continue; // don't try upwards direction, because the variable is already fixed
      }

      // apply probing for fixing the variable to 1
      SCIP_CALL( applyProbing(scip, presoldata->proprounds, vars, n, probvar, TRUE,
            oneimpllbs, oneimplubs, oneproplbs, onepropubs, &localcutoff) );

      if ( localcutoff )
      {
         SCIP_Bool fixed;

         // the variable can be fixed to 0
         SCIPdebugMsg(scip, "fixing stableprobing variable <%s> to 0.0\n", SCIPvarGetName(probvar));
         SCIP_CALL( SCIPfixVar(scip, probvar, 0.0, &cutoff, &fixed) );
         assert( fixed );
         ++(*nfixedvars);
         ntotaluseless = 0;
         continue; // don't analyze probing deductions, because the variable is already fixed
      }

      // analyze probing deductions
      for (size_t k = 0; k < n; ++k)
      {
         SCIP_VAR* var = vars[k];

         // skip already fixed variables
         if ( ! SCIPvarIsActive(var) || SCIPvarIsDeleted(var) || SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5 )
            continue;

         // new bounds of the variable are the union of the propagated bounds of the zero and one case
         SCIP_Real newlb = MIN(zeroproplbs[k], oneproplbs[k]);
         SCIP_Real newub = MAX(zeropropubs[k], onepropubs[k]);

         // check for fixed variables
         if ( SCIPisFeasEQ(scip, newlb, newub) )
         {
            SCIP_Real fixval;

            // in both probings, variable is deduced to 0: fix variable to 0, otherwise to 1
            if ( newlb < 0.5 )
               fixval = 0.0;
            else
               fixval = 1.0;

            SCIP_Bool fixed = FALSE;
            SCIP_CALL( SCIPfixVar(scip, var, fixval, &cutoff, &fixed) );
            if ( fixed )
            {
               SCIPdebugMsg(scip, "fixed variable <%s> to %d due to stableprobing on <%s>\n",
                  SCIPvarGetName(var), (int) fixval, SCIPvarGetName(probvar));
               ++(*nfixedvars);
            }
            continue;
         }

         int nboundchanges;
         if ( zeropropubs[k] < 0.5 && zeroimplubs[k] > 0.5 )
         {
            // insert implication: var == 0  =>  var == 0
            // SCIPdebugMsg(scip, "found implication <%s> == 0  =>  <%s> == 0\n", SCIPvarGetName(probvar), SCIPvarGetName(var));
            SCIP_CALL( SCIPaddVarImplication(scip, probvar, FALSE, var, SCIP_BOUNDTYPE_UPPER, 0.0, &cutoff, &nboundchanges) );
            ++nimplications;
            (*nchgbds) += nboundchanges;
            ntotaluseless = 0;
         }
         else if ( zeroproplbs[k] > 0.5 && zeroimpllbs[k] < 0.5 )
         {
            // insert implication: probvar == 0  =>  var == 1
            // SCIPdebugMsg(scip, "found implication <%s> == 0  =>  <%s> == 1\n", SCIPvarGetName(probvar), SCIPvarGetName(var));
            SCIP_CALL( SCIPaddVarImplication(scip, probvar, FALSE, var, SCIP_BOUNDTYPE_LOWER, 1.0, &cutoff, &nboundchanges) );
            ++nimplications;
            (*nchgbds) += nboundchanges;
            ntotaluseless = 0;
         }
         else if ( onepropubs[k] < 0.5 && oneimplubs[k] > 0.5 )
         {
            // insert implication: probvar == 1  =>  var == 0
            // SCIPdebugMsg(scip, "found implication <%s> == 1  =>  <%s> == 0\n", SCIPvarGetName(probvar), SCIPvarGetName(var));
            SCIP_CALL( SCIPaddVarImplication(scip, probvar, TRUE, var, SCIP_BOUNDTYPE_UPPER, 0.0, &cutoff, &nboundchanges) );
            ++nimplications;
            (*nchgbds) += nboundchanges;
            ntotaluseless = 0;
         }
         else if ( oneproplbs[k] > 0.5 && oneimpllbs[k] < 0.5 )
         {
            // insert implication: probvar == 1  =>  var == 1
            // SCIPdebugMsg(scip, "found implication <%s> == 1  =>  <%s> == 1\n", SCIPvarGetName(probvar), SCIPvarGetName(var));
            SCIP_CALL( SCIPaddVarImplication(scip, probvar, TRUE, var, SCIP_BOUNDTYPE_LOWER, 1.0, &cutoff, &nboundchanges) );
            ++nimplications;
            (*nchgbds) += nboundchanges;
            ntotaluseless = 0;
         }
      }

      if ( SCIPisStopped(scip) || ntotaluseless >= maxtotaluseless )
         break;
   }

   // free temporary memory
   SCIPfreeBufferArray(scip, &onepropubs);
   SCIPfreeBufferArray(scip, &oneproplbs);
   SCIPfreeBufferArray(scip, &oneimplubs);
   SCIPfreeBufferArray(scip, &oneimpllbs);
   SCIPfreeBufferArray(scip, &zeropropubs);
   SCIPfreeBufferArray(scip, &zeroproplbs);
   SCIPfreeBufferArray(scip, &zeroimplubs);
   SCIPfreeBufferArray(scip, &zeroimpllbs);

   // adjust result code
   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( *nfixedvars > oldnfixedvars || *nchgbds > oldnchgbds || nimplications > 0 )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeStableprobing)
{
   assert( scip != nullptr );
   assert( presol != nullptr );

   SCIP_PRESOLDATA* presoldata = SCIPpresolGetData(presol);
   assert( presoldata != nullptr );
   SCIPfreeBlockMemory(scip, &presoldata);

   return SCIP_OKAY;
}



/*
 * presolver specific interface methods
 */

//! creates the neighborhoods presolver and includes it in SCIP
SCIP_RETCODE BACSincludePresolStableprobing(
   SCIP*                 scip                //! SCIP data structure
   )
{
   SCIP_PRESOLDATA* presoldata = nullptr;
   SCIP_PRESOL* presol = nullptr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   // create neighborhoods presolver data
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecStableprobing, presoldata) );
   assert( presol != nullptr );

   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeStableprobing) );

   // add probing presolver parameters
   SCIP_CALL_ABORT( SCIPaddIntParam(scip,
         "presolving/stableprobing/maxruns",
         "maximal number of runs, probing participates in (-1: no limit)",
         &presoldata->maxruns, FALSE, DEFAULT_MAXRUNS, -1, INT_MAX, nullptr, nullptr) );
   SCIP_CALL_ABORT( SCIPaddIntParam(scip,
         "presolving/stableprobing/proprounds",
         "maximal number of propagation rounds in probing subproblems (-1: no limit, 0: auto)",
         &presoldata->proprounds, TRUE, DEFAULT_PROPROUNDS, -1, INT_MAX, nullptr, nullptr) );
   SCIP_CALL_ABORT( SCIPaddIntParam(scip,
         "presolving/stableprobing/maxtotaluseless",
         "maximal number of succ. probings without fixings, bound changes, and implications, until probing is aborted (0: don't abort)",
         &presoldata->maxtotaluseless, TRUE, DEFAULT_MAXTOTALUSELESS, 0, INT_MAX, nullptr, nullptr) );

   return SCIP_OKAY;
}
