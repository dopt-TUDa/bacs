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

/**@file   transfer_stats.c
 * @brief  code to transfer statistics from one SCIP to another one
 * @author Marc Pfetsch
 */

#include <transfer_stats.h>

#include "scip/clock.h"
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/struct_presol.h"
#include "scip/struct_prop.h"
#include "scip/struct_cons.h"
#include "scip/conflict.h"
#include "scip/struct_conflict.h"
#include "scip/struct_sepa.h"
#include "scip/struct_cutpool.h"
#include "scip/struct_cutsel.h"
#include "scip/struct_pricer.h"
#include "scip/struct_pricestore.h"
#include "scip/pricestore.h"
#include "scip/struct_branch.h"
#include "scip/struct_heur.h"
#include "scip/struct_stat.h"
#include "scip/scip_compr.h"
#include "scip/struct_compr.h"
#include "scip/struct_nlp.h"
#include "scip/scip_relax.h"
#include "scip/struct_relax.h"


/** transfer presolver data */
static
void SCIPtransferPresolStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIP_PRESOL** sourcepresols;
   int nsourcepresols;
   int i;

   sourcepresols = SCIPgetPresols(sourcescip);
   nsourcepresols = SCIPgetNPresols(sourcescip);

   for (i = 0; i < nsourcepresols; ++i)
   {
      SCIP_PRESOL* sourcepresol;
      SCIP_PRESOL* targetpresol;

      assert( sourcepresols != NULL );
      sourcepresol = sourcepresols[i];
      assert( sourcepresol != NULL );

      targetpresol = SCIPfindPresol(targetscip, SCIPpresolGetName(sourcepresol));
      if ( targetpresol != NULL )
      {
         SCIPclockSetTime(targetpresol->presolclock, SCIPpresolGetTime(sourcepresol) + SCIPpresolGetTime(targetpresol));
         SCIPclockSetTime(targetpresol->setuptime, SCIPpresolGetSetupTime(sourcepresol) + SCIPpresolGetSetupTime(targetpresol));
         targetpresol->ncalls += SCIPpresolGetNCalls(sourcepresol);
         targetpresol->nfixedvars += SCIPpresolGetNFixedVars(sourcepresol);
         targetpresol->naggrvars += SCIPpresolGetNAggrVars(sourcepresol);
         targetpresol->nchgvartypes += SCIPpresolGetNChgVarTypes(sourcepresol);
         targetpresol->nchgbds += SCIPpresolGetNChgBds(sourcepresol);
         targetpresol->naddholes += SCIPpresolGetNAddHoles(sourcepresol);
         targetpresol->ndelconss += SCIPpresolGetNDelConss(sourcepresol);
         targetpresol->naddconss += SCIPpresolGetNAddConss(sourcepresol);
         targetpresol->nchgsides += SCIPpresolGetNChgSides(sourcepresol);
         targetpresol->nchgcoefs += SCIPpresolGetNChgCoefs(sourcepresol);
      }
   }

   targetscip->stat->nrootintfixings += sourcescip->stat->nrootintfixings;
   targetscip->stat->nrootboundchgs += sourcescip->stat->nrootboundchgs;
}


/** transfer propagator data */
static
void SCIPtransferPropStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIP_PROP** sourceprops;
   int nsourceprops;
   int i;

   sourceprops = SCIPgetProps(sourcescip);
   nsourceprops = SCIPgetNProps(sourcescip);

   for (i = 0; i < nsourceprops; ++i)
   {
      SCIP_PROP* sourceprop;
      SCIP_PROP* targetprop;

      assert( sourceprops != NULL );
      sourceprop = sourceprops[i];
      assert( sourceprop != NULL );

      targetprop = SCIPfindProp(targetscip, SCIPpropGetName(sourceprop));

      if ( targetprop != NULL )
      {
         if( SCIPpropDoesPresolve(sourceprop) )
         {
            SCIPclockSetTime(targetprop->presoltime, SCIPpropGetPresolTime(sourceprop) + SCIPpropGetPresolTime(targetprop));

            targetprop->npresolcalls += SCIPpropGetNPresolCalls(sourceprop);
            targetprop->nfixedvars += SCIPpropGetNFixedVars(sourceprop);
            targetprop->naggrvars += SCIPpropGetNAggrVars(sourceprop);
            targetprop->nchgvartypes += SCIPpropGetNChgVarTypes(sourceprop);
            targetprop->nchgbds += SCIPpropGetNChgBds(sourceprop);
            targetprop->naddholes += SCIPpropGetNAddHoles(sourceprop);
            targetprop->ndelconss += SCIPpropGetNDelConss(sourceprop);
            targetprop->naddconss += SCIPpropGetNAddConss(sourceprop);
            targetprop->nchgsides += SCIPpropGetNChgSides(sourceprop);
            targetprop->nchgcoefs += SCIPpropGetNChgCoefs(sourceprop);
         }

         targetprop->ncalls += SCIPpropGetNCalls(sourceprop);
         targetprop->nrespropcalls += SCIPpropGetNRespropCalls(sourceprop);
         targetprop->ncutoffs += SCIPpropGetNCutoffs(sourceprop);
         targetprop->ndomredsfound += SCIPpropGetNDomredsFound(sourceprop);

         SCIPclockSetTime(targetprop->setuptime, SCIPpropGetSetupTime(sourceprop) + SCIPpropGetSetupTime(targetprop));
         SCIPclockSetTime(targetprop->proptime, SCIPpropGetTime(sourceprop) + SCIPpropGetTime(targetprop));
         SCIPclockSetTime(targetprop->resproptime, SCIPpropGetRespropTime(sourceprop) + SCIPpropGetRespropTime(targetprop));
         SCIPclockSetTime(targetprop->sbproptime, SCIPpropGetStrongBranchPropTime(sourceprop) + SCIPpropGetStrongBranchPropTime(targetprop));
      }
   }
}


/** transfer constraint handler data */
static
void SCIPtransferConshdlrStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIP_CONSHDLR** sourceconshdlrs;
   int nsourceconshdlrs;
   int i;

   sourceconshdlrs = SCIPgetConshdlrs(sourcescip);
   nsourceconshdlrs = SCIPgetNConshdlrs(sourcescip);

   for (i = 0; i < nsourceconshdlrs; ++i)
   {
      SCIP_CONSHDLR* sourceconshdlr;
      SCIP_CONSHDLR* targetconshdlr;

assert( sourceconshdlrs != NULL );
      sourceconshdlr = sourceconshdlrs[i];
      assert( sourceconshdlr != NULL );

      targetconshdlr = SCIPfindConshdlr(targetscip, SCIPconshdlrGetName(sourceconshdlr));

      if ( targetconshdlr != NULL )
      {
         if ( SCIPconshdlrDoesPresolve(targetconshdlr) || ! SCIPconshdlrNeedsCons(targetconshdlr) )
         {
            SCIPclockSetTime(targetconshdlr->presoltime, SCIPconshdlrGetPresolTime(sourceconshdlr) + SCIPconshdlrGetPresolTime(targetconshdlr));

            targetconshdlr->npresolcalls += SCIPconshdlrGetNPresolCalls(sourceconshdlr);
            targetconshdlr->nfixedvars += SCIPconshdlrGetNFixedVars(sourceconshdlr);
            targetconshdlr->naggrvars += SCIPconshdlrGetNAggrVars(sourceconshdlr);
            targetconshdlr->nchgvartypes += SCIPconshdlrGetNChgVarTypes(sourceconshdlr);
            targetconshdlr->nchgbds += SCIPconshdlrGetNChgBds(sourceconshdlr);
            targetconshdlr->naddholes += SCIPconshdlrGetNAddHoles(sourceconshdlr);
            targetconshdlr->ndelconss += SCIPconshdlrGetNDelConss(sourceconshdlr);
            targetconshdlr->naddconss += SCIPconshdlrGetNAddConss(sourceconshdlr);
            targetconshdlr->nchgsides += SCIPconshdlrGetNChgSides(sourceconshdlr);
            targetconshdlr->nchgcoefs += SCIPconshdlrGetNChgCoefs(sourceconshdlr);
         }

         SCIPclockSetTime(targetconshdlr->setuptime, SCIPconshdlrGetSetupTime(sourceconshdlr) + SCIPconshdlrGetSetupTime(targetconshdlr));
         SCIPclockSetTime(targetconshdlr->sepatime, SCIPconshdlrGetSepaTime(sourceconshdlr) + SCIPconshdlrGetSepaTime(targetconshdlr));
         SCIPclockSetTime(targetconshdlr->proptime, SCIPconshdlrGetPropTime(sourceconshdlr) + SCIPconshdlrGetPropTime(targetconshdlr));
         SCIPclockSetTime(targetconshdlr->enfolptime, SCIPconshdlrGetEnfoLPTime(sourceconshdlr) + SCIPconshdlrGetEnfoLPTime(targetconshdlr));
         SCIPclockSetTime(targetconshdlr->enfopstime, SCIPconshdlrGetEnfoPSTime(sourceconshdlr) + SCIPconshdlrGetEnfoPSTime(targetconshdlr));
         SCIPclockSetTime(targetconshdlr->enforelaxtime, SCIPconshdlrGetEnfoRelaxTime(sourceconshdlr) + SCIPconshdlrGetEnfoRelaxTime(targetconshdlr));
         SCIPclockSetTime(targetconshdlr->checktime, SCIPconshdlrGetCheckTime(sourceconshdlr) + SCIPconshdlrGetCheckTime(targetconshdlr));
         SCIPclockSetTime(targetconshdlr->resproptime, SCIPconshdlrGetRespropTime(sourceconshdlr) + SCIPconshdlrGetRespropTime(targetconshdlr));
         SCIPclockSetTime(targetconshdlr->sbproptime, SCIPconshdlrGetStrongBranchPropTime(sourceconshdlr) + SCIPconshdlrGetStrongBranchPropTime(targetconshdlr));

         targetconshdlr->nsepacalls += SCIPconshdlrGetNSepaCalls(sourceconshdlr);
         targetconshdlr->npropcalls += SCIPconshdlrGetNPropCalls(sourceconshdlr);
         targetconshdlr->nenfolpcalls += SCIPconshdlrGetNEnfoLPCalls(sourceconshdlr);
         targetconshdlr->nenforelaxcalls += SCIPconshdlrGetNEnfoRelaxCalls(sourceconshdlr);
         targetconshdlr->nenfopscalls += SCIPconshdlrGetNEnfoPSCalls(sourceconshdlr);
         targetconshdlr->ncheckcalls += SCIPconshdlrGetNCheckCalls(sourceconshdlr);
         targetconshdlr->nrespropcalls += SCIPconshdlrGetNRespropCalls(sourceconshdlr);
         targetconshdlr->ncutoffs += SCIPconshdlrGetNCutoffs(sourceconshdlr);
         targetconshdlr->ndomredsfound += SCIPconshdlrGetNDomredsFound(sourceconshdlr);
         targetconshdlr->ncutsfound += SCIPconshdlrGetNCutsFound(sourceconshdlr);
         targetconshdlr->ncutsapplied += SCIPconshdlrGetNCutsApplied(sourceconshdlr);
         targetconshdlr->nconssfound += SCIPconshdlrGetNConssFound(sourceconshdlr);
         targetconshdlr->nchildren += SCIPconshdlrGetNChildren(sourceconshdlr);
      }
   }
}


/** transfer conflict data */
static
void SCIPtransferConflictStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIP_CONFLICT* sourceconflict;
   SCIP_CONFLICT* targetconflict;

   sourceconflict = sourcescip->conflict;
   targetconflict = targetscip->conflict;

   targetconflict->npropcalls += SCIPconflictGetNPropCalls(sourceconflict);
   targetconflict->npropsuccess += SCIPconflictGetNPropSuccess(sourceconflict);
   targetconflict->npropconfconss += SCIPconflictGetNPropConflictConss(sourceconflict);
   targetconflict->npropconfliterals += SCIPconflictGetNPropConflictLiterals(sourceconflict);
   targetconflict->npropreconvconss += SCIPconflictGetNPropReconvergenceConss(sourceconflict);
   targetconflict->npropreconvliterals += SCIPconflictGetNPropReconvergenceLiterals(sourceconflict);

   SCIPclockSetTime(targetconflict->propanalyzetime, SCIPconflictGetPropTime(sourceconflict) + SCIPconflictGetPropTime(targetconflict));
   SCIPclockSetTime(targetconflict->inflpanalyzetime, SCIPconflictGetInfeasibleLPTime(sourceconflict) + SCIPconflictGetInfeasibleLPTime(targetconflict));

   targetconflict->ninflpcalls += SCIPconflictGetNInfeasibleLPCalls(sourceconflict);
   targetconflict->ninflpsuccess += SCIPconflictGetNInfeasibleLPSuccess(sourceconflict);
   targetconflict->ninflpconfconss += SCIPconflictGetNInfeasibleLPConflictConss(sourceconflict);
   targetconflict->ninflpconfliterals += SCIPconflictGetNInfeasibleLPConflictLiterals(sourceconflict);
   targetconflict->ninflpreconvconss += SCIPconflictGetNInfeasibleLPReconvergenceConss(sourceconflict);
   targetconflict->ninflpreconvliterals += SCIPconflictGetNInfeasibleLPReconvergenceLiterals(sourceconflict);

   targetconflict->ndualproofsinfsuccess += SCIPconflictGetNDualproofsInfSuccess(sourceconflict);
   targetconflict->dualproofsinfnnonzeros += SCIPconflictGetNDualproofsInfNonzeros(sourceconflict);
   targetconflict->ninflpiterations += SCIPconflictGetNInfeasibleLPIterations(sourceconflict);

   SCIPclockSetTime(targetconflict->boundlpanalyzetime, SCIPconflictGetBoundexceedingLPTime(sourceconflict) + SCIPconflictGetBoundexceedingLPTime(targetconflict));

   targetconflict->nboundlpcalls += SCIPconflictGetNBoundexceedingLPCalls(sourceconflict);
   targetconflict->nboundlpsuccess += SCIPconflictGetNBoundexceedingLPSuccess(sourceconflict);
   targetconflict->nboundlpconfconss += SCIPconflictGetNBoundexceedingLPConflictConss(sourceconflict);
   targetconflict->nboundlpconfliterals += SCIPconflictGetNBoundexceedingLPConflictLiterals(sourceconflict);
   targetconflict->nboundlpreconvconss += SCIPconflictGetNBoundexceedingLPReconvergenceConss(sourceconflict);
   targetconflict->nboundlpreconvliterals += SCIPconflictGetNBoundexceedingLPReconvergenceLiterals(sourceconflict);

   targetconflict->ndualproofsbndsuccess += SCIPconflictGetNDualproofsBndSuccess(sourceconflict);
   targetconflict->nboundlpiterations += SCIPconflictGetNBoundexceedingLPIterations(sourceconflict);

   SCIPclockSetTime(targetconflict->sbanalyzetime, SCIPconflictGetStrongbranchTime(sourceconflict) + SCIPconflictGetStrongbranchTime(targetconflict));

   targetconflict->nsbcalls += SCIPconflictGetNStrongbranchCalls(sourceconflict);
   targetconflict->nsbsuccess += SCIPconflictGetNStrongbranchSuccess(sourceconflict);
   targetconflict->nsbconfconss += SCIPconflictGetNStrongbranchConflictConss(sourceconflict);
   targetconflict->nsbconfliterals += SCIPconflictGetNStrongbranchConflictLiterals(sourceconflict);
   targetconflict->nsbreconvconss += SCIPconflictGetNStrongbranchReconvergenceConss(sourceconflict);
   targetconflict->nsbreconvliterals += SCIPconflictGetNStrongbranchReconvergenceLiterals(sourceconflict);
   targetconflict->nsbiterations += SCIPconflictGetNStrongbranchIterations(sourceconflict);

   SCIPclockSetTime(targetconflict->pseudoanalyzetime, SCIPconflictGetPseudoTime(sourceconflict) + SCIPconflictGetPseudoTime(targetconflict));

   targetconflict->npseudocalls += SCIPconflictGetNPseudoCalls(sourceconflict);
   targetconflict->npseudosuccess += SCIPconflictGetNPseudoSuccess(sourceconflict);
   targetconflict->npseudoconfconss += SCIPconflictGetNPseudoConflictConss(sourceconflict);
   targetconflict->npseudoconfliterals += SCIPconflictGetNPseudoConflictLiterals(sourceconflict);
   targetconflict->npseudoreconvconss += SCIPconflictGetNPseudoReconvergenceConss(sourceconflict);
   targetconflict->npseudoreconvliterals += SCIPconflictGetNPseudoReconvergenceLiterals(sourceconflict);

   SCIPclockSetTime(targetconflict->dIBclock, SCIPconflictGetGlobalApplTime(sourceconflict) + SCIPconflictGetGlobalApplTime(targetconflict));

   targetconflict->nglbchgbds += SCIPconflictGetNGlobalChgBds(sourceconflict);
   targetconflict->nappliedglbconss += SCIPconflictGetNAppliedGlobalConss(sourceconflict);
   targetconflict->nappliedglbliterals += SCIPconflictGetNAppliedGlobalLiterals(sourceconflict);
   targetconflict->ndualproofsinfglobal += SCIPconflictGetNDualproofsInfGlobal(sourceconflict);
   targetconflict->ndualproofsbndglobal += SCIPconflictGetNDualproofsBndGlobal(sourceconflict);

   targetconflict->nlocchgbds += SCIPconflictGetNLocalChgBds(sourceconflict);
   targetconflict->nappliedlocconss += SCIPconflictGetNAppliedLocalConss(sourceconflict);
   targetconflict->nappliedlocliterals += SCIPconflictGetNAppliedLocalLiterals(sourceconflict);
   targetconflict->ndualproofsinflocal += SCIPconflictGetNDualproofsInfLocal(sourceconflict);
   targetconflict->ndualproofsbndlocal += SCIPconflictGetNDualproofsBndLocal(sourceconflict);
}


/** transfer separator handler data */
static
void SCIPtransferSeparatorStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIP_SEPA** sourcesepas;
   int nsourcesepas;
   int i;

   sourcesepas = SCIPgetSepas(sourcescip);
   nsourcesepas = SCIPgetNSepas(sourcescip);

   for (i = 0; i < nsourcesepas; ++i)
   {
      SCIP_SEPA* sourcesepa;
      SCIP_SEPA* targetsepa;

assert( sourcesepas != NULL );
      sourcesepa = sourcesepas[i];
      assert( sourcesepa != NULL );

      targetsepa = SCIPfindSepa(targetscip, SCIPsepaGetName(sourcesepa));

      if ( targetsepa != NULL )
      {
         /* only output data for separators without parent separator */
         if ( SCIPsepaGetParentsepa(targetsepa) == NULL )
         {
            SCIPclockSetTime(targetsepa->sepaclock, SCIPsepaGetTime(sourcesepa) + SCIPsepaGetTime(targetsepa));
            SCIPclockSetTime(targetsepa->setuptime, SCIPsepaGetSetupTime(sourcesepa) + SCIPsepaGetSetupTime(targetsepa));

            targetsepa->ncalls += SCIPsepaGetNCalls(sourcesepa);
            targetsepa->nrootcalls += SCIPsepaGetNRootCalls(sourcesepa);
            targetsepa->ncutoffs += SCIPsepaGetNCutoffs(sourcesepa);
            targetsepa->ndomredsfound += SCIPsepaGetNDomredsFound(sourcesepa);
            targetsepa->ncutsfound += SCIPsepaGetNCutsFound(sourcesepa);
            targetsepa->ncutsaddedviapool += SCIPsepaGetNCutsAddedViaPool(sourcesepa);
            targetsepa->ncutsaddeddirect += SCIPsepaGetNCutsAddedDirect(sourcesepa);
            targetsepa->ncutsappliedviapool += SCIPsepaGetNCutsAppliedViaPool(sourcesepa);
            targetsepa->ncutsapplieddirect += SCIPsepaGetNCutsAppliedDirect(sourcesepa);
            targetsepa->nconssfound += SCIPsepaGetNConssFound(sourcesepa);
         }
      }
   }

   if ( targetscip->cutpool != NULL )
   {
      SCIPclockSetTime(targetscip->cutpool->poolclock, SCIPcutpoolGetTime(sourcescip->cutpool) + SCIPcutpoolGetTime(targetscip->cutpool));
      targetscip->cutpool->ncalls += SCIPcutpoolGetNCalls(sourcescip->cutpool);
      targetscip->cutpool->nrootcalls += SCIPcutpoolGetNRootCalls(sourcescip->cutpool);
      targetscip->cutpool->ncutsfound += SCIPcutpoolGetNCutsFound(sourcescip->cutpool);
      targetscip->cutpool->ncutsadded += SCIPcutpoolGetNCutsAdded(sourcescip->cutpool);
      targetscip->cutpool->maxncuts += SCIPcutpoolGetMaxNCuts(sourcescip->cutpool);
   }
}


/** transfer cut selector data */
static
void SCIPtransferCutSelectorStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIP_CUTSEL** sourcecutsels;
   int nsourcecutsels;
   int i;

   sourcecutsels = SCIPgetCutsels(sourcescip);
   nsourcecutsels = SCIPgetNCutsels(sourcescip);

   for (i = 0; i < nsourcecutsels; ++i)
   {
      SCIP_CUTSEL* sourcecutsel;
      SCIP_CUTSEL* targetcutsel;

assert( sourcecutsels != NULL );
      sourcecutsel = sourcecutsels[i];
      assert( sourcecutsel != NULL );

      targetcutsel = SCIPfindCutsel(targetscip, SCIPcutselGetName(sourcecutsel));

      if ( targetcutsel != NULL )
      {
         SCIPclockSetTime(targetcutsel->cutseltime, SCIPcutselGetTime(sourcecutsel) + SCIPcutselGetTime(targetcutsel));
         SCIPclockSetTime(targetcutsel->setuptime, SCIPcutselGetSetupTime(sourcecutsel) + SCIPcutselGetSetupTime(targetcutsel));

         targetcutsel->ncalls += SCIPcutselGetNCalls(sourcecutsel);
         targetcutsel->nrootcalls += SCIPcutselGetNRootCalls(sourcecutsel);
         targetcutsel->nrootcutsselected += SCIPcutselGetNRootCuts(sourcecutsel);
         targetcutsel->nlocalcutsselected += SCIPcutselGetNLocalCuts(sourcecutsel);
         targetcutsel->nrootcutsforced += SCIPcutselGetNRootForcedCuts(sourcecutsel);
         targetcutsel->nlocalcutsforced += SCIPcutselGetNLocalForcedCuts(sourcecutsel);
         targetcutsel->nrootcutsfiltered += SCIPcutselGetNRootCutsFiltered(sourcecutsel);
         targetcutsel->nlocalcutsfiltered += SCIPcutselGetNLocalCutsFiltered(sourcecutsel);
      }
   }
}


/** transfer pricer data */
static
void SCIPtransferPricerStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIP_PRICER** sourcepricers;
   int nsourcepricers;
   int i;

   sourcepricers = SCIPgetPricers(sourcescip);
   nsourcepricers = SCIPgetNPricers(sourcescip);

   for (i = 0; i < nsourcepricers; ++i)
   {
      SCIP_PRICER* sourcepricer;
      SCIP_PRICER* targetpricer;

assert( sourcepricers != NULL );
      sourcepricer = sourcepricers[i];
      assert( sourcepricer != NULL );

      targetpricer = SCIPfindPricer(targetscip, SCIPpricerGetName(sourcepricer));

      if ( targetpricer != NULL )
      {
         SCIPclockSetTime(targetpricer->pricerclock, SCIPpricerGetTime(sourcepricer) + SCIPpricerGetTime(targetpricer));
         SCIPclockSetTime(targetpricer->setuptime, SCIPpricerGetSetupTime(sourcepricer) + SCIPpricerGetSetupTime(targetpricer));

         targetpricer->ncalls += SCIPpricerGetNCalls(sourcepricer);
         targetpricer->nvarsfound += SCIPpricerGetNVarsFound(sourcepricer);
      }
   }

   if ( targetscip->pricestore != NULL )
   {
      SCIPclockSetTime(targetscip->pricestore->probpricingtime, SCIPpricestoreGetProbPricingTime(sourcescip->pricestore) + SCIPpricestoreGetProbPricingTime(targetscip->pricestore));
      targetscip->pricestore->nprobpricings += SCIPpricestoreGetNProbPricings(sourcescip->pricestore);
      targetscip->pricestore->nprobvarsfound += SCIPpricestoreGetNProbvarsFound(sourcescip->pricestore);
   }
}


/** transfer branching rule data */
static
void SCIPtransferBranchruleStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIP_BRANCHRULE** sourcebranchrules;
   int nsourcebranchrules;
   int i;

   sourcebranchrules = SCIPgetBranchrules(sourcescip);
   nsourcebranchrules = SCIPgetNBranchrules(sourcescip);

   for (i = 0; i < nsourcebranchrules; ++i)
   {
      SCIP_BRANCHRULE* sourcebranchrule;
      SCIP_BRANCHRULE* targetbranchrule;

assert( sourcebranchrules != NULL );
      sourcebranchrule = sourcebranchrules[i];
      assert( sourcebranchrule != NULL );

      targetbranchrule = SCIPfindBranchrule(targetscip, SCIPbranchruleGetName(sourcebranchrule));

      if ( targetbranchrule != NULL )
      {
         SCIPclockSetTime(targetbranchrule->branchclock, SCIPbranchruleGetTime(sourcebranchrule) + SCIPbranchruleGetTime(targetbranchrule));
         SCIPclockSetTime(targetbranchrule->setuptime, SCIPbranchruleGetSetupTime(sourcebranchrule) + SCIPbranchruleGetSetupTime(targetbranchrule));

         targetbranchrule->nlpcalls += SCIPbranchruleGetNLPCalls(sourcebranchrule);
         targetbranchrule->nexterncalls += SCIPbranchruleGetNExternCalls(sourcebranchrule);
         targetbranchrule->npseudocalls += SCIPbranchruleGetNPseudoCalls(sourcebranchrule);
         targetbranchrule->ncutoffs += SCIPbranchruleGetNCutoffs(sourcebranchrule);
         targetbranchrule->ndomredsfound += SCIPbranchruleGetNDomredsFound(sourcebranchrule);
         targetbranchrule->ncutsfound += SCIPbranchruleGetNCutsFound(sourcebranchrule);
         targetbranchrule->nconssfound += SCIPbranchruleGetNConssFound(sourcebranchrule);
         targetbranchrule->nchildren += SCIPbranchruleGetNChildren(sourcebranchrule);
      }
   }
}


/** transfer heuristics data */
static
void SCIPtransferHeuristicsStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIP_HEUR** sourceheurs;
   int nsourceheurs;
   int i;

   sourceheurs = SCIPgetHeurs(sourcescip);
   nsourceheurs = SCIPgetNHeurs(sourcescip);

   for (i = 0; i < nsourceheurs; ++i)
   {
      SCIP_HEUR* sourceheur;
      SCIP_HEUR* targetheur;

assert( sourceheurs != NULL );
      sourceheur = sourceheurs[i];
      assert( sourceheur != NULL );

      targetheur = SCIPfindHeur(targetscip, SCIPheurGetName(sourceheur));

      if ( targetheur != NULL )
      {
         SCIPclockSetTime(targetheur->heurclock, SCIPheurGetTime(sourceheur) + SCIPheurGetTime(targetheur));
         SCIPclockSetTime(targetheur->setuptime, SCIPheurGetSetupTime(sourceheur) + SCIPheurGetSetupTime(targetheur));

         targetheur->ncalls += SCIPheurGetNCalls(sourceheur);
         targetheur->nsolsfound += SCIPheurGetNSolsFound(sourceheur);
         targetheur->nbestsolsfound += SCIPheurGetNBestSolsFound(sourceheur);
      }
   }

   SCIPclockSetTime(targetscip->stat->lpsoltime, SCIPclockGetTime(sourcescip->stat->lpsoltime) + SCIPclockGetTime(targetscip->stat->lpsoltime));
   targetscip->stat->nlpsolsfound += sourcescip->stat->nlpsolsfound;
   targetscip->stat->nlpbestsolsfound += sourcescip->stat->nlpbestsolsfound;

   SCIPclockSetTime(targetscip->stat->relaxsoltime, SCIPclockGetTime(sourcescip->stat->relaxsoltime) + SCIPclockGetTime(targetscip->stat->relaxsoltime));
   targetscip->stat->nrelaxsolsfound += sourcescip->stat->nrelaxsolsfound;
   targetscip->stat->nrelaxbestsolsfound += sourcescip->stat->nrelaxbestsolsfound;

   SCIPclockSetTime(targetscip->stat->pseudosoltime, SCIPclockGetTime(sourcescip->stat->pseudosoltime) + SCIPclockGetTime(targetscip->stat->pseudosoltime));
   targetscip->stat->npssolsfound += sourcescip->stat->npssolsfound;
   targetscip->stat->npsbestsolsfound += sourcescip->stat->npsbestsolsfound;

   SCIPclockSetTime(targetscip->stat->sbsoltime, SCIPclockGetTime(sourcescip->stat->sbsoltime) + SCIPclockGetTime(targetscip->stat->sbsoltime));
   targetscip->stat->nsbsolsfound += sourcescip->stat->nsbsolsfound;
   targetscip->stat->nsbbestsolsfound += sourcescip->stat->nsbbestsolsfound;

   targetscip->stat->nexternalsolsfound += sourcescip->stat->nexternalsolsfound;

   /* do not transfer diving heuristic information due to lazyness */
}


/** transfer compression data */
static
void SCIPtransferCompressionStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIP_COMPR** sourcecomprs;
   int nsourcecomprs;
   int i;

   sourcecomprs = SCIPgetComprs(sourcescip);
   nsourcecomprs = SCIPgetNCompr(sourcescip);

   for (i = 0; i < nsourcecomprs; ++i)
   {
      SCIP_COMPR* sourcecompr;
      SCIP_COMPR* targetcompr;

assert( sourcecomprs != NULL );
      sourcecompr = sourcecomprs[i];
      assert( sourcecompr != NULL );

      targetcompr = SCIPfindCompr(targetscip, SCIPcomprGetName(sourcecompr));

      if ( targetcompr != NULL )
      {
         SCIPclockSetTime(targetcompr->comprclock, SCIPcomprGetTime(sourcecompr) + SCIPcomprGetTime(targetcompr));
         SCIPclockSetTime(targetcompr->setuptime, SCIPcomprGetSetupTime(sourcecompr) + SCIPcomprGetSetupTime(targetcompr));

         targetcompr->ncalls += SCIPcomprGetNCalls(sourcecompr);
         targetcompr->nfound += SCIPcomprGetNFound(sourcecompr);
      }
   }
}


/** transfer LP data */
static
void SCIPtransferLPStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIPclockSetTime(targetscip->stat->primallptime, SCIPclockGetTime(sourcescip->stat->primallptime) + SCIPclockGetTime(targetscip->stat->primallptime));
   targetscip->stat->nprimallps += sourcescip->stat->nprimallps;
   targetscip->stat->nprimalzeroitlps += sourcescip->stat->nprimalzeroitlps;
   targetscip->stat->nprimallpiterations += sourcescip->stat->nprimallpiterations;
   targetscip->stat->primalzeroittime += sourcescip->stat->primalzeroittime;
   targetscip->stat->nprimalzeroitlps += sourcescip->stat->nprimalzeroitlps;

   SCIPclockSetTime(targetscip->stat->duallptime, SCIPclockGetTime(sourcescip->stat->duallptime) + SCIPclockGetTime(targetscip->stat->duallptime));
   targetscip->stat->nduallps += sourcescip->stat->nduallps;
   targetscip->stat->ndualzeroitlps += sourcescip->stat->ndualzeroitlps;
   targetscip->stat->nduallpiterations += sourcescip->stat->nduallpiterations;
   targetscip->stat->dualzeroittime += sourcescip->stat->dualzeroittime;
   targetscip->stat->ndualzeroitlps += sourcescip->stat->ndualzeroitlps;

   SCIPclockSetTime(targetscip->stat->lexduallptime, SCIPclockGetTime(sourcescip->stat->lexduallptime) + SCIPclockGetTime(targetscip->stat->lexduallptime));
   targetscip->stat->nlexduallps += sourcescip->stat->nlexduallps;
   targetscip->stat->nlexduallpiterations += sourcescip->stat->nlexduallpiterations;

   SCIPclockSetTime(targetscip->stat->barrierlptime, SCIPclockGetTime(sourcescip->stat->barrierlptime) + SCIPclockGetTime(targetscip->stat->barrierlptime));
   targetscip->stat->nbarrierlps += sourcescip->stat->nbarrierlps;
   targetscip->stat->nbarrierlpiterations += sourcescip->stat->nbarrierlpiterations;
   targetscip->stat->barrierzeroittime += sourcescip->stat->barrierzeroittime;
   targetscip->stat->nbarrierzeroitlps += sourcescip->stat->nbarrierzeroitlps;

   SCIPclockSetTime(targetscip->stat->resolveinstablelptime, SCIPclockGetTime(sourcescip->stat->resolveinstablelptime) + SCIPclockGetTime(targetscip->stat->resolveinstablelptime));
   targetscip->stat->nresolveinstablelps += sourcescip->stat->nresolveinstablelps;
   targetscip->stat->nresolveinstablelpiters += sourcescip->stat->nresolveinstablelpiters;

   SCIPclockSetTime(targetscip->stat->divinglptime, SCIPclockGetTime(sourcescip->stat->divinglptime) + SCIPclockGetTime(targetscip->stat->divinglptime));
   targetscip->stat->ndivinglps += sourcescip->stat->ndivinglps;
   targetscip->stat->ndivinglpiterations += sourcescip->stat->ndivinglpiterations;

   SCIPclockSetTime(targetscip->stat->strongbranchtime, SCIPclockGetTime(sourcescip->stat->strongbranchtime) + SCIPclockGetTime(targetscip->stat->strongbranchtime));
   targetscip->stat->nstrongbranchs += sourcescip->stat->nstrongbranchs;
   targetscip->stat->nsblpiterations += sourcescip->stat->nsblpiterations;

   targetscip->stat->nrootstrongbranchs += sourcescip->stat->nrootstrongbranchs;
   targetscip->stat->nrootsblpiterations += sourcescip->stat->nrootsblpiterations;

   SCIPclockSetTime(targetscip->stat->conflictlptime, SCIPclockGetTime(sourcescip->stat->conflictlptime) + SCIPclockGetTime(targetscip->stat->conflictlptime));
   targetscip->stat->nconflictlps += sourcescip->stat->nconflictlps;
   targetscip->stat->nconflictlpiterations += sourcescip->stat->nconflictlpiterations;
}


/** transfer NLP data */
static
SCIP_RETCODE SCIPtransferNLPStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   int nnlrowlinear;
   int nnlrowconvexineq;
   int nnlrownonconvexineq;
   int nnlrownonlineareq;

   SCIPclockSetTime(targetscip->stat->nlpsoltime, SCIPclockGetTime(sourcescip->stat->nlpsoltime) + SCIPclockGetTime(targetscip->stat->nlpsoltime));
   targetscip->stat->nnlps += sourcescip->stat->nnlps;

   if ( sourcescip->nlp != NULL )
   {
      SCIP_CALL( SCIPgetNLPNlRowsStat(sourcescip, &nnlrowlinear, &nnlrowconvexineq, &nnlrownonconvexineq, &nnlrownonlineareq) );

      if ( targetscip->nlp != NULL )
      {
         targetscip->nlp->nnlrowlinear += nnlrowlinear;
         targetscip->nlp->nnlrowconvexineq += nnlrowconvexineq;
         targetscip->nlp->nnlrownonconvexineq += nnlrownonconvexineq;
         targetscip->nlp->nnlrownonlineareq += nnlrownonlineareq;
      }
   }

   return SCIP_OKAY;
}


/** transfer relaxator data */
static
void SCIPtransferRelaxStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIP_RELAX** sourcerelaxs;
   int nsourcerelaxs;
   int i;

   sourcerelaxs = SCIPgetRelaxs(sourcescip);
   nsourcerelaxs = SCIPgetNRelaxs(sourcescip);

   for (i = 0; i < nsourcerelaxs; ++i)
   {
      SCIP_RELAX* sourcerelax;
      SCIP_RELAX* targetrelax;

assert( sourcerelaxs != NULL );
      sourcerelax = sourcerelaxs[i];
      assert( sourcerelax != NULL );

      targetrelax = SCIPfindRelax(targetscip, SCIPrelaxGetName(sourcerelax));

      if ( targetrelax != NULL )
      {
         SCIPclockSetTime(targetrelax->relaxclock, SCIPrelaxGetTime(sourcerelax) + SCIPrelaxGetTime(targetrelax));
         SCIPclockSetTime(targetrelax->setuptime, SCIPrelaxGetSetupTime(sourcerelax) + SCIPrelaxGetSetupTime(targetrelax));

         targetrelax->ncalls += SCIPrelaxGetNCalls(sourcerelax);
         targetrelax->ncutoffs += SCIPrelaxGetNCutoffs(sourcerelax);
         targetrelax->nimprbounds += SCIPrelaxGetNImprovedLowerbound(sourcerelax);
         targetrelax->imprtime += SCIPrelaxGetImprovedLowerboundTime(sourcerelax);
         targetrelax->nreduceddom += SCIPrelaxGetNReducedDomains(sourcerelax);
         targetrelax->nseparated += SCIPrelaxGetNSeparatedCuts(sourcerelax);
         targetrelax->naddedconss += SCIPrelaxGetNAddedConss(sourcerelax);
      }
   }
}


/** transfer tree data */
static
void SCIPtransferTreeStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   /* do not update nnodes. ninternalnodes because we have the total counts */

   targetscip->stat->nruns += sourcescip->stat->nruns;
   targetscip->stat->nfeasleaves += sourcescip->stat->nfeasleaves;
   targetscip->stat->ninfeasleaves += sourcescip->stat->ninfeasleaves;
   targetscip->stat->nobjleaves += sourcescip->stat->nobjleaves;
   targetscip->stat->ntotalnodes += sourcescip->stat->ntotalnodes;
   targetscip->stat->ntotalinternalnodes += sourcescip->stat->ntotalinternalnodes;
   targetscip->stat->maxdepth = MAX(targetscip->stat->maxdepth, sourcescip->stat->maxdepth);
   targetscip->stat->maxtotaldepth = MAX(targetscip->stat->maxtotaldepth, sourcescip->stat->maxtotaldepth);
   targetscip->stat->nbacktracks += sourcescip->stat->nbacktracks;
   targetscip->stat->nearlybacktracks += sourcescip->stat->nearlybacktracks;
   targetscip->stat->nnodesaboverefbound += sourcescip->stat->nnodesaboverefbound;
   targetscip->stat->ndelayedcutoffs += sourcescip->stat->ndelayedcutoffs;
   targetscip->stat->nreprops += sourcescip->stat->nreprops;
   targetscip->stat->nrepropboundchgs += sourcescip->stat->nrepropboundchgs;
   targetscip->stat->nrepropcutoffs += sourcescip->stat->nrepropcutoffs;
   targetscip->stat->nactivatednodes += sourcescip->stat->nactivatednodes;

   SCIPclockSetTime(targetscip->stat->nodeactivationtime, SCIPclockGetTime(targetscip->stat->nodeactivationtime) + SCIPclockGetTime(sourcescip->stat->nodeactivationtime));
}


/** transfers statistics from sourcescip to targetscip */
SCIP_RETCODE SCIPtransferStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIPtransferPresolStatistics(sourcescip, targetscip);
   SCIPtransferPropStatistics(sourcescip, targetscip);
   SCIPtransferConshdlrStatistics(sourcescip, targetscip);
   SCIPtransferConflictStatistics(sourcescip, targetscip);
   SCIPtransferSeparatorStatistics(sourcescip, targetscip);
   SCIPtransferCutSelectorStatistics(sourcescip, targetscip);
   SCIPtransferPricerStatistics(sourcescip, targetscip);
   SCIPtransferBranchruleStatistics(sourcescip, targetscip);
   SCIPtransferHeuristicsStatistics(sourcescip, targetscip);
   SCIPtransferCompressionStatistics(sourcescip, targetscip);
   SCIPtransferLPStatistics(sourcescip, targetscip);
   SCIP_CALL( SCIPtransferNLPStatistics(sourcescip, targetscip) );
   SCIPtransferRelaxStatistics(sourcescip, targetscip);
   SCIPtransferTreeStatistics(sourcescip, targetscip);

   /* do not update solution and root statistics as well as primal/dual integral, since we are not interested in the partial solutions */
   /* do not change timing statistics, since it already measures the total time */

   /* currently skip: benders, concurrent. expression, NLPI, reoptimization, branching, LNS statistics, diving heuristics */

   return SCIP_OKAY;
}
