// -*- C++ -*-
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

/**@file   presol_persistence.cpp
 * @ingroup PRESOLVERS
 * @brief  persistence presolver based on Nemhauser-Trotter
 * @author Erik Jansen
 * @author Annika Jaeger
 * @author Max Glaeser
 * @author Jonas Alker
 * @author Marc Pfetsch
 *
 * This presolver solves the LP arising from the edge relaxation. Nemhauser and Trotter proved that the vertices have
 * components with values in {0, 0.5, 1}. Moreover, there is a persistence property: if the components have value 0 or
 * 1, then there exits an integer optimal solution that has the same values. Thus, one can fix variables to these
 * values. Moreover, we can perform probing by tentatively fixing variables to 1 and repeating the process.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "presol_persistence.h"
#include "graph.h"
#include "struct_probdata.h"
#include "lpi/lpi.h"


#define PRESOL_NAME            "persistence"
#define PRESOL_DESC            "persistence presolver based on Nemhauser-Trotter"
#define PRESOL_PRIORITY         1000000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS              0 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_RUNPROBING       FALSE    //!< whether Nemhauser-Trotter presolving should be applied as probing


/*
 * Callback methods of presolver
 */

//! presolver data
struct SCIP_PresolData
{
   SCIP_Bool             runprobing;         //!< whether Nemhauser-Trotter presolving should be applied as probing
};


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecPersistence)
{
   assert( result != nullptr );
   *result = SCIP_DIDNOTRUN;

   // do not run, if not all variables are explicitly known
   if ( SCIPgetNActivePricers(scip) > 0 )
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

   SCIPdebugMsg(scip, "Running Nemhauser-Trotter presolving ...\n");
   *result = SCIP_DIDNOTFIND;

   // prepare edge formulation as LP
   SCIP_LPI* lpi;
   SCIP_MESSAGEHDLR* messagehdlr =SCIPgetMessagehdlr(scip);
   SCIP_CALL( SCIPlpiCreate(&lpi, messagehdlr, "edge_relaxation", SCIP_OBJSEN_MAXIMIZE) );

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );
   const Graph* G = probdata->G;
   assert( G != nullptr );

   // add variables - one for each node
   for (size_t i = 0; i < probdata->n; ++i)
   {
      const SCIP_Real obj = SCIPvarGetObj(probdata->vars[i]);
      const SCIP_Real lb = SCIPvarGetLbLocal(probdata->vars[i]);
      const SCIP_Real ub = SCIPvarGetUbLocal(probdata->vars[i]);

      SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, nullptr, 0, nullptr, nullptr, nullptr) );
   }

   // add edge constraints
   {
      const SCIP_Real lhs = -SCIPlpiInfinity(lpi);
      constexpr SCIP_Real rhs = 1.0;
      constexpr int beg = 0;
      constexpr SCIP_Real val[2] = {1.0, 1.0};
      int ind[2];

      EdgeIterator eit, eend;
      for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
      {
         Edge e = *eit;
         ind[0] = (int) boost::source(e, *G);
         ind[1] = (int) boost::target(e, *G);

         SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, nullptr, 2, &beg, ind, val) );
      }
   }

   // transfer time limit
   SCIP_Real timelimit;
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if ( ! SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      if ( timelimit > 0.0 )
      {
         SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_LPTILIM, timelimit) );
      }
   }

   // solve LP
   SCIPdebugMsg(scip, "Solving edge formulation ...\n");
   SCIP_CALL( SCIPlpiSolveDual(lpi) );
   SCIP_Real lpopt = SCIP_INVALID;
#ifdef SCIP_DEBUG
   int noldfixedvars = *nfixedvars;
#endif

   // variables needed for changing the problem
   constexpr SCIP_Real bound0 = 0.0;
   constexpr SCIP_Real bound1 = 1.0;

   // check solution status
   if ( SCIPlpiIsTimelimExc(lpi) )
   {
      SCIPwarningMessage(scip, "Time limit exceeded when solving the edge formulation in presolver <%s>.\n", PRESOL_NAME);
      SCIP_CALL( SCIPlpiFree(&lpi) );
      return SCIP_OKAY;
   }
   else if ( ! SCIPlpiIsOptimal(lpi) )
   {
      SCIPwarningMessage(scip, "Could not solve edge formulation to optimality in presolver <%s>.\n", PRESOL_NAME);
      SCIP_CALL( SCIPlpiFree(&lpi) );
      return SCIP_OKAY;
   }
   else
   {
      assert( SCIPlpiIsOptimal(lpi) );

      SCIP_CALL( SCIPlpiGetObjval(lpi, &lpopt) );

      // analyze solution by using results of Nemhauser-Trotter
      SCIP_Real* primsol;
      SCIP_CALL( SCIPallocBufferArray(scip, &primsol, probdata->n) );
      SCIP_CALL( SCIPlpiGetSol(lpi, nullptr, primsol, nullptr, nullptr, nullptr) );

      for (size_t i = 0; i < probdata->n; ++i)
      {
         // fix unfixed variables with value 0 in LP to 0
         if ( SCIPvarGetUbLocal(probdata->vars[i]) > 0.5 && SCIPisFeasEQ(scip, primsol[i], 0.0) )
         {
            const int ind = (int) i;
            SCIP_CALL( SCIPchgVarUb(scip, probdata->vars[i], 0.0) );
            SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &bound0, &bound0) );
            ++(*nfixedvars);
            *result = SCIP_SUCCESS;
            SCIPdebugMsg(scip, "fixed <%s> to 0.\n", SCIPvarGetName(probdata->vars[i]));
         }

         // fix unfixed variables with value 1 in LP to 1
         if ( SCIPvarGetLbLocal(probdata->vars[i]) < 0.5 && SCIPisFeasEQ(scip, primsol[i], 1.0) )
         {
            const int ind = (int) i;
            SCIP_CALL( SCIPchgVarLb(scip, probdata->vars[i], 1.0) );
            SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &bound1, &bound1) );
            ++(*nfixedvars);
            *result = SCIP_SUCCESS;
            SCIPdebugMsg(scip, "fixed <%s> to 1.\n", SCIPvarGetName(probdata->vars[i]));
         }
      }
      SCIPfreeBufferArray(scip, &primsol);
   }
   assert( lpopt != SCIP_INVALID );

   // perform probing: tentatively fix variables and see whether this does not change the objective
   SCIP_PRESOLDATA* presoldata = SCIPpresolGetData(presol);
   assert( presoldata != nullptr );
   if ( presoldata->runprobing )
   {
      for (size_t i = 0; i < probdata->n; ++i)
      {
         // skip already fixed variables
         if ( SCIPvarGetLbLocal(probdata->vars[i]) > 0.5 || SCIPvarGetUbLocal(probdata->vars[i]) < 0.5 )
            continue;

         SCIPdebugMsg(scip, "Tentatively fix <%s> to 1 ...\n", SCIPvarGetName(probdata->vars[i]));

         // fix variable to 1
         const int ind = (int) i;
         SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &bound1, &bound1) );

         // solve problem
         SCIP_CALL( SCIPlpiSolveDual(lpi) );

         // check solution status
         if ( SCIPlpiIsTimelimExc(lpi) )
         {
            SCIPwarningMessage(scip, "Time limit exceeded when solving the edge formulation in presolver <%s>.\n", PRESOL_NAME);
            break;
         }
         else if ( ! SCIPlpiIsOptimal(lpi) )
         {
            SCIPwarningMessage(scip, "Could not solve edge formulation to optimality in presolver <%s>.\n", PRESOL_NAME);
            break;
         }
         else
         {
            assert( SCIPlpiIsOptimal(lpi) );

            SCIP_Real lpfixedopt;
            SCIP_CALL( SCIPlpiGetObjval(lpi, &lpfixedopt) );

            // if the objective has not changed to original LP
            if ( SCIPisFeasEQ(scip, lpopt, lpfixedopt) )
            {
               // analyze solution
               SCIP_Real* primsol;
               SCIP_CALL( SCIPallocBufferArray(scip, &primsol, probdata->n) );

               SCIP_CALL( SCIPlpiGetSol(lpi, nullptr, primsol, nullptr, nullptr, nullptr) );
               for (size_t k = 0; k < probdata->n; ++k)
               {
                  // fix unfixed variables with value 0 in LP to 0
                  if ( SCIPvarGetUbLocal(probdata->vars[k]) > 0.5 && SCIPisFeasEQ(scip, primsol[k], 0.0) )
                  {
                     const int indd = (int) k;
                     SCIP_CALL( SCIPchgVarUb(scip, probdata->vars[k], 0.0) );
                     SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &indd, &bound0, &bound0) );
                     ++(*nfixedvars);
                     *result = SCIP_SUCCESS;
                     SCIPdebugMsg(scip, "fixed <%s> to 0.\n", SCIPvarGetName(probdata->vars[k]));
                  }

                  if ( SCIPvarGetLbLocal(probdata->vars[k]) < 0.5 && SCIPisFeasEQ(scip, primsol[k], 1.0) )
                  {
                     const int indd = (int) k;
                     SCIP_CALL( SCIPchgVarLb(scip, probdata->vars[k], 1.0) );
                     SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &indd, &bound1, &bound1) );
                     ++(*nfixedvars);
                     *result = SCIP_SUCCESS;
                     SCIPdebugMsg(scip, "fixed <%s> to 1.\n", SCIPvarGetName(probdata->vars[k]));
                  }
               }
               SCIPfreeBufferArray(scip, &primsol);
            }
            else
            {
               // reset bounds in lpi
               constexpr double lb_old = 0.0;
               constexpr double ub_old = 1.0;
               SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &lb_old, &ub_old) );
            }
         }
      }
   }
#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "Nemhauser-Trotter presolving fixed %d variables.\n", *nfixedvars - noldfixedvars);
#endif
   SCIP_CALL( SCIPlpiFree(&lpi) );

   return SCIP_OKAY;
}


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreePersistence)
{
   assert( scip != nullptr );
   assert( presol != nullptr );

   SCIP_PRESOLDATA* presoldata = SCIPpresolGetData(presol);
   assert( presoldata != nullptr );
   SCIPfreeBlockMemory(scip, &presoldata);

   return SCIP_OKAY;
}


/** creates the persistence presolver and includes it in SCIP */
SCIP_RETCODE BACSincludePresolPersistence(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata = nullptr;
   SCIP_PRESOL* presol = nullptr;

   // create comp data
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   /* create persistence presolver data */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecPersistence, presoldata) );
   assert( presol != nullptr );

   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreePersistence) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/runprobing",
         "whether Nemhauser-Trotter presolving should be applied as probing",
         &presoldata->runprobing, FALSE, DEFAULT_RUNPROBING, nullptr, nullptr) );

   return SCIP_OKAY;
}
