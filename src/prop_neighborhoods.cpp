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

/**@file   prop_neighborhoods.cpp
 * @ingroup PROPAGATORS
 * @brief  propagator to fix nodes to zero if neighborhood includes neighborhood of an neighbor
 * @author Erik Jansen
 * @author Annika Jaeger
 * @author Jonas Alker
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "prop_neighborhoods.h"
#include "graph.h"
#include "struct_probdata.h"

#include <scip/prop_symmetry.h>


#define PROP_NAME                  "neighborhoods"
#define PROP_DESC                  "propagator to fix nodes to zero if neighborhood includes neighborhood of a neighbor"
#define PROP_TIMING                 SCIP_PROPTIMING_BEFORELP
#define PROP_PRIORITY              -1000000  //!< priority of the propagator (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PROP_FREQ                         2  //!< propagation frequency
#define PROP_DELAY                     TRUE  //!< should propagation method be delayed, if other propagators found reductions?
#define PROP_PRESOL_PRIORITY       -9999999  //!< priority of the propagator (>= 0: before, < 0: after constraint handlers)
#define PROP_PRESOL_MAXROUNDS             1  //!< maximal number of propving rounds the propagator participates in (-1: no limit)
#define PROP_PRESOLTIMING           SCIP_PRESOLTIMING_MEDIUM //!< timing of the presolving method (fast, medium, or exhaustive)

#define DEFAULT_PROBINGONLY           FALSE  //!< whether the propagator should only run in probing
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
};


/*
 * Local methods
 */

//! perform propagation
static
SCIP_RETCODE propagateNeighborhoods(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Bool&            cutoff,             /**< whether we detected a cutoff */
   int&                  nfixed              /**< number of variables fixed */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );

   cutoff = FALSE;
   nfixed = 0;

   const Graph* G = probdata->G;
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex s = *vit;
      SCIP_VAR* svar = SCIPvarGetTransVar(probdata->vars[s]);

      // skip fixed node
      if ( SCIPvarGetUbLocal(svar) < 0.5 || SCIPvarGetLbLocal(svar) > 0.5 )
         continue;

      // get objective coefficient (use original variable, because aggregated variables have objective 0)
      SCIP_Real sobj = SCIPvarGetObj(probdata->vars[s]);

      // loop through neighbors
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(s, *G); ait != aend; ++ait)
      {
         Vertex t = *ait;
         SCIP_VAR* tvar = SCIPvarGetTransVar(probdata->vars[t]);

         // only check s < t, break as neighborhood of s is sorted
         if ( t > s )
            break;

         // skip fixed node
         if ( SCIPvarGetUbLocal(tvar) < 0.5 || SCIPvarGetLbLocal(tvar) > 0.5 )
            continue;

         // loop over neighborhood of s node
         AdjacencyIterator sit, send;
         boost::tie(sit, send) = boost::adjacent_vertices(s, *G);

         // loop over neighborhood of t node
         AdjacencyIterator tit, tend;
         boost::tie(tit, tend) = boost::adjacent_vertices(t, *G);

         // setting up Booleans deciding whether neighborhood of s is subset of t and vice versa
         SCIP_Bool sint = TRUE;
         SCIP_Bool tins = TRUE;

         // checking objective criteria (we have to take the original variables here, because aggregated variables have objective 0)
         SCIP_Real tobj = SCIPvarGetObj(probdata->vars[t]);

         // Booleans indicating whether a neighbor of s or t is already fixed to one
         SCIP_Bool sneighfixedone = FALSE;
         SCIP_Bool tneighfixedone = FALSE;

         if ( SCIPisGT(scip, sobj, tobj) )
            tins = FALSE;
         if ( SCIPisGT(scip, tobj, sobj) )
            sint = FALSE;

         // loop over both adjacency lists to determine inclusion in either direction
         while ( sit != send || tit != tend )
         {
            // if we already have non-inclusion in both directions
            if ( ! sint && ! tins )
               break;

            // s iterator at end
            if ( sit == send )
            {
               assert( tit != tend );
               while ( tit != tend && ( SCIPvarGetUbLocal(probdata->vars[*tit]) < 0.5 || *tit == s ) )
                  ++tit;

               if ( tit != tend )
               {
                  tins = FALSE;
                  if ( SCIPvarGetLbLocal(probdata->vars[*tit]) > 0.5 )
                  {
                     tneighfixedone = TRUE;
                     break;
                  }
               }
               break;
            }

            // t iterator at end
            if ( tit == tend )
            {
               assert( sit != send );
               while ( sit != send && ( SCIPvarGetUbLocal(probdata->vars[*sit]) < 0.5 || *sit == t ) )
                  ++sit;

               if ( sit != send )
               {
                  sint = FALSE;
                  if ( SCIPvarGetLbLocal(probdata->vars[*sit]) > 0.5 )
                  {
                     sneighfixedone = TRUE;
                     break;
                  }
               }
               break;
            }

            // neighbor of s already fixed to one
            if ( SCIPvarGetLbLocal(probdata->vars[*sit]) > 0.5 )
            {
               sneighfixedone = TRUE;
               break;
            }

            // neighbor of t already fixed to one
            if ( SCIPvarGetLbLocal(probdata->vars[*tit]) > 0.5 )
            {
               tneighfixedone = TRUE;
               break;
            }

            // advance iterators over shared neighborhood
            if ( *sit == *tit )
            {
               ++sit;
               ++tit;
               continue;
            }

            // s has neighbor which t does not have
            if ( *sit < *tit )
            {
               if ( SCIPvarGetUbLocal(probdata->vars[*sit]) > 0.5 && *sit != t )
                  sint = FALSE;
               ++sit;
               continue;
            }

            // t has neighbor which s does not have
            if ( *sit > *tit )
            {
               if ( SCIPvarGetUbLocal(probdata->vars[*tit]) > 0.5 && *tit != s )
                  tins = FALSE;
               ++tit;
               continue;
            }

            // previous cases are exhaustive
            abort();
         }

         // fix s to zero if neighbor fixed to one has been found
         if ( sneighfixedone )
         {
            SCIPdebugMsg(scip, "Fixing node <%s> to 0.\n", SCIPvarGetName(svar));

            // We currently directly fix the variable without conflict analysis, since this case will rarely happen and
            // it would make reverse propagation more complicated.
            SCIP_Bool tightened = FALSE;
            SCIP_CALL( SCIPtightenVarUb(scip, svar, 0.0, FALSE, &cutoff, &tightened) );
            assert( ! cutoff );
            assert( tightened );
            ++nfixed;
            break; // exit loop for t, since s is fixed now
         }

         // fix t to zero if neighbor fixed to one has been found
         if ( tneighfixedone )
         {
            SCIPdebugMsg(scip, "Fixing node <%s> to 0.\n", SCIPvarGetName(tvar));

            // We currently directly fix the variable without conflict analysis, since this case will rarely happen and
            // it would make reverse propagation more complicated.
            SCIP_Bool tightened = FALSE;
            SCIP_CALL( SCIPtightenVarUb(scip, tvar, 0.0, FALSE, &cutoff, &tightened) );
            assert( ! cutoff );
            assert( tightened );
            ++nfixed;
            continue;
         }

         // fix one variable if condition is satisfied
         if ( sint )
         {
            SCIPdebugMsg(scip, " -> can fix <%s> to 0, because its neighborhood contains the one of <%s>.\n",
               SCIPvarGetName(tvar), SCIPvarGetName(svar));

            // neighborhood of s is subset of neighborhood of t and t is not better for obj -> fix t to 0.0
            SCIP_Bool tightened = FALSE;
            SCIP_CALL( SCIPinferBinvarProp(scip, tvar, FALSE, prop, (int) s, &cutoff, &tightened) );
            ++nfixed;
            if ( cutoff )
            {
               SCIPdebugMsg(scip, " -> detected infeasibility with nodes <%s> and <%s>.\n",
                  SCIPvarGetName(svar), SCIPvarGetName(tvar));
               return SCIP_OKAY;
            }
            assert( tightened );
         }
         else if ( tins )
         {
            SCIPdebugMsg(scip, " -> can fix <%s> to 0, because its neighborhood contains the one of <%s>.\n",
               SCIPvarGetName(svar), SCIPvarGetName(tvar));

            // neighborhood of t is subset of neighborhood of s and s is not better for obj -> fix s to 0.0
            SCIP_Bool tightened = FALSE;
            SCIP_CALL( SCIPinferBinvarProp(scip, svar, FALSE, prop, (int) t, &cutoff, &tightened) );
            ++nfixed;
            if ( cutoff )
            {
               SCIPdebugMsg(scip, " -> detected infeasibility with nodes <%s> and <%s>.\n",
                  SCIPvarGetName(svar), SCIPvarGetName(tvar));
               return SCIP_OKAY;
            }
            assert( tightened );
            break; // exit loop for t since s is fixed now
         }
      }
   }

   return SCIP_OKAY;
}



/*
 * Callback methods of propagator
 */

//! presolving initialization method of propagator (called when presolving is about to begin)
static
SCIP_DECL_PROPINITPRE(propInitpreNeighborhoods)
{
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
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "Turned off neighborhoods presolver, because the problem does not contain exaclty one clique constraint.\n");

   return SCIP_OKAY;
}

//! solving process initialization method of propagator (called when branch and bound process is about to begin)
static
SCIP_DECL_PROPINITSOL(propInitsolNeighborhoods)
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
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "Turning off symmetry handling to avoid conflicts with neighborhood propagation.\n");
         else
            assert( retcode == SCIP_PARAMETERUNKNOWN );
      }
   }

   return SCIP_OKAY;
}


//! destructor of propagator to free user data (called when SCIP is exiting)
static
SCIP_DECL_PROPFREE(propFreeNeighborhoods)
{
   assert( prop != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   SCIPfreeBlockMemory(scip, &propdata);

   return SCIP_OKAY;
}


//! presolving method
static
SCIP_DECL_PROPPRESOL(propPresolNeighborhoods)
{
   assert( scip != nullptr );
   assert( result != nullptr );
   assert( nfixedvars != nullptr );

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

#ifndef NDEBUG
   // assert AdjacencyIterator gives nodes in ascending order as this is crucial for following neighborhood search
   const Graph* G = probdata->G;
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      AdjacencyIterator nit, nend;
      boost::tie(nit, nend) = boost::adjacent_vertices(*vit, *G);
      Vertex nprev;
      while ( nit != nend )
      {
         nprev = *nit;
         ++nit;
         assert( nit == nend || nprev < *nit );
      }
   }
#endif

   // do not run if disabled
   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );
   if ( ! propdata->enabled )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "running neighborhood presolving.\n");

   SCIP_Bool cutoff;
   int nfixed;
   SCIP_CALL( propagateNeighborhoods(scip, prop, probdata, cutoff, nfixed) );

   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( nfixed > 0 )
   {
      *result = SCIP_SUCCESS;
      *nfixedvars += nfixed;
   }
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


//! execution method of propagator
static
SCIP_DECL_PROPEXEC(propExecNeighborhoods)
{  /*lint --e{715}*/
   assert( scip != nullptr );
   assert( result != nullptr );
   assert( prop != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   // if propagation should only run in probing
   if ( propdata->probingonly )
   {
      if ( ! SCIPinProbing(scip) || SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING )
      {
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }
   }

   // do not run if disabled
   if ( ! propdata->enabled )
      return SCIP_OKAY;

   // possibly only run in the tree if there are no symmetries
   if ( propdata->onlyfornosymmetry && SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING && SCIPgetSymmetryNGenerators(scip) > 0 )
      return SCIP_OKAY;

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   if ( SCIPinProbing(scip) )
      SCIPdebugMsg(scip, "running neighborhood propagation in probing.\n");
   else
      SCIPdebugMsg(scip, "running neighborhood propagation.\n");

   SCIP_Bool cutoff;
   int nfixed;
   SCIP_CALL( propagateNeighborhoods(scip, prop, probdata, cutoff, nfixed) );

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
SCIP_DECL_PROPRESPROP(propRespropNeighborhoods)
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
   SCIPdebugMsg(scip, "neighborhood of <%s> was contained in <%s>.\n", SCIPvarGetName(svar), SCIPvarGetName(tvar));

   // check preconditions for original propagation
   assert( SCIPvarGetTransVar(probdata->vars[t]) == tvar );
   assert( boost::edge(s, t, *G).second );
   // The following two asserts might fail, because conflict analysis can determine that the variables could have been fixed
   // to 0 higher up in the tree.
   // assert( SCIPgetVarUbAtIndex(scip, svar, bdchgidx, FALSE) > 0.5 );
   // assert( SCIPgetVarUbAtIndex(scip, tvar, bdchgidx, FALSE) > 0.5 );
   assert( SCIPgetVarLbAtIndex(scip, svar, bdchgidx, FALSE) < 0.5 );
   assert( SCIPgetVarLbAtIndex(scip, tvar, bdchgidx, FALSE) < 0.5 );
   assert( SCIPisLE(scip, SCIPvarGetObj(probdata->vars[t]), SCIPvarGetObj(probdata->vars[s])) );

   // set up iterator for neighborhood of t
   AdjacencyIterator tit, tend;
   boost::tie(tit, tend) = boost::adjacent_vertices(t, *G);

   // loop over neighborhood of s
   AdjacencyIterator sit, send;
   for (boost::tie(sit, send) = boost::adjacent_vertices(s, *G); sit != send; ++sit)
   {
      if ( tit != tend )
      {
         // ignore common neighbors
         if ( *sit == *tit )
            continue;

         if ( *sit == t )
            continue;

         // s has neighbor which t does not have - thus this variable needs to be fixed to 0 for a valid reduction
         if ( *sit < *tit )
         {
            SCIP_VAR* var = probdata->vars[*sit];
            assert( SCIPvarGetUbLocal(var) < 0.5 && *sit != t );
            SCIPdebugMsg(scip, "fixings in neighborhood of <%s>: <%s> = 0.\n", SCIPvarGetName(svar), SCIPvarGetName(var));
            SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
            *result = SCIP_SUCCESS;
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

//! creates the neighborhoods propagator and includes it in SCIP
SCIP_RETCODE BACSincludePropNeighborhoods(
   SCIP*                 scip                //! SCIP data structure
   )
{
   SCIP_PROPDATA* propdata = nullptr;
   SCIP_PROP* prop = nullptr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );
   propdata->enabled = TRUE;

   // include propagator
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecNeighborhoods, propdata) );
   assert( prop != nullptr );

   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolNeighborhoods, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS,
         PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropInitpre(scip, prop, propInitpreNeighborhoods) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropNeighborhoods) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeNeighborhoods) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolNeighborhoods) );

   // add parameters
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/probingonly",
         "whether the propagator should only run in probing",
         &propdata->probingonly, FALSE, DEFAULT_PROBINGONLY, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/onlyfornosymmetry",
         "whether the propagator should only be run in the tree if there are no symmetries",
         &propdata->onlyfornosymmetry, FALSE, DEFAULT_ONLYFORNOSYMMETRY, nullptr, nullptr) );

   return SCIP_OKAY;
}
