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

/**@file   prop_lowdegreenodes.cpp
 * @ingroup PROPAGATORS
 * @brief  propagator to fix nodes with low local degree
 * @author Jonas Alker
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>

#include "prop_lowdegreenodes.h"
#include "graph.h"
#include "struct_probdata.h"

#include <scip/prop_symmetry.h>

#define PROP_NAME                  "lowdegreenodes"
#define PROP_DESC                  "propagator to fix nodes with low local degree"
#define PROP_TIMING                 SCIP_PROPTIMING_BEFORELP
#define PROP_PRIORITY               -300000  //!< priority of the propagator (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PROP_FREQ                         2  //!< propagation frequency
#define PROP_DELAY                    FALSE  //!< should propagation method be delayed, if other propagators found reductions?
#define PROP_PRESOL_PRIORITY        -300000  //!< priority of the propagator (>= 0: before, < 0: after constraint handlers)
#define PROP_PRESOL_MAXROUNDS            -1  //!< maximal number of presolving rounds the propagator participates in (-1: no limit)
#define PROP_PRESOLTIMING           SCIP_PRESOLTIMING_FAST //!< timing of the presolving method (fast, medium, or exhaustive)


#define DEFAULT_ONLYFORNOSYMMETRY      TRUE  //!< whether the propagator should only be run in the tree if there are no symmetries
#define DEFAULT_CONFLICTANALYSIS      FALSE  //!< whether conflict analysis should be applied

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_Bool             enabled;            //!< whether the propagator is enabled
   SCIP_Bool             onlyfornosymmetry;  //!< whether the propagator should only be run in the tree if there are no symmetries
   SCIP_Bool             conflictanalysis;   //!< whether conflict analysis should be applied
   size_t                nfixeddegzero;      //!< number of isolated nodes fixed
   size_t                nfixeddegone;       //!< number of nodes fixed because of degree one node
   size_t                nfixeddegtwo;       //!< number of nodes fixed because of degree two node
   size_t                nfixeddegthree;     //!< number of nodes fixed because of degree three node
   size_t                nfixedneighbor;     //!< number of nodes fixed because neighbor is already fixed to one
};

/*
 * Local methods
 */

//! fixing variable to zero and incrementing nfixed
static
SCIP_RETCODE fixVarToZero(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   Vertex                v,                  /**< node to fix */
   int                   inferinfo,          /**< inferinfo for conflict analysis */
   size_t&               nfixed              /**< fixed nodes counter to increment */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( v < probdata->n );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   SCIP_VAR* var = probdata->vars[v];
   assert( var != nullptr );
   assert( SCIPvarGetLbLocal(var) < 0.5 );

   if ( SCIPvarGetUbLocal(var) < 0.5 )
      return SCIP_OKAY;

   SCIP_Bool cutoff;
   SCIP_Bool tightened;

   if ( propdata->conflictanalysis )
      SCIP_CALL( SCIPinferBinvarProp(scip, var, FALSE, prop, inferinfo, &cutoff, &tightened) );
   else
      SCIP_CALL( SCIPfixVar(scip, var, 0.0, &cutoff, &tightened ) );

   SCIPdebugMsg(scip, "fixing node <%s> to 0.0.\n", SCIPvarGetName(var));
   assert( ! cutoff );
   assert( tightened );
   ++nfixed;

   return SCIP_OKAY;
}

//! fixing variable to one and incrementing nfixed
static
SCIP_RETCODE fixVarToOne(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   Vertex                v,                  /**< node to fix */
   int                   inferinfo,          /**< inferinfo for conflict analysis */
   size_t&               nfixed              /**< fixed nodes counter to increment */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( v < probdata->n );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   SCIP_VAR* var = probdata->vars[v];
   assert( var != nullptr );
   assert( SCIPvarGetUbLocal(var) > 0.5 );

   if ( SCIPvarGetLbLocal(var) > 0.5 )
      return SCIP_OKAY;

   SCIP_Bool cutoff;
   SCIP_Bool tightened;

   if ( propdata->conflictanalysis )
      SCIP_CALL( SCIPinferBinvarProp(scip, var, TRUE, prop, inferinfo, &cutoff, &tightened) );
   else
      SCIP_CALL( SCIPfixVar(scip, var, 1.0, &cutoff, &tightened ) );

   SCIPdebugMsg(scip, "fixing node <%s> to 1.0.\n", SCIPvarGetName(var));
   assert( ! cutoff );
   assert( tightened );
   ++nfixed;

   return SCIP_OKAY;
}

//! perform propagation
static
SCIP_RETCODE propagateLowdegreenodes(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   size_t&               nfixedzero,         /**< number of isolated nodes fixed */
   size_t&               nfixedone,          /**< number of nodes fixed because of degree one node */
   size_t&               nfixedtwo,          /**< number of nodes fixed because of degree two node */
   size_t&               nfixedthree,        /**< number of nodes fixed because of degree three node */
   size_t&               nfixedneigh         /**< number of nodes fixed because neighbor is already fixed to one */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );

   nfixedzero = 0;
   nfixedone = 0;
   nfixedtwo = 0;
   nfixedthree = 0;
   nfixedneigh = 0;

   const Graph* G = probdata->G;

   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      SCIP_VAR* var = probdata->vars[v];
      assert( probdata->localdegrees[v] < 0 || SCIPisGT(scip, SCIPvarGetObj(var), 0.0) );

      if ( probdata->localdegrees[v] < 0 )
      {
         assert( SCIPvarGetUbLocal(var) < 0.5 || SCIPvarGetLbLocal(var) > 0.5 );
      }
      else if ( probdata->localdegrees[v] == 0 )
      {
         SCIP_Bool fix = TRUE;
         AdjacencyIterator ait, aend;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            assert( probdata->localdegrees[*ait] < 0 );

            // check for neighbor fixed to 1
            if ( SCIPvarGetLbLocal(probdata->vars[*ait]) > 0.5 )
            {
               SCIP_CALL( fixVarToZero(scip, prop, probdata, v, 0, nfixedneigh) );
               fix = FALSE;

               break;
            }
         }

         if ( fix )
            SCIP_CALL( fixVarToOne(scip, prop, probdata, v, 0, nfixedzero) );
      }
      else if ( probdata->localdegrees[v] == 1 )
      {
         Vertex w = Graph::null_vertex();
         SCIP_Bool fix = TRUE;

         AdjacencyIterator ait, aend;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            // check for neighbor fixed to 1
            if ( SCIPvarGetLbLocal(probdata->vars[*ait]) > 0.5 )
            {
               SCIP_CALL( fixVarToZero(scip, prop, probdata, v, 0, nfixedneigh) );
               fix = FALSE;
               break;
            }

            if ( SCIPvarGetUbLocal(probdata->vars[*ait]) > 0.5 )
            {
               assert( w == Graph::null_vertex() );
               w = *ait;
            }
         }
         if ( ! fix )
            continue;

         assert( w != Graph::null_vertex() );

         // w is better than v
         if ( SCIPisLT(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[w])) )
            continue;

         SCIP_CALL( fixVarToOne(scip, prop, probdata, v, 1, nfixedone) );
         SCIP_CALL( fixVarToZero(scip, prop, probdata, w, 1, nfixedone) );
      }
      else if ( probdata->localdegrees[v] == 2 )
      {
         Vertex u = Graph::null_vertex();
         Vertex w = Graph::null_vertex();
         SCIP_Bool fix = TRUE;

         AdjacencyIterator ait, aend;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            // check for neighbor fixed to 1
            if ( SCIPvarGetLbLocal(probdata->vars[*ait]) > 0.5 )
            {
               SCIP_CALL( fixVarToZero(scip, prop, probdata, v, 0, nfixedneigh) );
               fix = FALSE;
               break;
            }

            if ( SCIPvarGetUbLocal(probdata->vars[*ait]) > 0.5 )
            {
               // third unfixed var has been found
               assert( u == Graph::null_vertex() || w == Graph::null_vertex() );

               if ( u != Graph::null_vertex() )
                  w = *ait;
               else
                  u = *ait;
            }
         }
         if ( ! fix )
            continue;

         assert( u != Graph::null_vertex() && w != Graph::null_vertex() );
         assert( u != w );

         // u or w is better than v:
         if ( SCIPisLT(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[u])) || SCIPisLT(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[w])) )
            continue;

         // check if u, v, w are a triangle:
         SCIP_Bool triangle = FALSE;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(u, *G); ait != aend; ++ait)
         {
            if ( *ait == w )
            {
               triangle = TRUE;
               break;
            }
         }

         if ( ! triangle && SCIPisLT(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[u]) + SCIPvarGetObj(probdata->vars[w])) )
            continue;

         SCIP_CALL( fixVarToOne(scip, prop, probdata, v, 2, nfixedtwo) );
         SCIP_CALL( fixVarToZero(scip, prop, probdata, u, 2, nfixedtwo) );
         SCIP_CALL( fixVarToZero(scip, prop, probdata, w, 2, nfixedtwo) );
      }
      else if ( probdata->localdegrees[v] == 3 )
      {
         Vertex x = Graph::null_vertex();
         Vertex y = Graph::null_vertex();
         Vertex z = Graph::null_vertex();
         SCIP_Bool fix = TRUE;

         AdjacencyIterator ait, aend;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            // neighbor fixed to one
            if ( SCIPvarGetLbLocal(probdata->vars[*ait]) > 0.5 )
            {
               SCIP_CALL( fixVarToZero(scip, prop, probdata, v, 0, nfixedneigh) );
               fix = FALSE;
               break;
            }

            if ( SCIPvarGetUbLocal(probdata->vars[*ait]) > 0.5 )
            {
               // forth unfixed var has been found
               assert( x == Graph::null_vertex() || y == Graph::null_vertex() || z == Graph::null_vertex() );

               if ( x != Graph::null_vertex() && y != Graph::null_vertex() )
                  z = *ait;
               else if ( x != Graph::null_vertex() )
                  y = *ait;
               else
                  x = *ait;
            }
         }
         if ( ! fix )
            continue;

         assert( x != Graph::null_vertex() && y != Graph::null_vertex() && z != Graph::null_vertex() );
         assert( x != y && x != z && y != z );

         // x, y or z is better than v:
         if ( SCIPisLT(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[x])) \
            || SCIPisLT(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[y])) \
            || SCIPisLT(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[z])) )
            continue;

         // check connectivity of x, y and z:
         SCIP_Bool xandy = FALSE;
         SCIP_Bool xandz = FALSE;
         SCIP_Bool yandz = FALSE;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(x, *G); ait != aend; ++ait)
         {
            if ( *ait == y )
               xandy = TRUE;
            else if ( *ait == z )
               xandz = TRUE;

            if ( xandy && xandz )
               break;
         }

         for (boost::tie(ait, aend) = boost::adjacent_vertices(y, *G); ait != aend; ++ait)
         {
            if ( *ait == z )
            {
               yandz = TRUE;
               break;
            }
         }

         SCIP_Bool fixv = FALSE; // should v be fixed to one
         SCIP_Bool fixx = FALSE; // should x be fixed to zero
         SCIP_Bool fixy = FALSE; // should y be fixed to zero
         SCIP_Bool fixz = FALSE; // should z be fixed to zero
         if ( xandy && xandz && yandz )
         {  // v is simplicial, everything can be fixed
               fixv = TRUE;
               fixx = TRUE;
               fixy = TRUE;
               fixz = TRUE;
         }
         else if ( xandy && xandz )
         {  // fix node that is connected to both other neighbors of v to zero
            // then check whether v is better than both remaining neighbors
            assert( ! yandz );

            fixx = TRUE;
            if ( SCIPisGE(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[y]) + SCIPvarGetObj(probdata->vars[z])) )
            {
               fixv = TRUE;
               fixy = TRUE;
               fixz = TRUE;
            }
         }
         else if ( xandy && yandz )
         {
            assert( ! xandz );

            fixy = TRUE;
            if ( SCIPisGE(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[x]) + SCIPvarGetObj(probdata->vars[z])) )
            {
               fixv = TRUE;
               fixx = TRUE;
               fixz = TRUE;
            }
         }
         else if ( xandz && yandz )
         {
            assert( ! xandy );

            fixz = TRUE;
            if ( SCIPisGE(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[x]) + SCIPvarGetObj(probdata->vars[y])) )
            {
               fixv = TRUE;
               fixx = TRUE;
               fixy = TRUE;
            }
         }
         else if ( xandy || xandz || yandz )
         {  // everything can be fixed if v is better than the best stable set in its neighborhood

            assert( ! ( xandy && xandz ) && ! ( xandy && yandz ) && ! ( xandz && yandz ) );

            if ( (xandy && SCIPisGE(scip, SCIPvarGetObj(var), MAX( SCIPvarGetObj(probdata->vars[x]), SCIPvarGetObj(probdata->vars[y]) ) + SCIPvarGetObj(probdata->vars[z]))) \
               || (xandz && SCIPisGE(scip, SCIPvarGetObj(var), MAX( SCIPvarGetObj(probdata->vars[x]), SCIPvarGetObj(probdata->vars[z]) ) + SCIPvarGetObj(probdata->vars[y]))) \
               || (yandz && SCIPisGE(scip, SCIPvarGetObj(var), MAX( SCIPvarGetObj(probdata->vars[y]), SCIPvarGetObj(probdata->vars[z]) ) + SCIPvarGetObj(probdata->vars[x]))) )
            {
               fixv = TRUE;
               fixx = TRUE;
               fixy = TRUE;
               fixz = TRUE;
            }
         }
         else
         {  // everything can be fixed if v is better than its neighborhood
            assert( ! xandy );
            assert( ! xandz );
            assert( ! yandz );
            if ( SCIPisGE(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[x]) + SCIPvarGetObj(probdata->vars[y]) + SCIPvarGetObj(probdata->vars[z])) )
            {
               fixv = TRUE;
               fixx = TRUE;
               fixy = TRUE;
               fixz = TRUE;
            }
         }

         // apply fixings
         if ( fixv )
            SCIP_CALL( fixVarToOne(scip, prop, probdata, v, 3, nfixedthree) );

         if ( fixx )
            SCIP_CALL( fixVarToZero(scip, prop, probdata, x, 3, nfixedthree) );

         if ( fixy )
            SCIP_CALL( fixVarToZero(scip, prop, probdata, y, 3, nfixedthree) );

         if ( fixz )
            SCIP_CALL( fixVarToZero(scip, prop, probdata, z, 3, nfixedthree) );
      }
   }

   return SCIP_OKAY;
}



/*
 * Callback methods of propagator
 */

//! presolving initialization method of propagator (called when presolving is about to begin)
static
SCIP_DECL_PROPINITPRE(propInitpreLowdegreenodes)
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
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "Turned off lowdegreenodes presolver, because the problem does not contain exaclty one clique constraint.\n");

   return SCIP_OKAY;
}

//! solving process initialization method of propagator (called when branch and bound process is about to begin)
static
SCIP_DECL_PROPINITSOL(propInitsolLowdegreenodes)
{
   assert( prop != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   // if propagation should be performed in the tree, we currently turn off symmetry handling to avoid conflicts, unless we turn it off below anyways
   if ( SCIPpropGetFreq(prop) >= 0 && ! propdata->onlyfornosymmetry )
   {
      SCIP_RETCODE retcode = SCIPsetIntParam(scip, "propagating/symmetry/freq", -1);
      if ( retcode == SCIP_OKAY )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "Turning off symmetry handling to avoid conflicts with lowdegreenodes propagation.\n");
      else
         assert( retcode == SCIP_PARAMETERUNKNOWN );
   }

   return SCIP_OKAY;
}

//! destructor of propagator to free user data (called when SCIP is exiting)
static
SCIP_DECL_PROPFREE(propFreeLowdegreenodes)
{
   assert( prop != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "\nFixed %lu, %lu, %lu, %lu variables because of nodes with local degree 0, 1, 2, 3.\n", \
      propdata->nfixeddegzero, propdata->nfixeddegone, propdata->nfixeddegtwo, propdata->nfixeddegthree);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "Fixed %lu vars to zero because of fixed neighbor.\n", propdata->nfixedneighbor);

   SCIPfreeBlockMemory(scip, &propdata);

   return SCIP_OKAY;
}


//! presolving method
static
SCIP_DECL_PROPPRESOL(propPresolLowdegreenodes)
{
   assert( scip != nullptr );
   assert( result != nullptr );
   assert( nfixedvars != nullptr );

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   // do not run if disabled
   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );
   if ( ! propdata->enabled )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "running lowdegreenodes presolving.\n");

   size_t nfixedzero;
   size_t nfixedone;
   size_t nfixedtwo;
   size_t nfixedthree;
   size_t nfixedneigh;
   SCIP_CALL( propagateLowdegreenodes(scip, prop, probdata, nfixedzero, nfixedone, nfixedtwo, nfixedthree, nfixedneigh) );

   if ( nfixedzero + nfixedone + nfixedtwo + nfixedthree + nfixedneigh > 0 )
   {
      *result = SCIP_SUCCESS;
      if ( ! SCIPinProbing(scip) )
      {
         propdata->nfixeddegzero += nfixedzero;
         propdata->nfixeddegone += nfixedone;
         propdata->nfixeddegtwo += nfixedtwo;
         propdata->nfixeddegthree += nfixedthree;
         propdata->nfixedneighbor += nfixedneigh;
      }
      *nfixedvars += (int) (nfixedzero + nfixedone + nfixedtwo + nfixedthree + nfixedneigh);
      if ( nfixedzero + nfixedone + nfixedtwo + nfixedthree > 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "   (%.1fs) Fixed %lu, %lu, %lu, %lu variables because of nodes with degree 0, 1, 2, 3.\n", \
            SCIPgetSolvingTime(scip), nfixedzero, nfixedone, nfixedtwo, nfixedthree);
      }
      if ( nfixedneigh > 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "   (%.1fs) Fixed %lu vars to zero because of neighbor fixed to 1.\n", \
            SCIPgetSolvingTime(scip), nfixedneigh);
      }
   }
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


//! execution method of propagator
static
SCIP_DECL_PROPEXEC(propExecLowdegreenodes)
{  /*lint --e{715}*/
   assert( scip != nullptr );
   assert( result != nullptr );
   assert( prop != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   *result = SCIP_DIDNOTRUN;

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
      SCIPdebugMsg(scip, "running lowdegreenodes propagation in probing.\n");
   else
      SCIPdebugMsg(scip, "running lowdegreenodes propagation.\n");

   size_t nfixedzero;
   size_t nfixedone;
   size_t nfixedtwo;
   size_t nfixedthree;
   size_t nfixedneigh;
   SCIP_CALL( propagateLowdegreenodes(scip, prop, probdata, nfixedzero, nfixedone, nfixedtwo, nfixedthree, nfixedneigh) );

   if ( nfixedzero + nfixedone + nfixedtwo + nfixedthree + nfixedneigh > 0 )
   {
      *result = SCIP_REDUCEDDOM;
      if ( ! SCIPinProbing(scip) )
      {
         propdata->nfixeddegzero += nfixedzero;
         propdata->nfixeddegone += nfixedone;
         propdata->nfixeddegtwo += nfixedtwo;
         propdata->nfixeddegthree += nfixedthree;
         propdata->nfixedneighbor += nfixedneigh;
      }
   }
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropLowdegreenodes)
{
   assert( scip != nullptr );
   assert( prop != nullptr );
   assert( infervar != nullptr );
   assert( bdchgidx != nullptr );
   assert( result != nullptr );

   *result = SCIP_DIDNOTFIND;

   assert ( SCIPisConflictAnalysisApplicable(scip) );

   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

#ifndef NDEBUG
   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   assert( propdata->conflictanalysis );
#endif

   // inferinfos refer to localdegree of fixed variable when it was fixed
   assert( inferinfo >= 0 );
   assert( inferinfo <= 3 );

   const Graph* G = probdata->G;
   Vertex v = (Vertex) SCIPvarGetData(infervar);
   SCIP_VAR* var = infervar;

   assert( v < probdata->n );

   if ( boundtype == SCIP_BOUNDTYPE_UPPER )
   {  // find neighbor fixed to one
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;
         if ( SCIPgetVarLbAtIndex(scip, probdata->vars[w], bdchgidx, TRUE) > 0.5 )
         {
            // if v was fixed to 0 not because w was fixed to 1 already, w has to be at least as good as v
            if ( inferinfo > 0 && SCIPisGT(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[w])) )
               continue;

            SCIP_CALL( SCIPaddConflictBinvar(scip, probdata->vars[w]) );
            *result = SCIP_SUCCESS;
            break;
         }
      }

      if ( *result == SCIP_DIDNOTFIND )
      {  // in case 3 we might fix a variable to zero, even though no neighbor will be set to one
         assert( inferinfo == 3 );

         SCIPdebugMsg(scip, "Node at depth: %d, lowdegreenodes reverse propagation method at var %lu: not implemented case occurs\n", SCIPgetDepth(scip), v);

         // this case rarely happens:
         // v has three neighbors {u, x, y} \subset N(v) such that
         // u has exactly three neighbors {v, x, y} = N(u)
         // u is at least as good as all of his neighbors: c(u) \geq c(v), c(x), c(y)
         // x and y are not connected and c(x) + c(y) > c(u)

         // therefore we would never choose v over u but we can not fix u to 1 either due to x and y
      }
   }
   else if ( boundtype == SCIP_BOUNDTYPE_LOWER )
   {  // add all neighbors fixed to zero to conflict cons, assert that inferinfo vertices remain
      int unfixed = 0;
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;
         assert( SCIPgetVarLbAtIndex(scip, probdata->vars[w], bdchgidx, FALSE) < 0.5 );

         if ( SCIPgetVarUbAtIndex(scip, probdata->vars[w], bdchgidx, FALSE) < 0.5 )
            SCIP_CALL( SCIPaddConflictBinvar(scip, probdata->vars[w]) );
         else
         {
            ++unfixed;
            assert( SCIPisGE(scip, SCIPvarGetObj(var), SCIPvarGetObj(probdata->vars[w])) );
         }
      }
      // it may happen that other fixings are valid earlier leading here to fewer unfixed neighbors
      assert( unfixed <= inferinfo );

      *result = SCIP_SUCCESS;
   }
   assert( *result == SCIP_SUCCESS || ( inferinfo == 3 && boundtype == SCIP_BOUNDTYPE_UPPER ) );

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

//! creates the lowdegreenodess propagator and includes it in SCIP
SCIP_RETCODE BACSincludePropLowdegreenodes(
   SCIP*                 scip                //! SCIP data structure
   )
{
   SCIP_PROPDATA* propdata = nullptr;
   SCIP_PROP* prop = nullptr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );
   propdata->enabled = TRUE;
   propdata->nfixeddegzero = 0;
   propdata->nfixeddegone = 0;
   propdata->nfixeddegtwo = 0;
   propdata->nfixeddegthree = 0;
   propdata->nfixedneighbor = 0;

   // include propagator
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecLowdegreenodes, propdata) );
   assert( prop != nullptr );

   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolLowdegreenodes, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS,
         PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropInitpre(scip, prop, propInitpreLowdegreenodes) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeLowdegreenodes) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolLowdegreenodes) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropLowdegreenodes) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "propagating/" PROP_NAME "/onlyfornosymmetry",
      "whether the propagator should only be run in the tree if there are no symmetries",
      &propdata->onlyfornosymmetry, FALSE, DEFAULT_ONLYFORNOSYMMETRY, nullptr, nullptr) );

      SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/conflictanalysis",
         "whether conflict analysis should be applied",
         &propdata->conflictanalysis, FALSE, DEFAULT_CONFLICTANALYSIS, nullptr, nullptr) );

   return SCIP_OKAY;
}
