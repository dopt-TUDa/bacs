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

/**@file   prop_simplicial.cpp
 * @ingroup PROPAGATORS
 * @brief  propagator to fix simplicial nodes
 * @author Marc Pfetsch
 *
 * A node v in the graph is simplicial if its neighborhood is a clique. Assume that the graph is unweighted. Then there
 * exists an optimal solution that contains v. We then fix v to 1.
 *
 * There exists a generalization to weighted graphs, see
 * Reductions for the Stable Set Problem@n
 * E. C. Sewell, S. H. Jacobson, Hemanshu Kaul
 * Algorithmic Operations Research Vol.6 (2011) 40–55.
 *
 * We also check for nodes \f$v\f$ of degree two with neighbors \f$u\f$, \f$w\f$. Then the following constraint can be
 * added \f$2 x_v + x_u + x_w = 2\f$.
 *
 * The priority is set such that this propagator runs before prop_neighborhoods, because a simplicial node always would
 * trigger a neighborhood reduction, but fixing the node to 1 is stronger.
 *
 * @todo Implement the generalization.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include <scip/cons_linear.h>
#include <scip/cons_logicor.h>
#include "prop_simplicial.h"
#include "graph.h"
#include "struct_probdata.h"


#define PROP_NAME                  "simplicial"
#define PROP_DESC                  "propagator to fix simplicial nodes"
#define PROP_TIMING                 SCIP_PROPTIMING_BEFORELP
#define PROP_PRIORITY               -800000  //!< priority of the propagator (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PROP_FREQ                        -1  //!< propagation frequency
#define PROP_DELAY                    FALSE  //!< should propagation method be delayed, if other propagators found reductions?
#define PROP_PRESOL_PRIORITY        2000000  //!< priority of the propagator (>= 0: before, < 0: after constraint handlers)
#define PROP_PRESOL_MAXROUNDS             0  //!< maximal number of presolving rounds the propagator participates in (-1: no limit)
#define PROP_PRESOLTIMING           SCIP_PRESOLTIMING_MEDIUM //!< timing of the presolving method (fast, medium, or exhaustive)

#define DEFAULT_CHECKDEGREETWO  FALSE        //!< whether constraints for degree 2 nodes should be added
#define DEFAULT_ADDLOGICOR      FALSE        //!< whether logicor constraints should be added for degree 2 nodes
#define DEFAULT_PROBINGONLY     FALSE        //!< whether the propagator should only run in probing

//! propagator data
struct SCIP_PropData
{
   SCIP_Bool             checkdegreetwo;     //!< whether constraints for degree 2 nodes should be added
   SCIP_Bool             addlogicor;         //!< whether logicor constraints should be added for degree 2 nodes
   SCIP_Bool             probingonly;        //!< whether the propagator should only run in probing
   SCIP_Bool             enabled;            //!< whether the propagator is enabled
};


//! presolving
static
SCIP_RETCODE presolSimplicial(
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
   assert( G != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );

   // prepare data
   size_t n = probdata->n;
   int* cand;
   SCIP_Bool* covered;
   SCIP_CALL( SCIPallocBufferArray(scip, &cand, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &covered, n) );
   for (size_t j = 0; j < n; ++j)
   {
      cand[j] = -1;
      covered[j] = FALSE;
   }

   // loop through nodes
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;

      // skip fixed nodes
      if ( SCIPvarGetUbLocal(probdata->vars[v]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[v]) > 0.5 )
         continue;

      // mark node
      cand[v] = (int) v;
      SCIP_Real weightv = boost::get(vertex_weight_t(), *G, v);

      // determine degree of v w.r.t. nonfixed nodes
      size_t deg = 0;
      bool neighfixedone = false;
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         // if for some reason a neighbor of v is already fixed to 1, we stop
         if ( SCIPvarGetLbLocal(probdata->vars[w]) > 0.5 )
         {
            neighfixedone = true;
            break;
         }

         // skip nodes fixed to 0
         if ( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 )
            continue;

         cand[w] = (int) v;
         ++deg;
      }
      assert( deg == (size_t) probdata->localdegrees[v] );

      // treat case in which neighbor is fixed to 1
      if ( neighfixedone )
      {
         SCIP_Bool tightened = FALSE;

         SCIPdebugMsg(scip, "Fixing node %lu to 0.\n", v);
         SCIP_CALL( SCIPtightenVarUb(scip, probdata->vars[v], 0.0, FALSE, &cutoff, &tightened) );
         assert( ! cutoff );
         assert( tightened );
         ++nfixed;
         continue;
      }

      // loop over neighborhood of node
      size_t ncliquenodes = 0;
      bool weightslarger = false;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         if ( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[w]) > 0.5 )
            continue;

         // mark if the weight of w is larger than that of v
         if ( SCIPisGT(scip, boost::get(vertex_weight_t(), *G, w), weightv) )
            weightslarger = true;

         // count how many neighboring nodes are adjacent to v
         unsigned int nneigh = 0;
         AdjacencyIterator ait2, aend2;
         for (boost::tie(ait2, aend2) = boost::adjacent_vertices(w, *G); ait2 != aend2; ++ait2)
         {
            Vertex x = *ait2;

            if ( SCIPvarGetUbLocal(probdata->vars[x]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[x]) > 0.5 )
               continue;

            if ( cand[x] == (int) v )
            {
               ++nneigh;
               if ( nneigh == deg )
                  break;
            }
         }

         // if the neighbors of w cover all nodes, we have a clique
         if ( nneigh == deg )
         {
            ++ncliquenodes;
            if ( ncliquenodes == deg )
               break;
         }
         else
            break;
      }

      // we found a simplicial node if all neighbors are adjacent to all neighbors
      if ( ncliquenodes == deg )
      {
         SCIP_Bool tightened = FALSE;

         // fix node to 1 if all weights in neighborhood are at most as large
         if ( ! weightslarger )
         {
            SCIPdebugMsg(scip, "Fixing simplicial node %lu to 1.\n", v);
            SCIP_CALL( SCIPtightenVarLb(scip, probdata->vars[v], 1.0, FALSE, &cutoff, &tightened) );
            assert( ! cutoff );
            assert( tightened );
            ++nfixed;
         }

         // explicitly fix all neighbors to 0 in order to avoid a complete round of propagating
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            Vertex w = *ait;

            if ( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[w]) > 0.5 )
               continue;

            SCIP_Real weightw = boost::get(vertex_weight_t(), *G, w);
            if ( SCIPisLE(scip, weightw, weightv) )
            {
               SCIPdebugMsg(scip, "Fixing adjacent node %lu of simplicial node %lu to 0.\n", w, v);
               SCIP_CALL( SCIPtightenVarUb(scip, probdata->vars[w], 0.0, FALSE, &cutoff, &tightened) );
               assert( ! cutoff );
               assert( tightened );
               ++nfixed;
            }
         }
      }
      else
      {
         // special treatment of degree 2 nodes
         if ( deg == 2 && propdata->checkdegreetwo && ! covered[v] )
         {
            SCIPdebugMsg(scip, "Found degree 2 node <%s>.\n", SCIPvarGetName(probdata->vars[v]));
            if ( propdata->addlogicor )
            {
               SCIP_VAR* vars[2];

               vars[0] = probdata->vars[v];
               covered[v] = TRUE;

               int nneigh = 0;
               for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
               {
                  Vertex w = *ait;

                  if ( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[w]) > 0.5 )
                     continue;

                  vars[1] = probdata->vars[w];
                  covered[w] = TRUE;
                  ++nneigh;

                  SCIP_CONS* cons;
                  char name[SCIP_MAXSTRLEN];
                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "degree2#%d#%lu", nneigh, v);
                  SCIP_CALL( SCIPcreateConsBasicLogicor(scip, &cons, name, 2, vars) );
                  SCIP_CALL( SCIPaddCons(scip, cons) );
#ifdef SCIP_DEBUG
                  SCIPprintCons(scip, cons, nullptr);
                  SCIPinfoMessage(scip, nullptr, "\n");
#endif
                  SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               }
               assert( nneigh == 2 );
            }
            else
            {
               SCIP_VAR* vars[3];
               SCIP_Real vals[3];

               vars[0] = probdata->vars[v];
               vals[0] = 2.0;
               covered[v] = TRUE;

               int nvars = 1;
               for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
               {
                  Vertex w = *ait;

                  if ( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[w]) > 0.5 )
                     continue;

                  vars[nvars] = probdata->vars[w];
                  vals[nvars++] = 1.0;
                  covered[w] = TRUE;
               }
               assert( nvars == 3 );

               // generate constraint that either v or the two neighbors of v are in the solution
               SCIP_CONS* cons;
               char name[SCIP_MAXSTRLEN];
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "degree2#%lu", v);
               SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, nvars, vars, vals, 2.0, 2.0) );
               SCIP_CALL( SCIPaddCons(scip, cons) );
#ifdef SCIP_DEBUG
               SCIPprintCons(scip, cons, nullptr);
               SCIPinfoMessage(scip, nullptr, "\n");
#endif
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &covered);
   SCIPfreeBufferArray(scip, &cand);

   return SCIP_OKAY;
}


/** perform propagation
 *
 *  There is some code duplication to presolving above, but we separated the functions for clarity.
 */ /*lint --e{715}*/
static
SCIP_RETCODE propagateSimplicial(
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
   assert( G != nullptr );

   // prepare data
   size_t n = probdata->n;
   int* cand;
   SCIP_CALL( SCIPallocBufferArray(scip, &cand, n) );
   for (size_t j = 0; j < n; ++j)
      cand[j] = -1;

   // loop through nodes
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;

      // skip fixed nodes
      if ( SCIPvarGetUbLocal(probdata->vars[v]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[v]) > 0.5 )
         continue;

      // mark node
      cand[v] = (int) v;
      SCIP_Real weightv = boost::get(vertex_weight_t(), *G, v);

      // determine degree of v w.r.t. nonfixed nodes
      size_t deg = 0;
      bool neighfixedone = false;
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         // if for some reason a neighbor of v is already fixed to 1, we stop
         if ( SCIPvarGetLbLocal(probdata->vars[w]) > 0.5 )
         {
            neighfixedone = true;
            break;
         }

         // skip nodes fixed to 0
         if ( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 )
            continue;

         cand[w] = (int) v;
         ++deg;
      }

      // treat case in which neighbor is fixed to 1
      if ( neighfixedone )
      {
         SCIP_Bool tightened = FALSE;

         SCIPdebugMsg(scip, "Fixing node %lu to 0.\n", v);
         SCIP_CALL( SCIPtightenVarUb(scip, probdata->vars[v], 0.0, FALSE, &cutoff, &tightened) );
#if 0
         SCIP_CALL( SCIPinferBinvarProp(scip, probdata->vars[v], TRUE, prop, (int) v, &cutoff, &tightened) );
#endif
         assert( ! cutoff );
         assert( tightened );
         ++nfixed;
         continue;
      }

      assert( deg == (size_t) probdata->localdegrees[v] );

      // loop over neighborhood of node
      size_t ncliquenodes = 0;
      bool weightslarger = false;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         if ( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[w]) > 0.5 )
            continue;

         // mark if the weight of w is larger than that of v
         if ( SCIPisGT(scip, boost::get(vertex_weight_t(), *G, w), weightv) )
            weightslarger = true;

         // count how many neighboring nodes are adjacent to v
         unsigned int nneigh = 0;
         AdjacencyIterator ait2, aend2;
         for (boost::tie(ait2, aend2) = boost::adjacent_vertices(w, *G); ait2 != aend2; ++ait2)
         {
            Vertex x = *ait2;

            if ( SCIPvarGetUbLocal(probdata->vars[x]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[x]) > 0.5 )
               continue;

            if ( cand[x] == (int) v )
            {
               ++nneigh;
               if ( nneigh == deg )
                  break;
            }
         }

         // if the neighbors of w cover all nodes, we have a clique
         if ( nneigh == deg )
         {
            ++ncliquenodes;
            if ( ncliquenodes == deg )
               break;
         }
         else
            break;
      }

      // we found a simplicial node if all neighbors are adjacent to all neighbors
      if ( ncliquenodes == deg )
      {
         SCIP_Bool tightened = FALSE;

         // fix node to 1 if all weights in neighborhood are at most as large
         if ( ! weightslarger )
         {
            SCIPdebugMsg(scip, "Fixing simplicial node %lu to 1.\n", v);
            SCIP_CALL( SCIPtightenVarLb(scip, probdata->vars[v], 1.0, FALSE, &cutoff, &tightened) );
#if 0
            SCIP_CALL( SCIPinferBinvarProp(scip, probdata->vars[v], TRUE, prop, (int) v, &cutoff, &tightened) );
#endif
            assert( ! cutoff );
            assert( tightened );
            ++nfixed;
         }

         // explicitly fix all neighbors to 0 in order to avoid a complete round of propagating
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            Vertex w = *ait;

            if ( SCIPvarGetUbLocal(probdata->vars[w]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[w]) > 0.5 )
               continue;

            SCIP_Real weightw = boost::get(vertex_weight_t(), *G, w);
            if ( SCIPisLE(scip, weightw, weightv) )
            {
               SCIPdebugMsg(scip, "Fixing adjacent node %lu of simplicial node %lu to 0.\n", w, v);
               SCIP_CALL( SCIPtightenVarUb(scip, probdata->vars[w], 0.0, FALSE, &cutoff, &tightened) );
#if 0
               SCIP_CALL( SCIPinferBinvarProp(scip, probdata->vars[w], FALSE, prop, (int) v, &cutoff, &tightened) );
#endif
               assert( ! cutoff );
               assert( tightened );
               ++nfixed;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &cand);

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */

//! presolving initialization method of propagator (called when presolving is about to begin)
static
SCIP_DECL_PROPINITPRE(propInitpreSimplicial)
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
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "Turned off simplicial presolver, because the problem does not contain exaclty one clique constraint.\n");

   return SCIP_OKAY;
}

//! solving process initialization method of propagator (called when branch and bound process is about to begin)
static
SCIP_DECL_PROPINITSOL(propInitsolSimplicial)
{
   assert( prop != nullptr );

   // if propagation should be performed in the tree, we currently turn off symmetry handling to avoid conflicts
   if ( SCIPpropGetFreq(prop) >= 0 )
   {
      SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
      assert( propdata != nullptr );

      if ( ! propdata->probingonly )
      {
         SCIP_RETCODE retcode = SCIPsetIntParam(scip, "propagating/symmetry/freq", -1);
         if ( retcode == SCIP_OKAY )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, nullptr, "Turning off symmetry handling to avoid conflicts with simplicial propagation.\n");
         else
            assert( retcode == SCIP_PARAMETERUNKNOWN );
      }
   }

   return SCIP_OKAY;
}


//! presolving method
static
SCIP_DECL_PROPPRESOL(propPresolSimplicial)
{
   assert( scip != nullptr );
   assert( result != nullptr );
   assert( nfixedvars != nullptr );

   *result = SCIP_DIDNOTRUN;

   // the propagator performs dual reductions
   if ( ! SCIPallowStrongDualReds(scip) || ! SCIPallowWeakDualReds(scip) )
      return SCIP_OKAY;

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   // currently the propagator only works for unweighted cases
   if ( ! probdata->unweighted )
      return SCIP_OKAY;

   // do not run if disabled
   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );
   if ( ! propdata->enabled )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Running presolving of propagator <%s> ...\n", PROP_NAME);
   SCIP_Bool cutoff;
   int nfixed;
   SCIP_CALL( presolSimplicial(scip, prop, probdata, cutoff, nfixed) );

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
SCIP_DECL_PROPEXEC(propExecSimplicial)
{
   assert( scip != nullptr );
   assert( result != nullptr );

   *result = SCIP_DIDNOTRUN;

   // the propagator performs dual reductions
   if ( ! SCIPallowStrongDualReds(scip) || ! SCIPallowWeakDualReds(scip) )
      return SCIP_OKAY;

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   // currently the propagator only works for unweighted cases
   if ( ! probdata->unweighted )
      return SCIP_OKAY;

   // do not run if disabled
   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );
   if ( ! propdata->enabled )
      return SCIP_OKAY;

   // if propagation should only run in probing
   if ( propdata->probingonly )
   {
      if ( ! SCIPinProbing(scip) || SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING )
      {
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }
   }

   SCIPdebugMsg(scip, "Running propagator <%s> ...\n", PROP_NAME);
   SCIP_Bool cutoff;
   int nfixed;
   SCIP_CALL( propagateSimplicial(scip, prop, probdata, cutoff, nfixed) );

   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( nfixed > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeSimplicial)
{
   assert( scip != nullptr );
   assert( prop != nullptr );

   SCIP_PROPDATA* propdata = SCIPpropGetData(prop);
   assert( propdata != nullptr );
   SCIPfreeBlockMemory(scip, &propdata);


   return SCIP_OKAY;
}



/*
 * propagator specific interface methods
 */

//! creates the simplicial propagator and includes it in SCIP
SCIP_RETCODE BACSincludePropSimplicial(
   SCIP*                 scip                //! SCIP data structure
   )
{
   SCIP_PROPDATA* propdata = nullptr;
   SCIP_PROP* prop = nullptr;

   // create comp data
   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );
   propdata->enabled = TRUE;

   // create simplicial propagator data
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecSimplicial, propdata) );
   assert( prop != nullptr );

   SCIP_CALL( SCIPsetPropInitpre(scip, prop, propInitpreSimplicial) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolSimplicial, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeSimplicial) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolSimplicial) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/checkdegreetwo",
         "whether constraints for degree 2 nodes should be added",
         &propdata->checkdegreetwo, FALSE, DEFAULT_CHECKDEGREETWO, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/addlogicor",
         "whether logicor constraints should be added for degree 2 nodes",
         &propdata->addlogicor, FALSE, DEFAULT_ADDLOGICOR, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/probingonly",
         "whether the propagator should only run in probing",
         &propdata->probingonly, FALSE, DEFAULT_PROBINGONLY, nullptr, nullptr) );

   return SCIP_OKAY;
}
