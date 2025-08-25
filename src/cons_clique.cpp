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
//#define SCIP_DEBUG
/**@file   cons_clique.cpp
 * @brief  constraint handler for clique inequalities
 * @author Marc Pfetsch
 */

#include "cons_clique.h"
#include "vectorpool.hpp"
#include <sys/times.h>
#include <cassert>

#if SCIP_VERSION >= 900
#include "scip/symmetry_graph.h"
#include "symmetry/struct_symmetry.h"
#endif

#define CONSHDLR_NAME          "clique"
#define CONSHDLR_DESC          "clique constraints"
#define CONSHDLR_SEPAPRIORITY    +50000   // after setppc/linear
#define CONSHDLR_ENFOPRIORITY  -1000100   // < 0 - processed after integrality, setppc, linear
#define CONSHDLR_CHECKPRIORITY -1000100   // < 0 - processed after integrality, setppc, linear
#define CONSHDLR_SEPAFREQ             5
#define CONSHDLR_PROPFREQ             1
#define CONSHDLR_EAGERFREQ           -1
#define CONSHDLR_MAXPREROUNDS        -1
#define CONSHDLR_DELAYSEPA        FALSE
#define CONSHDLR_DELAYPROP        FALSE
#define CONSHDLR_NEEDSCONS         TRUE
#define CONSHDLR_PROPTIMING        SCIP_PROPTIMING_ALWAYS
#define CONSHDLR_PRESOLTIMING      SCIP_PRESOLTIMING_MEDIUM

// default values for parameters
#define DEFAULT_MAXNSEPA                     1000    //!< maximal number of clique inequalities separated
#define DEFAULT_MAXNINITIAL                  2000    //!< maximal number of clique inequalities added initially
#define DEFAULT_MAXCLIQUEBBNODES             10000   //!< maximal number of nodes in the clique B&B

#define MAX_CLIQUEPOOL_GEN                   500     //!< maximal number of cuts generated from the pool

#define DEFAULT_INITIALCLIQUES           FALSE       //!< Should initial cliques be generated?
#define DEFAULT_INITIALSELECTION         FALSE       //!< Should initial cliques be selected from an edge covering?
#define DEFAULT_INITIALPOOL              FALSE       //!< Should initial cliques be added to the pool?
#define DEFAULT_EXACTCLIQUES             TRUE        //!< Should exact clique generation be called if no other cuts have been found?
#define DEFAULT_GREEDYCLIQUES            TRUE        //!< Should greedy clique generation be called before exact separation?
#define DEFAULT_EDGEENFORCE              TRUE        //!< Run enforcement based on edges with clique extension?


// parameters for score computation in greedy clique construction
#define COEF_LP_VAL                          1       //!< coeff of current LP value (x) in the value of covering a vertex
#define COEF_OBJ_VAL                         1       //!< coeff of objective value (c) in the value of covering a vertex
#define COEF_ABS_OBJ_VAL                     1       //!< coeff of absolute objective value (i.e. x*c) in the value of covering a vertex
#define COEF_LP_VAL_NEIGH                    1       //!< coeff of current LP value (x) of neighbours in the value of covering a vertex
#define COEF_OBJ_VAL_NEIGH                   1       //!< coeff of objective value (c) of neighbours in the value of covering a vertex
#define COEF_ABS_OBJ_VAL_NEIGH               1       //!< coeff of absolute objective value (i.e. x*c) of neighbours in the value of covering a vertex


//! data for clique constraint handler
struct SCIP_ConshdlrData
{
   unsigned int          ninitialcliques;    //!< number of initially generated cliques
   unsigned int          nsepacliques;       //!< total number of separated clique inequalities
   unsigned int          nsepapoolcliques;   //!< number of clique cuts separated from pool
   unsigned int          nsepaexactcliques;  //!< number of clique cuts separated by exact algorithm
   unsigned int          nsepagreedycliques; //!< number of clique cuts separated by greedy algorithm
   unsigned int          nsepabbcliquerounds;//!< number of separation rounds with B&B clique separation
   unsigned int          nbbnodelimit;       //!< number of times the B&B code hit the node limit
   SCIP_CLOCK*           timepoolsepa;       //!< time for pool separation
   SCIP_CLOCK*           timeexactsepa;      //!< time for exact clique separation
   SCIP_CLOCK*           timegreedysepa;     //!< time for greedy clique separation
   int                   maxnsepa;           //!< maximal number of clique inequalities to be separated
   int                   maxninitial;        //!< maximal number of clique inequalities added initially
   int                   maxcliquebbnodes;   //!< maximal number of nodes in the B&B algorithm to separate cliques
   SCIP_EVENTHDLR*       eventhdlr;          //!< event handler pointer
   SCIP_Bool             initialcliques;     //!< Should initial cliques be generated?
   SCIP_Bool             initialselection;   //!< Should initial cliques be selected from an edge covering?
   SCIP_Bool             initialpool;        //!< Should initial cliques be added to the pool?
   SCIP_Bool             exactcliques;       //!< Should exact clique generation be called if no other cuts have been found?
   SCIP_Bool             greedycliques;      //!< Should greedy clique generation be called before exact separation?
   SCIP_Bool             edgeenforce;        //!< Run enforcement based on edges with clique extension?
};


//! data for clique constraints
struct SCIP_ConsData
{
   const Graph*          G;                  //!< underlying graph
   const std::vector<bool>* isolated;        //!< marks whether a node is isolated
   TCLIQUE_GRAPH*        TCG;                //!< graph for clique separation
   unsigned int          n;                  //!< number of nodes
   SCIP_VAR**            vars;               //!< variables for constraint
   SCIP_Real*            vals;               //!< array for storing LP-solution
   SCIP_Bool             islpsol;            //!< whether the stored solution is an LP-solution
   unsigned int          nvarschanged;       //!< number of changed variables
   SCIP_VAR**            varschanged;        //!< array to store changed variables
   SCIP_Real             scaleclique;        //!< scaling parameter for clique separation
   SCIP_Bool             cutoff;             //!< store whether a cutoff occured in the tclique callback
   vectorpool*           pool;               //!< clique pool
   unsigned int          lastclique;         //!< clique of last pool separation
   SCIP_Bool             trianglefree;       //!< stores whether we know that the graph is triangle-free

   SCIP_VAR**            tmpvars;            //!< temporary storage for variables
   int*                  cliquenodes;        //!< temporary storage for clique
   SCIP*                 scip;               //!< copy of SCIP instance - needed for tclique callback
};




//--------------------------------------------------------------------------------------------
//--------------------------------- event handler --------------------------------------------
//--------------------------------------------------------------------------------------------

#define EVENTHDLR_NAME         "clique"
#define EVENTHDLR_DESC         "report bound changes to clique constraint handler"

//! event execution
SCIP_DECL_EVENTEXEC(BACSeventExecClique)
{
   assert( scip != nullptr );
   assert( eventhdlr != nullptr );
   assert( event != nullptr );
   assert( eventdata != nullptr );

   assert( SCIPeventGetType(event) == SCIP_EVENTTYPE_LBTIGHTENED || SCIPeventGetType(event) == SCIP_EVENTTYPE_GLBCHANGED );

   SCIP_VAR* var = SCIPeventGetVar(event);
   assert( var != nullptr );
   assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );

   // take care of negated variables
   if ( SCIPvarIsNegated(var) )
   {
      var = SCIPvarGetNegatedVar(var);
      assert( var != nullptr );
   }
   assert( ! SCIPvarIsNegated(var) );
   SCIPdebugMsg(scip, "Mark variable <%s> as changed.\n", SCIPvarGetName(var));

   SCIP_CONSDATA* consdata = (SCIP_CONSDATA*) eventdata;
   assert( consdata != nullptr );
   assert( consdata->varschanged != nullptr );
   assert( consdata->nvarschanged <= consdata->n );

   if ( consdata->nvarschanged < consdata->n )
      consdata->varschanged[consdata->nvarschanged++] = var;

   return SCIP_OKAY;
}


//--------------------------------------------------------------------------------------------
//--------------------------------- tclique callback -----------------------------------------
//--------------------------------------------------------------------------------------------

/** callback function to handle a newly found clique
 *
 *  This function is called when the tclique code finds a new clique. We set the parameters such that tclique returns
 *  with a clique if its value is at most 10% below the best clique found so far. This tends to produce enough cliques,
 *  but not too much.
 *
 *  In this function we check whether the corresponding clique inequality is violated and we generate a cut if this is
 *  the case. The cut is entered into the pool of the constraint.
 *
 *  Note that the cuts are only used if they are efficacious (i.e., their violation is larger than the given parameter
 *  (usually around 0.05). Hence, there might be clique (and edge) inequalities, which are slightly violated, but are
 *  not generated.
 */
TCLIQUE_NEWSOL(tcliqueNewsolClique)
{
   assert( acceptsol != nullptr );
   assert( stopsolving != nullptr );
   assert( tcliquedata != nullptr );

   // we do not accept the solution as new incumbent, because we want to find many violated clique inequalities
   *acceptsol = FALSE;
   *stopsolving = FALSE;

   // slightly increase the minimal weight for additional cliques
   TCLIQUE_WEIGHT minweightinc = (cliqueweight - *minweight)/10;
   minweightinc = MAX(minweightinc, 1);
   *minweight += minweightinc;

   // get constraint data
   SCIP_CONS* cons = (SCIP_CONS*) tcliquedata;
   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   assert( consdata != nullptr );
   assert( consdata->n > 0 );
   assert( consdata->vals != nullptr );
   assert( consdata->vars != nullptr );
   assert( consdata->tmpvars != nullptr );
   assert( consdata->scip != nullptr );
   assert( consdata->pool != nullptr );

   // adds cut if weight of the clique is greater than the rhs (1.0 * scaling)
   if ( cliqueweight > consdata->scaleclique )
   {
      SCIP_VAR** vars = consdata->vars;
      SCIP_Real* vals = consdata->vals;
      SCIP* scip = consdata->scip;

      // calculate the weight of the clique in unscaled fractional variable space
      assert( ncliquenodes > 0 );
      SCIP_Real unscaledweight = 0.0;
      for (int i = 0; i < ncliquenodes; ++i)
      {
         assert( 0 <= cliquenodes[i] && cliquenodes[i] < (int) consdata->n );
         unscaledweight += vals[cliquenodes[i]];
      }

      // if the weight of the clique is larger than 1
      if ( SCIPisEfficacious(scip, unscaledweight - 1.0) )
      {
         std::vector<unsigned int> C;
         C.reserve((unsigned) ncliquenodes);

         SCIP_VAR** tmpvars = consdata->tmpvars;
         for (int i = 0; i < ncliquenodes; ++i)
         {
            tmpvars[i] = vars[cliquenodes[i]];
            C.push_back((unsigned int) cliquenodes[i]);
         }

         SCIP_CONSHDLR* conshdlr = SCIPconsGetHdlr(cons);
         assert( conshdlr != nullptr );
         SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert( conshdlrdata != nullptr );

         SCIP_ROW* row;
         SCIP_RETCODE retcode;
#ifndef NDEBUG
         char name[SCIP_MAXSTRLEN];
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "clique%u", conshdlrdata->nsepacliques);
         retcode = SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE);
#else
         retcode = SCIPcreateEmptyRowCons(scip, &row, cons, "clique", -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE);
#endif

         if ( retcode != SCIP_OKAY )
         {
            *stopsolving = TRUE;
            return;
         }
         retcode = SCIPaddVarsToRowSameCoef(scip, row, ncliquenodes, tmpvars, 1.0);
         if ( retcode != SCIP_OKAY )
         {
            (void) SCIPreleaseRow(scip, &row);
            *stopsolving = TRUE;
            return;
         }

#ifdef SHOW_CLIQUES
         retcode =  SCIPprintRow(scip, row, nullptr);
         if ( retcode != SCIP_OKAY )
         {
            (void) SCIPreleaseRow(scip, &row);
            *stopsolving = TRUE;
            return;
         }
#endif
         // note: the cut may not be violated if we are separating a solution different from the LP solution!
         if ( consdata->islpsol || SCIPisFeasNegative(scip, SCIPgetRowLPFeasibility(scip, row)) )
         {
            SCIP_Bool cutoff;
            retcode = SCIPaddRow(scip, row, FALSE, &cutoff);
            if ( retcode != SCIP_OKAY )
            {
               (void) SCIPreleaseRow(scip, &row);
               *stopsolving = TRUE;
               return;
            }

            // mark whether we have a cutoff
            if ( cutoff )
            {
               consdata->cutoff = TRUE;
               *stopsolving = TRUE;
            }

            ++conshdlrdata->nsepacliques;
         }

         retcode = SCIPreleaseRow(scip, &row);
         if ( retcode != SCIP_OKAY )
         {
            *stopsolving = TRUE;
            return;
         }

         // only insert cliques of size at least three
         if ( C.size() > 2 )
            (void) consdata->pool->insert(C);

         return;
      }
   }
}


/** create tclique graph */
static
SCIP_RETCODE computeTCliqueGraph(
   SCIP*                 scip,               //!< SCIP pointer
   SCIP_CONSDATA*        consdata            //!< constraint data
   )
{
   assert( scip != nullptr );
   assert( consdata != nullptr );

   // init TCG
   if ( ! tcliqueCreate(&consdata->TCG) )
      return SCIP_NOMEMORY;

   // create nodes of tclique-graph and initialize numbers for each node
   const Graph* G = consdata->G;
   int nodenumber = 0;
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      if ( ! tcliqueAddNode(consdata->TCG, nodenumber++, 0) )
      {
         SCIPABORT();
         return SCIP_INVALIDDATA;
      }
   }

   // create edges
   EdgeIterator eit, eend;
   for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
   {
      Vertex s = boost::source(*eit, *G);
      Vertex t = boost::target(*eit, *G);

      assert( s < boost::num_vertices(*G) );
      assert( t < boost::num_vertices(*G) );

      if ( ! tcliqueAddEdge(consdata->TCG, (int) s, (int) t) )
      {
         SCIPABORT();
         return SCIP_INVALIDDATA;
      }
   }
   if ( ! tcliqueFlush(consdata->TCG) )
   {
      SCIPABORT();
      return SCIP_INVALIDDATA;
   }

   assert( tcliqueGetNEdges(consdata->TCG) == (int) (2 * boost::num_edges(*G)) );
   assert( tcliqueGetNNodes(consdata->TCG) == (int) boost::num_vertices(*G) );

   return SCIP_OKAY;
}

//--------------------------------------------------------------------------------------------
//------------------------------------- local functions --------------------------------------
//--------------------------------------------------------------------------------------------


/** get values of SCIP solution @p sol
 *
 *  We store the solution in the data of the constraint.
 *
 *  The values are stored in the constraint because of efficiency reasons: During pool separation the solution values of
 *  variables are often accessed more than one time. Furthermore, if the exact separation method is executed, the
 *  variables are also accessed.
 */
static
SCIP_RETCODE getSolValues(
   SCIP*                 scip,               //!< SCIP pointer
   SCIP_SOL*             sol,                //!< solution (or nullptr)
   SCIP_CONSDATA*        consdata            //!< constraint data
   )
{
   assert( scip != nullptr );
   assert( consdata != nullptr );
   assert( consdata->vars != nullptr );
   assert( consdata->vals != nullptr );

   SCIP_CALL( SCIPgetSolVals(scip, sol, (int) consdata->n, consdata->vars, consdata->vals) );

   if ( sol == nullptr )
      consdata->islpsol = TRUE;
   else
      consdata->islpsol = FALSE;

   return SCIP_OKAY;
}


/** set initial solution values
 *
 *  All values are set to 1, except for fixed variables.
 */
static
SCIP_RETCODE setupInitialSolVals(
   SCIP*                 scip,               //!< SCIP pointer
   SCIP_CONSDATA*        consdata            //!< constraint data
   )
{
   assert( scip != nullptr );
   assert( consdata != nullptr );
   assert( consdata->vars != nullptr );
   assert( consdata->vals != nullptr );

   unsigned int n = consdata->n;
   SCIP_VAR** vars = consdata->vars;
   SCIP_Real* vals = consdata->vals;

   for (unsigned int i = 0; i < n; ++i)
   {
      // set values to 0 for fixed variables (cliques with vars fixed to 1 are always fulfilled)
      if ( SCIPvarGetUbLocal(vars[i]) < 0.5 || SCIPvarGetLbLocal(vars[i]) > 0.5 )
         vals[i] = 0.0;
      else
         vals[i] = 1.0;
   }
   consdata->islpsol = FALSE;

   return SCIP_OKAY;
}


/** Scoring function approximating value of covering a node */
static
SCIP_Real score(
   const Graph*          G,                  //!< graph
   SCIP_Real*            vals,               //!< values of nodes
   Vertex                v                   //!< node to be scored
   )
{
   SCIP_Real val = vals[v];
   SCIP_Real weight = boost::get(vertex_weight_t(), *G, v);

   SCIP_Real self_val = COEF_LP_VAL * val + COEF_OBJ_VAL * weight + COEF_ABS_OBJ_VAL * weight * val;

   AdjacencyIterator ait, aend;
   SCIP_Real neighbor_val = 1e-3;
   for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
   {
      val = vals[*ait];
      weight = boost::get(vertex_weight_t(), *G, *ait);

      neighbor_val += COEF_LP_VAL_NEIGH * val + COEF_OBJ_VAL_NEIGH * weight + COEF_ABS_OBJ_VAL_NEIGH * weight * val;
   }
   return self_val * (3.0 + std::max(neighbor_val, 6.0));
}


/** generate greedy cliques
 *
 *  We use a greedy heuristic to compute a clique covering and insert the corresponding cliques into the pool.
 */
static
SCIP_RETCODE greedyCliquesHeur(
   SCIP*                 scip,               //!< SCIP pointer
   SCIP_CONSHDLRDATA*    conshdlrdata,       //!< constraint handler data
   SCIP_CONS*            cons,               //!< clique constraint
   unsigned int&         ngen,               //!< number of generated cliques
   SCIP_Bool&            cutoff              //!< has clique been found, that cuts off node
   )
{
   assert( scip != nullptr );
   assert( conshdlrdata != nullptr );
   assert( cons != nullptr );

   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   assert( consdata != nullptr );
   assert( consdata->cliquenodes );
   assert( consdata->pool != nullptr );

   ngen = 0;
   cutoff = FALSE;

   // nothing to do if there are no nodes
   unsigned int n = consdata->n;
   if ( n <= 1 )
      return SCIP_OKAY;

   const Graph* G = consdata->G;
   vectorpool* pool = consdata->pool;

   // sort nodes
   SCIP_Real* scores;
   int* nodes;       // use int because we sort using SCIP functions below
   int* neighbors;   // use int because we sort using SCIP functions below
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &neighbors, n) );
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      scores[v] = score(consdata->G, consdata->vals, v);
      nodes[v] = (int) v;
   }
   SCIPsortDownRealInt(scores, nodes, (int) n);

   // prepare loop
   std::vector<int> cand(n, -1);
   std::vector<unsigned int> C;
   C.reserve(n);

   // loop through sorted nodes
   for (unsigned int k = 0; k < n; ++k)
   {
      assert( nodes[k] >= 0 );
      Vertex v = (Vertex) nodes[k]; /*lint !e571*/

      // skip nodes already covered
      if ( cand[v] >= 0 )
         continue;
      cand[v] = (int) v;

      SCIP_Real cval = consdata->vals[v];

      // skip nodes with 0 LP value
      if ( SCIPisFeasZero(scip, cval) )
         continue;

      // init clique
      C.clear();
      C.push_back((unsigned int) v);
      std::size_t csize = 1;

      // sorted adjacency list
      size_t nneighbors = 0;
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;
         if ( cand[w] >= 0 )
            continue;
         scores[nneighbors] = score(G, consdata->vals, w);
         neighbors[nneighbors++] = (int) w;
      }

      if ( nneighbors == 0 )
         continue;

      SCIPsortDownRealInt(scores, neighbors, (int) nneighbors);

      for(size_t i = 0; i < nneighbors; i++)
      {
         assert( neighbors[i] >= 0 );
         Vertex w = (Vertex) neighbors[i]; /*lint !e571*/
         assert( cand[w] < 0 );

         if ( csize == 1 )
         {
            C.push_back((unsigned int) w);
            ++csize;
            cval += consdata->vals[w];
            cand[w] = (int) v;
            continue;
         }

         // check whether we can extend the clique if we are not triangle-free
         if ( ! consdata->trianglefree )
         {
            unsigned int nneigh = 0;
            for (boost::tie(ait, aend) = boost::adjacent_vertices(w, *G); ait != aend; ++ait)
            {
               Vertex x = *ait;
               if ( cand[x] == (int) v )
               {
                  ++nneigh;
                  if ( nneigh == csize )
                     break;
               }
            }

            // if the neighbors of w cover all nodes in the current clique, we can extend it
            if ( nneigh == csize )
            {
               C.push_back((unsigned int) w);
               ++csize;
               cval += consdata->vals[w];
               cand[w] = (int) v;
            }
         }
      }

      // post processing for maximality of cliques if we are not triangle-free
      if ( ! consdata->trianglefree )
      {
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            Vertex w = *ait;

            // already included in clique
            if ( cand[w] == (int) v )
               continue;

            // check whether we can extend the clique
            unsigned int nneigh = 0;
            AdjacencyIterator ait2, aend2;
            for (boost::tie(ait2, aend2) = boost::adjacent_vertices(w, *G); ait2 != aend2; ++ait2)
            {
               Vertex x = *ait2;

               if ( cand[x] == (int) v )
               {
                  ++nneigh;
                  if ( nneigh == csize )
                     break;
               }
            }

            // if the neighbors of w cover all nodes in the current clique, we can extend it
            if ( nneigh == csize )
            {
               C.push_back((unsigned int) w);
               ++csize;
               cval += consdata->vals[w];
               cand[w] = (int) v;
            }
         }
      }

      // insert cliques of size larger than 2 into pool and possibly add rows
      if ( SCIPisEfficacious(scip, cval - 1.0) && csize >= 2 )
      {
         if ( csize > 2 )
            (void) pool->insert(C);
         ++ngen;

         SCIP_VAR** tmpvars = consdata->tmpvars;
         for (size_t i = 0; i < csize; ++i)
         {
            v = C[i];
            tmpvars[i] = consdata->vars[v];
         }

         // create row and add to SCIP
         SCIP_ROW* row;
#ifndef NDEBUG
         char name[SCIP_MAXSTRLEN];
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "clique%u", conshdlrdata->nsepacliques);
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
#else
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "clique", -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
#endif
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, row, (int) csize, tmpvars, 1.0) );

#ifdef SHOW_CLIQUES
         SCIP_CALL( SCIPprintRow(scip, row, nullptr) );
#endif

         // note: the cut may not be violated if we are separating a solution different from the LP solution!
         if ( consdata->islpsol || SCIPisFeasNegative(scip, SCIPgetRowLPFeasibility(scip, row)) )
         {
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &cutoff) );
            ++conshdlrdata->nsepacliques;
         }
         SCIP_CALL( SCIPreleaseRow(scip, &row) );

         if ( cutoff )
            break;
      }
   }
   SCIPfreeBufferArray(scip, &neighbors);
   SCIPfreeBufferArray(scip, &nodes);
   SCIPfreeBufferArray(scip, &scores);

   return SCIP_OKAY;
}


/** separate clique inequalities from the pool
 *
 *  We check for all cliques in the pool whether they lead to violated inequalities.
 */
static
SCIP_RETCODE separatePoolCliques(
   SCIP*                 scip,               //!< SCIP pointer
   SCIP_CONSHDLRDATA*    conshdlrdata,       //!< constraint handler data
   SCIP_CONS*            cons,               //!< constraint
   unsigned int          maxgen,             //!< maximal number of inequalities to generate
   unsigned int&         ngen,               //!< number of generated clique inequalities
   SCIP_Bool&            cutoff              //!< return whether a cutoff has been detected
   )
{
   assert( scip != nullptr );
   assert( cons != nullptr );
   assert( conshdlrdata != nullptr );

   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   assert( consdata != nullptr );
   assert( consdata->vars != nullptr );
   assert( consdata->vals != nullptr );
   assert( consdata->tmpvars != nullptr );
   assert( consdata->pool != nullptr );

   SCIP_VAR** vars = consdata->vars;
   SCIP_Real* vals = consdata->vals;
   SCIP_VAR** tmpvars = consdata->tmpvars;
   vectorpool* pool = consdata->pool;

   ngen = 0;
   cutoff = FALSE;
   SCIP_CALL( SCIPstartClock(scip, conshdlrdata->timepoolsepa) );

   SCIP_Real minefficacy = SCIPgetSepaMinEfficacy(scip);

   // perform first part of the actual pool separation part (last clique to end)
   unsigned int size = (unsigned int) pool->size();
   assert( size == 0 || consdata->lastclique < size );
   for (unsigned int i = consdata->lastclique; i < size; ++i)
   {
      const std::vector<unsigned int>* C = (*pool)[i];

      // there should be no clique of size 2 in the pool....
      assert( (*C).size() > 2 );

      SCIP_Real value = -1.0;
      std::vector<unsigned int>::const_iterator vit, vend = C->end();
      for (vit = C->begin(); vit != vend; ++vit)
      {
         assert( *vit < consdata->n );
         value += vals[*vit];
         if ( value > minefficacy ) // equivalent to SCIPisEfficacious(scip, value)
            break;
      }

      // if the weight of the clique is larger than 1.0
      if ( value > minefficacy )  // equivalent to SCIPisEfficacious(scip, value)
      {
         int nvars = 0;
         for (vit = C->begin(); vit != vend; ++vit)
         {
            assert( *vit < consdata->n );
            tmpvars[nvars++] = vars[*vit];
         }

         SCIP_ROW* row;
#ifndef NDEBUG
         char name[SCIP_MAXSTRLEN];
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "poolclique%u", conshdlrdata->nsepapoolcliques);
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
#else
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "poolclique", -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
#endif
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, row, nvars, tmpvars, 1.0) );

#ifdef SHOW_CLIQUES
         SCIP_CALL( SCIPprintRow(scip, row, nullptr) );
#endif

         // note: cut may not be violated if we are separating a solution different from the LP solution!
         if ( consdata->islpsol || SCIPisFeasNegative(scip, SCIPgetRowLPFeasibility(scip, row)) )
         {
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &cutoff) );
            ++conshdlrdata->nsepapoolcliques;
            ++ngen;
         }

         SCIP_CALL( SCIPreleaseRow(scip, &row) );

         if ( ngen >= maxgen || cutoff )
         {
            consdata->lastclique = i;
            SCIP_CALL( SCIPstopClock(scip, conshdlrdata->timepoolsepa) );
            return SCIP_OKAY;
         }
      }
   }

   // second part (beginning to last clique)
   for (unsigned int i = 0; i < consdata->lastclique; ++i)
   {
      const std::vector<unsigned int>* C = (*pool)[i];

      // there should be no clique of size 2 in the pool....
      assert( (*C).size() > 2 );

      SCIP_Real value = -1.0;
      std::vector<unsigned int>::const_iterator vit, vend = C->end();
      for (vit = C->begin(); vit != vend; ++vit)
      {
         assert( *vit < consdata->n );
         value += vals[*vit];
         if ( value > minefficacy ) // equivalent to SCIPisEfficacious(scip, value)
            break;
      }

      // if the weight of the clique is larger than 1.0
      if ( value > minefficacy )  // equivalent to SCIPisEfficacious(scip, value)
      {
         int nvars = 0;
         for (vit = C->begin(); vit != vend; ++vit)
         {
            assert( *vit < consdata->n );
            tmpvars[nvars++] = vars[*vit];
         }

         SCIP_ROW* row;
#ifndef NDEBUG
         char name[SCIP_MAXSTRLEN];
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "poolclique%u", conshdlrdata->nsepapoolcliques);
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
#else
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "poolclique", -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
#endif
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, row, nvars, tmpvars, 1.0) );

#ifdef SHOW_CLIQUES
         SCIP_CALL( SCIPprintRow(scip, row, nullptr) );
#endif

         // note: cut may not be violated if we are separating a solution different from the LP solution!
         if ( consdata->islpsol || SCIPisFeasNegative(scip, SCIPgetRowLPFeasibility(scip, row)) )
         {
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &cutoff) );
            ++conshdlrdata->nsepapoolcliques;
            ++ngen;
         }

         SCIP_CALL( SCIPreleaseRow(scip, &row) );

         if ( ngen >= maxgen || cutoff )
         {
            consdata->lastclique = i;
            SCIP_CALL( SCIPstopClock(scip, conshdlrdata->timepoolsepa) );
            return SCIP_OKAY;
         }
      }
   }
   SCIP_CALL( SCIPstopClock(scip, conshdlrdata->timepoolsepa) );
   consdata->lastclique = 0;

   return SCIP_OKAY;
}

/** generate exact cliques
 *
 * We call the function @c tcliqueMaxClique() (included in SCIP), which computes a maximum clique. We need to pass
 *  functions to access the graph. Here we use the functions already implemented in SCIP, except @c
 *  tcliqueNewsolClique() defined above.
 *
 *  The last four parameters are:
 *  - @c TCLIQUE_WEIGHT @c maxfirstnodeweight,  maximum weight of branching nodes in level 0; 0 if not used
 *                                              for cliques with at least one fractional node
 *  - @c TCLIQUE_WEIGHT @c minweight,           lower bound for weight of generated cliques
 *  - @c int            @c maxntreenodes,       maximal number of nodes of B&B tree
 *  - @c int            @c maxnzeroextensions   maximal number of zero-valued variables extending the clique
 *
 *  - We set @c minweight to @c scaleClique-1, since we want the clique to contain at least one fractional variable,
 *    and the maximal alue a fractional variable can have is exactly this value.
 *  - We set @c minweight to be @c scaleClique+1, since the size should be larger than 1.0.
 *  - The maximal number of tree nodes is @c maxCliqueIter (a large running time for some instances can result if this
 *    is too large).
 *  - The maximal number of variables extending a clique is 1000.
 */
static
SCIP_RETCODE exactCliques(
   SCIP*                 scip,               //!< SCIP pointer
   SCIP_CONSHDLRDATA*    conshdlrdata,       //!< constraint handler data
   SCIP_CONS*            cons,               //!< constraint
   unsigned int&         ngen,               //!< number of generated cuts
   SCIP_Bool&            cutoff              //!< store whether a cutoff occured
   )
{
   assert( scip != nullptr );
   assert( conshdlrdata != nullptr );

   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   assert( consdata != nullptr );
   assert( consdata->vals != nullptr );

   if ( consdata->TCG == nullptr )
   {
      assert( ! SCIPconsIsOriginal(cons) );
      SCIP_CALL( computeTCliqueGraph(scip, consdata) );
   }
   assert( consdata->TCG != nullptr );

   cutoff = FALSE;

   unsigned int n = consdata->n;
   SCIP_Real* vals = consdata->vals;
   TCLIQUE_GRAPH* TCG = consdata->TCG;

   // store previous number of separated cliques
   unsigned int noldsepacliques = conshdlrdata->nsepacliques;

   // initialize weights
   for (unsigned int i = 0; i < n; ++i)
   {
      // max with 0, because of numerical issues
      TCLIQUE_WEIGHT weight = MAX((TCLIQUE_WEIGHT) SCIPfeasFloor(scip, vals[i] * consdata->scaleclique), 0);
      tcliqueChangeWeight(TCG, (int) i, weight);
   }

   // call the clique algorithm ...
   int ncliquenodes = 0;
   TCLIQUE_WEIGHT weight = 0;
   TCLIQUE_STATUS status;
   consdata->cutoff = FALSE;

   // 3rd to last parameter: allow to generate cliques containing 0 weight nodes
   tcliqueMaxClique(tcliqueGetNNodes, tcliqueGetWeights, tcliqueIsEdge, tcliqueSelectAdjnodes,
      TCG, tcliqueNewsolClique, (TCLIQUE_DATA*) cons, consdata->cliquenodes, &ncliquenodes,
      &weight, (TCLIQUE_WEIGHT) consdata->scaleclique-1, (TCLIQUE_WEIGHT) consdata->scaleclique + 1,
      conshdlrdata->maxcliquebbnodes, 0, (int) n, -1, nullptr, &status);

   ++conshdlrdata->nsepabbcliquerounds;
   if ( status == TCLIQUE_NODELIMIT )
      ++conshdlrdata->nbbnodelimit;

   if ( consdata->cutoff )
      cutoff = TRUE;

   ngen = conshdlrdata->nsepacliques - noldsepacliques;

   return SCIP_OKAY;
}

/** separate clique inequalities
 *
 *  We call the function @c greedyCliquesHeur() if the greedycliques parameter is set to TRUE. If the greedy algorithm has
 *  not generated any cliques, the current node has not been cut off and exactcliques are activated, the @c exactCliques()
 *  function is called to generate cliques. 
 */
static
SCIP_RETCODE separateCliques(
   SCIP*                 scip,               //!< SCIP pointer
   SCIP_CONSHDLRDATA*    conshdlrdata,       //!< constraint handler data
   SCIP_CONS*            cons,               //!< constraint
   unsigned int&         ngen,               //!< number of generated cuts
   SCIP_Bool&            cutoff              //!< store whether a cutoff occured
   )
{
   assert( scip != nullptr );
   assert( conshdlrdata != nullptr );

   ngen = 0;
   cutoff = FALSE;

   if ( conshdlrdata->greedycliques )
   {
      SCIP_CALL( SCIPstartClock(scip, conshdlrdata->timegreedysepa) );
      SCIP_CALL( greedyCliquesHeur(scip, conshdlrdata, cons, ngen, cutoff) );
      SCIP_CALL( SCIPstopClock(scip, conshdlrdata->timegreedysepa) );
      conshdlrdata->nsepagreedycliques += ngen;
   }

   if ( ! cutoff && ngen == 0 && conshdlrdata->exactcliques )
   {
      SCIP_CALL( SCIPstartClock(scip, conshdlrdata->timeexactsepa) );
      SCIP_CALL( exactCliques(scip, conshdlrdata, cons, ngen, cutoff) );
      conshdlrdata->nsepaexactcliques += ngen;
      SCIP_CALL( SCIPstopClock(scip, conshdlrdata->timeexactsepa) );
   }

   return SCIP_OKAY;
}


/** enforce integral LP solution
 *
 *  We loop through the edges of the graph. If the solution is infeasible at least one edge must be violated. Then we
 *  possibly extend the edge to a maximal clique.
 */
static
SCIP_RETCODE enforceLPSol(
   SCIP*                 scip,               //!< SCIP pointer
   SCIP_CONSHDLRDATA*    conshdlrdata,       //!< constraint handler data
   SCIP_CONS*            cons,               //!< clique constraint
   SCIP_Bool             initial,            //!< whether we are creating initial cliques
   unsigned int&         ngen,               //!< number of generated cliques
   SCIP_Bool&            cutoff              //!< has clique been found, that cuts off node
   )
{
   assert( scip != nullptr );
   assert( conshdlrdata != nullptr );
   assert( cons != nullptr );

   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   assert( consdata != nullptr );

   ngen = 0;
   cutoff = FALSE;

   // nothing to do if there are no nodes
   unsigned int n = consdata->n;
   if ( n == 0 )
      return SCIP_OKAY;

   const Graph* G = consdata->G;
   SCIP_Real* vals = consdata->vals;
   SCIP_VAR** tmpvars = consdata->tmpvars;
   assert( vals != nullptr );
   assert( tmpvars != nullptr );

   // prepare clique marker
   std::vector<int> cand(n, -1);

   // loop through all edges
   SCIP_Bool extended = FALSE;
   size_t edgenum = 0;
   EdgeIterator eit, eend;
   for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
   {
      Edge e = *eit;
      Vertex source = boost::source(e, *G);
      Vertex target = boost::target(e, *G);
      assert( SCIPisFeasIntegral(scip, vals[source]) );
      assert( SCIPisFeasIntegral(scip, vals[target]) );

      // skip edge if nodes have already been covered by a clique
      if ( cand[source] >= 0 && cand[target] >= 0 )
         continue;

      // if the weight of the edge is larger than 1.0
      if ( vals[source] + vals[target] > 1.5 )
      {
         // initialize clique
         cand[source] = (int) edgenum;
         cand[target] = (int) edgenum;
         size_t csize = 2;
         tmpvars[0] = consdata->vars[source];
         tmpvars[1] = consdata->vars[target];

         // post processing for possibly extending clique if we are not triangle-free
         if ( ! consdata->trianglefree )
         {
            AdjacencyIterator ait, aend;
            for (boost::tie(ait, aend) = boost::adjacent_vertices(source, *G); ait != aend; ++ait)
            {
               Vertex w = *ait;

               // if already included in clique
               if ( cand[w] == (int) edgenum )
                  continue;

               // check whether we can extend the clique
               unsigned int nneigh = 0;
               AdjacencyIterator ait2, aend2;
               for (boost::tie(ait2, aend2) = boost::adjacent_vertices(w, *G); ait2 != aend2; ++ait2)
               {
                  Vertex x = *ait2;

                  if ( cand[x] == (int) edgenum )
                  {
                     ++nneigh;
                     if ( nneigh == csize )
                        break;
                  }
               }

               // if the neighbors of w cover all nodes in the current clique, we can extend it
               if ( nneigh == csize )
               {
                  assert( SCIPisFeasIntegral(scip, vals[w]) );
                  tmpvars[csize++] = consdata->vars[w];
                  cand[w] = (int) edgenum;
                  extended = TRUE;
               }
            }
         }

         // create row and add to SCIP
         SCIP_ROW* row;
#ifndef NDEBUG
         char name[SCIP_MAXSTRLEN];
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "clique%u", conshdlrdata->nsepacliques);
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
#else
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "clique", -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
#endif
         ++conshdlrdata->nsepacliques;

         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, row, (int) csize, tmpvars, 1.0) );
         SCIP_CALL( SCIPaddRow(scip, row, FALSE, &cutoff) );
#if 0
         SCIP_CALL( SCIPprintRow(scip, row, nullptr) );
#endif
         // cut should be violated:
         assert( ! consdata->islpsol || SCIPisFeasNegative(scip, SCIPgetRowLPFeasibility(scip, row)) );
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
         ++ngen;

         if ( cutoff )
            return SCIP_OKAY;

         ++edgenum;
      }
   }

   // if we are running initially and could not extend the clique, we are triangle-free
   if ( initial && ! extended )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Graph of constraint <%s> is triangle free.\n", SCIPconsGetName(cons));
      consdata->trianglefree = TRUE;
   }

   return SCIP_OKAY;
}


/** ensure size of cliquenodes array */
static
SCIP_RETCODE ensureCliquenodesSize(
   SCIP*                 scip,               //!< SCIP pointer
   unsigned int**        cliquenodes,        /**< array to be checked */
   int                   cliquenodessize,    /**< size to be ensured */
   int&                  cliquenodesmaxsize  /**< current size of cliquenodes */
   )
{
   assert( cliquenodes != nullptr );

   if ( cliquenodessize > cliquenodesmaxsize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, cliquenodessize);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, cliquenodes, cliquenodesmaxsize, newsize) );
      cliquenodesmaxsize = newsize;
   }
   return SCIP_OKAY;
}

/** compute initial clique cuts to add
 *
 * - We first compute a covering of the edges by cliques. This is done by possibly extending each edge of the graph to a
 *   maximal clique.
 * - Then we run a subgradient algorithm to compute an approximation of the corresponding LP-relaxation with the
 *   produced cliques. For this we take the Lagrangian relaxation of the cuts.
 * - Then we sort the dual values and create clique cuts for the largest values.
 *
 * The subgradient method works as follows. Let \f$Ax \leq b\f$ be the clique inequalities, i.e., \f$b\f$ is the
 * all-ones vector. Then the Lagrangian function is \f$L(\lambda) = c^T x + \lambda^T(b - A x) = (c^T - \lambda^T A) x +
 * \lambda^T b\f$. For its computation with given \f$\lambda\f$, one sets \f$x_j^*\f$ to 1 if the corresponding
 * objective value is positive and 0 otherwise. The corresponding subgradient is then \f$h^k = b - Ax^*\f$, which is the
 * slack (violation) of the constraints. We also use stabilization, i.e., the update formula is \f$\lambda^{k+1} =
 * \lambda^k + \sigma_k(\alpha h^K + (1 - \alpha) h^{k-1})\f$, where \f$\sigma_k\f$ is the step size and \f$\alpha \in
 * [0,1]\f$ is a convexity parameter. The step size is strongly decreased if the objective increases for several
 * iterations and is slowly decreased otherwise.
 */
static
SCIP_RETCODE generateInitialCliqueCuts(
   SCIP*                 scip,               //!< SCIP pointer
   SCIP_CONSHDLRDATA*    conshdlrdata,       //!< constraint handler data
   SCIP_CONS*            cons,               //!< clique constraint
   unsigned int&         ngen,               //!< number of generated cliques
   SCIP_Bool&            cutoff              //!< whether a clique has been found that cuts node off
   )
{
   assert( scip != nullptr );
   assert( conshdlrdata != nullptr );
   assert( cons != nullptr );

   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   assert( consdata != nullptr );

   cutoff = FALSE;
   ngen = 0;

   // nothing to do if there are no nodes
   unsigned int n = consdata->n;
   if ( n <= 1 )
      return SCIP_OKAY;

   // get graph
   const Graph* G = consdata->G;
   assert( G != nullptr );

   // estimate size needed for storing cliques
   int m = (int) boost::num_edges(*G);
   int cliquenodesmaxsize = m * 7;  // 7 is just an arbitrary guess for the average size of a clique

   // prepare space for cliques
   unsigned int* cliquenodes;   // use ints to save space
   int* begcliques;
   size_t numcliques = 0;
   int cliquenodessize = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cliquenodes, cliquenodesmaxsize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &begcliques, m + 1) );  // +1 to encode the end

   // prepare clique marker
   size_t edgenum = 0;
   std::vector<int> cand(n, -1);

   // first compute edge covering by cliques - this might be too large
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Computing cliques covering all edges ...\n");
   SCIP_Bool extended = FALSE;
   EdgeIterator eit, eend;
   for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
   {
      Edge e = *eit;
      Vertex source = boost::source(e, *G);
      Vertex target = boost::target(e, *G);

      // skip edge if it has already been part of a clique (if cand is equal for both end nodes)
      if ( cand[source] >= 0 && cand[source] == cand[target] )
         continue;

      // skip edge if one of its variables is fixed
      SCIP_VAR* sourcevar = consdata->vars[source];
      SCIP_VAR* targetvar = consdata->vars[target];
      if ( SCIPvarGetLbLocal(sourcevar) > 0.5 || SCIPvarGetUbLocal(sourcevar) < 0.5 )
         continue;
      if ( SCIPvarGetLbLocal(targetvar) > 0.5 || SCIPvarGetUbLocal(targetvar) < 0.5 )
         continue;

      // initialize marker for clique
      cand[source] = (int) edgenum;
      cand[target] = (int) edgenum;

      // start new clique
      SCIP_CALL( ensureCliquenodesSize(scip, &cliquenodes, cliquenodessize + 2, cliquenodesmaxsize) );
      begcliques[numcliques++] = cliquenodessize;
      cliquenodes[cliquenodessize++] = (unsigned int) source;
      cliquenodes[cliquenodessize++] = (unsigned int) target;
      size_t csize = 2;

      // check whether we can extend the clique
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(source, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         // if already included in clique
         if ( cand[w] == (int) edgenum )
            continue;

         // check whether we can extend the clique
         unsigned int nneigh = 0;
         AdjacencyIterator ait2, aend2;
         for (boost::tie(ait2, aend2) = boost::adjacent_vertices(w, *G); ait2 != aend2; ++ait2)
         {
            Vertex x = *ait2;

            if ( cand[x] == (int) edgenum )
            {
               ++nneigh;
               if ( nneigh == csize )
                  break;
            }
         }

         // if the neighbors of w cover all nodes in the current clique, we can extend it
         if ( nneigh == csize )
         {
            SCIP_CALL( ensureCliquenodesSize(scip, &cliquenodes, cliquenodessize + 1, cliquenodesmaxsize) );
            cliquenodes[cliquenodessize++] = (unsigned int) w;
            ++csize;
            cand[w] = (int) edgenum;
            extended = TRUE;
         }
      }
      ++edgenum;
   }

   // finish array
   begcliques[numcliques] = cliquenodessize;

   // if we could not extend the clique, we are triangle-free
   if ( ! extended )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Graph of constraint <%s> is triangle free.\n", SCIPconsGetName(cons));
      consdata->trianglefree = TRUE;
   }

#ifdef SCIP_MORE_OUTPUT
   // print found cliques
   for (size_t c = 0; c < numcliques; ++c)
   {
      SCIPinfoMessage(scip, nullptr, "clique %lu: ", c);
      for (int j = begcliques[c]; j < begcliques[c+1]; ++j)
         SCIPinfoMessage(scip, nullptr, "%d ", cliquenodes[j]);
      SCIPinfoMessage(scip, nullptr, "\n");
   }
#endif
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Found %lu cliques covering the edges.\n", numcliques);

   // prepare dual solution of Lagranian problem, initially 0
   SCIP_Real* solution;
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &solution, numcliques) );

   // run subgradient method if we have to many cliques
   if ( (int) numcliques > conshdlrdata->maxninitial )
   {
      // allocate space
      SCIP_Real* oldsolution;
      SCIP_Real* costs;
      SCIP_Real* currentcost;
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &oldsolution, numcliques) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &costs, n) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &currentcost, (int) n) );
      for (unsigned int i = 0; i < n; ++i)
         costs[i] = - SCIPvarGetObj(consdata->vars[i]);  // -1 because we maximize

      // run subgradient method
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Running subgradient method to select cliques ...\n");
      constexpr size_t maxiter = 20000;
#ifdef SCIP_DEBUG
      constexpr size_t iterfreq = 1000;
#endif
      SCIP_Real step = 10.0;
      constexpr SCIP_Real stepfactor = 0.5;  // fast decreasing factor
      constexpr SCIP_Real slowstepfactor = 0.9999;  // slow decreasing factor
      constexpr SCIP_Real stepbound = 1e-5;  // lower bound for step size
      constexpr SCIP_Real stabilization = 0.6; // convexity parameter
      SCIP_Real oldobj = 0.0;
      size_t numincobj = 0;
      constexpr size_t maxnumincobj = 2;
      for (size_t iter = 0; iter < maxiter; ++iter)
      {
         // init current cost
         for (unsigned int j = 0; j < n; ++j)
            currentcost[j] = costs[j];

         // compute current cost: each member of the clique has coefficient 1 and is multiplied with the dual solution
         for (size_t c = 0; c < numcliques; ++c)
         {
            for (int j = begcliques[c]; j < begcliques[c + 1]; ++j)
            {
               unsigned int pos = cliquenodes[j];
               assert( pos < n );
               currentcost[pos] -= solution[c];
            }
         }

         // compute objective and dual violation
         SCIP_Real obj = 0.0;
#ifdef SCIP_DEBUG
         SCIP_Real dualviol = 0.0;
#endif
         for (unsigned int j = 0; j < n; ++j)
         {
            // the solution is 1 if the current cost is positive (otherwise 0)
            if ( currentcost[j] > 0.0 )
            {
               obj += costs[j];
#ifdef SCIP_DEBUG
               dualviol += currentcost[j];
#endif
            }
         }

         // compute slack (= subgradient) and delta to previous solution
#ifdef SCIP_DEBUG
         SCIP_Real diffsolution = 0.0;
#endif
         for (size_t c = 0; c < numcliques; ++c)
         {
            obj += solution[c]; // for constant part \lambda\T b, b is always 1

            SCIP_Real slack = 1.0;
            for (int j = begcliques[c]; j < begcliques[c + 1]; ++j)
            {
               unsigned int pos = cliquenodes[j];

               // the solution is 1 if the current cost is positive (otherwise 0)
               if ( currentcost[pos] > 0.0 )
               {
                  slack -= 1.0;        // coefficient is 1
                  obj -= solution[c];  // for \lambda\T A
               }
            }

            // store old solution
            SCIP_Real oldsub = oldsolution[c];
            oldsolution[c] = solution[c];

            // compute new solution as convex combination with old solution
            solution[c] = solution[c] - step * (stabilization * slack + (1.0 - stabilization) * oldsub);

            // compute difference
#ifdef SCIP_DEBUG
            diffsolution += fabs(oldsolution[c] - solution[c]);
#endif

            // project on nonnegative space
            if ( solution[c] < 0.0 )
               solution[c] = 0.0;
         }

         // count number of consecutive worsenings in objective
         if ( obj >= oldobj )
            ++numincobj;
         else
            numincobj = 0;  // reset for improvements
         oldobj = obj;

         // update stepsize
         if ( numincobj >= maxnumincobj )
            step *= stepfactor;  // use fast decrease
         else
            step *= slowstepfactor;  // use slow decrease

         // possibly output
#ifdef SCIP_DEBUG
         if ( iter % iterfreq == 0 )
         {
            SCIPinfoMessage(scip, nullptr, "%6lu: obj: %8.2f, step: %5.4f, delta subgr: %6.3f, dual viol.: %6.3f\n", iter, obj, step, diffsolution, dualviol);
         }
#endif

         // stop if stepsize is too small
         if ( step < stepbound )
         {
#ifdef SCIP_DEBUG
            SCIPinfoMessage(scip, nullptr, "%6lu: obj: %8.2f, step: %5.4f, delta subgr: %6.3f, dual viol.: %6.3f\n", iter, obj, step, diffsolution, dualviol);
#endif
            break;
         }
      }
      SCIPfreeBlockMemoryArray(scip, &currentcost, (int) n);
      SCIPfreeBlockMemoryArray(scip, &costs, (int) n);
      SCIPfreeBlockMemoryArray(scip, &oldsolution, numcliques);
   }

   // sort solution
   int* idx = nullptr;
   // no sorting is needed for small number of cliques
   if ( (int) numcliques > conshdlrdata->maxninitial )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idx, numcliques) );
      for (size_t j = 0; j < numcliques; ++j)
         idx[j] = (int) j;

      SCIPsortDownRealInt(solution, idx, (int) numcliques);
   }

   // generate corresponding clique cuts
   SCIP_VAR** tmpvars = consdata->tmpvars;
   assert( tmpvars != nullptr );
   int ncliques = MIN((int) numcliques, conshdlrdata->maxninitial);
   for (int i = 0; i < ncliques; ++i)
   {
      // get clique index
      int c = i;
      if ( idx != nullptr )
         c = idx[i];
      assert( 0 <= c && c < (int) numcliques );

      int csize = 0;
      for (int j = begcliques[c]; j < begcliques[c+1]; ++j)
      {
         unsigned int pos = cliquenodes[j];
         tmpvars[csize++] = consdata->vars[pos];
      }

      // create row and add to SCIP
      SCIP_ROW* row;
#ifndef NDEBUG
      char name[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "clique%u", conshdlrdata->nsepacliques);
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
#else
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "clique", -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
#endif
      ++conshdlrdata->nsepacliques;
      SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, row, csize, tmpvars, 1.0) );
      SCIP_CALL( SCIPaddRow(scip, row, FALSE, &cutoff) );
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintRow(scip, row, nullptr) );
#endif
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
      ++ngen;

      if ( cutoff )
         break;
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Added %u initial clique cuts.\n", ngen);

   SCIPfreeBlockMemoryArrayNull(scip, &idx, numcliques);
   SCIPfreeBlockMemoryArray(scip, &solution, numcliques);

   SCIPfreeBlockMemoryArray(scip, &begcliques, m + 1);
   SCIPfreeBlockMemoryArray(scip, &cliquenodes, cliquenodesmaxsize);

   return SCIP_OKAY;
}


/** add extensions of edges to cliques into pool */
static
SCIP_RETCODE addEdgeExtensionsToPool(
   SCIP*                 scip,               //!< SCIP pointer
   SCIP_CONSHDLRDATA*    conshdlrdata,       //!< constraint handler data
   SCIP_CONS*            cons,               //!< clique constraint
   unsigned int&         ngen                //!< number of generated cliques
   )
{
   assert( scip != nullptr );
   assert( conshdlrdata != nullptr );
   assert( cons != nullptr );

   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   assert( consdata != nullptr );

   ngen = 0;

   // nothing to do if there are no nodes
   unsigned int n = consdata->n;
   if ( n <= 1 )
      return SCIP_OKAY;

   // get graph
   const Graph* G = consdata->G;
   assert( G != nullptr );
   vectorpool* pool = consdata->pool;

   // prepare clique marker
   size_t edgenum = 0;
   std::vector<int> cand(n, -1);
   std::vector<unsigned int> C;
   C.reserve(n);

   // first compute edge covering by cliques - this might be too large
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Computing cliques covering all edges ...\n");
   SCIP_Bool extended = FALSE;
   EdgeIterator eit, eend;
   for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
   {
      Edge e = *eit;
      Vertex source = boost::source(e, *G);
      Vertex target = boost::target(e, *G);

      // skip edge if it has already been part of a clique (if cand is equal for both end nodes)
      if ( cand[source] >= 0 && cand[source] == cand[target] )
         continue;

      // skip edge if one of its variables is fixed
      SCIP_VAR* sourcevar = consdata->vars[source];
      SCIP_VAR* targetvar = consdata->vars[target];
      if ( SCIPvarGetLbLocal(sourcevar) > 0.5 || SCIPvarGetUbLocal(sourcevar) < 0.5 )
         continue;
      if ( SCIPvarGetLbLocal(targetvar) > 0.5 || SCIPvarGetUbLocal(targetvar) < 0.5 )
         continue;

      // initialize marker for clique
      cand[source] = (int) edgenum;
      cand[target] = (int) edgenum;

      C.clear();
      C.push_back(source);
      C.push_back(target);
      size_t csize = 2;

      // check whether we can extend the clique
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(source, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         // if already included in clique
         if ( cand[w] == (int) edgenum )
            continue;

         // check whether we can extend the clique
         unsigned int nneigh = 0;
         AdjacencyIterator ait2, aend2;
         for (boost::tie(ait2, aend2) = boost::adjacent_vertices(w, *G); ait2 != aend2; ++ait2)
         {
            Vertex x = *ait2;

            if ( cand[x] == (int) edgenum )
            {
               ++nneigh;
               if ( nneigh == csize )
                  break;
            }
         }

         // if the neighbors of w cover all nodes in the current clique, we can extend it
         if ( nneigh == csize )
         {
            C.push_back(w);
            ++csize;
            cand[w] = (int) edgenum;
            extended = TRUE;
         }
      }
      if ( csize > 2 )
         (void) pool->insert(C);
      ++ngen;
      ++edgenum;
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Added %u initial clique cuts to pool.\n", ngen);

   // if we could not extend the clique, we are triangle-free
   if ( ! extended )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Graph of constraint <%s> is triangle free.\n", SCIPconsGetName(cons));
      consdata->trianglefree = TRUE;
   }

   return SCIP_OKAY;
}


/** propagate marked variables
 *
 *  We loop through the colors and check the following: If a variable is fixed to 1, all variables corresponding to
 *  adjacent vertices with this color are fixed to 0.
 *
 *  @pre We do not have isolated nodes.
 */
static
SCIP_RETCODE propMarkedVariables(
   SCIP*                 scip,               //!< SCIP pointer
   SCIP_CONS*            cons,               //!< constraint to be propagated
   SCIP_Bool&            cutoff,             //!< whether a cutoff was detected
   unsigned int&         ngen                //!< number of propagated variables
   )
{
   assert( scip != nullptr );
   assert( cons != nullptr );

   ngen = 0;
   cutoff = FALSE;

   // get data of constraint
   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   assert( consdata != nullptr );
   assert( consdata->vars != nullptr );
   assert( consdata->varschanged != nullptr );
   assert( consdata->G != nullptr );
   assert( consdata->nvarschanged <= consdata->n );

   // exit if no variable has been changed
   if ( consdata->nvarschanged == 0 )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Propagating marked %u variables of constraint <%s>.\n", consdata->nvarschanged, SCIPconsGetName(cons));

   unsigned int n = consdata->n;
   SCIP_VAR** vars = consdata->vars;
   SCIP_VAR** varschanged = consdata->varschanged;
   const Graph* G = consdata->G;

   // make sure to propagate all variables if list is full
   if ( consdata->nvarschanged >= n )
   {
      for (unsigned int i = 0; i < n; ++i)
         varschanged[i] = vars[i];
      consdata->nvarschanged = n;
   }

   /* Loop through list: The following does not cause conflicts, because we only fix variables to 0 and only variables
    * that are fixed to 1 are entered into the list - thus the list stays untouched. */
   for (unsigned l = 0; l < consdata->nvarschanged; ++l)
   {
      // get changed variable
      SCIP_VAR* changedvar = varschanged[l];
      assert( changedvar != nullptr );

      // possibly get non-negated variables (variables might have been aggregated in between)
      if ( SCIPvarIsNegated(changedvar) )
         changedvar = SCIPvarGetNegatedVar(changedvar);
      assert( ! SCIPvarIsNegated(changedvar) );

      // first check whether the variable was globally fixed to 1 (can also happen deeper in the tree)
      if ( SCIPvarGetLbGlobal(changedvar) > 0.5 )
      {
         SCIPdebugMsg(scip, "Check globally changed variable <%s>.\n", SCIPvarGetName(changedvar));

         // loop through neighbors of v and fix them to 0
         Vertex v = (Vertex) SCIPvarGetData(changedvar);
         assert( v < consdata->n );
         OutEdgeIterator eit, eend;
         for (boost::tie(eit, eend) = boost::out_edges(v, *G); eit != eend; ++eit)
         {
            assert( v == boost::source(*eit, *G) );
            Vertex w = boost::target(*eit, *G);

            // if variable is alredy fixed to 1, we are infeasible
            if ( SCIPvarGetLbGlobal(vars[w]) > 0.5 )
            {
               cutoff = TRUE;
               return SCIP_OKAY;
            }

            // fix variable to 0
            if ( SCIPvarGetUbGlobal(vars[w]) > 0.5 )
            {
               SCIP_CALL( SCIPchgVarUbGlobal(scip, vars[w], 0.0) );
               SCIPdebugMsg(scip, "Globally fix variable <%s> to 0, because variable <%s> is fixed to 1.\n", SCIPvarGetName(vars[w]), SCIPvarGetName(changedvar));
               ++ngen;
            }
         }
      }
      // if the variable is fixed to 1 -> adjacent nodes get fixed to 0
      else if ( SCIPvarGetLbLocal(changedvar) > 0.5 )
      {
         SCIPdebugMsg(scip, "Check changed variable <%s>.\n", SCIPvarGetName(changedvar));

         // loop through neighbors of v
         Vertex v = (Vertex) SCIPvarGetData(changedvar);
         assert( v < consdata->n );
         OutEdgeIterator eit, eend;
         for (boost::tie(eit, eend) = boost::out_edges(v, *G); eit != eend; ++eit)
         {
            assert( v == boost::source(*eit, *G) );
            Vertex w = boost::target(*eit, *G);

            // fix variable to 0
            SCIP_Bool tightened = FALSE;
            SCIP_CALL( SCIPinferBinvarCons(scip, vars[w], FALSE, cons, (int) v, &cutoff, &tightened) );
            if ( cutoff )
            {
               SCIPdebugMsg(scip, " -> node infeasible (variable <%s> is 1 and adjacent variable <%s> fixed to 1, too).\n",
                  SCIPvarGetName(changedvar), SCIPvarGetName(vars[w]));

               // perform conflict analysis
               if ( SCIPisConflictAnalysisApplicable(scip) )
               {
                  SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[v]) );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[w]) );
                  SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, nullptr) );
               }

               /* Keep variable as changed, in order to find cutoffs at future nodes.  This is necessary, because it may
                * happen that SCIP does not generate bound change events, if the bound changes corresponding to an
                * infeasible node cancel out. */
               varschanged[0] = changedvar;
               assert( ! SCIPvarIsNegated(varschanged[0]) );
               consdata->nvarschanged = 1;
               return SCIP_OKAY;
            }
            if ( tightened )
            {
               SCIPdebugMsg(scip, "Fix variable <%s> to 0, because variable <%s> is fixed to 1.\n", SCIPvarGetName(vars[w]), SCIPvarGetName(changedvar));
               ++ngen;
            }
         }
      }
   }

   SCIPdebugMsg(scip, "Propagated marked variables for constraint <%s>: %u  (totally marked: %u)\n",
      SCIPconsGetName(cons), ngen, consdata->nvarschanged);
   consdata->nvarschanged = 0;

   return SCIP_OKAY;
}


//--------------------------------------------------------------------------------------------
//--------------------------------- SCIP callback functions ----------------------------------
//--------------------------------------------------------------------------------------------

//! free constraint handler data
static
SCIP_DECL_CONSFREE(BACSconshdlrFreeClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   // free constraint handler data
   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != nullptr );

   // output statistics
   if ( SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_MINIMAL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\nNumber of initial cliques:\t\t\t\t%8u\n", conshdlrdata->ninitialcliques);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Number of clique cuts separated:\t\t\t%8u\n", conshdlrdata->nsepacliques);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Number of clique cuts separated from the pool:\t\t%8u\n", conshdlrdata->nsepapoolcliques);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Number of clique cuts separated by exact algorithm:\t%8u\n", conshdlrdata->nsepaexactcliques);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Number of clique cuts separated by Greedy algorithm:\t%8u\n", conshdlrdata->nsepagreedycliques);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\nNumber of separation rounds with B&B code:\t%u\n", conshdlrdata->nsepabbcliquerounds);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Number of times the B&B code hit node limit:\t%u\n", conshdlrdata->nbbnodelimit);

      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\nTime for exact clique separation:\t\t%8.2f\n",
         SCIPgetClockTime(scip, conshdlrdata->timeexactsepa));
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Time for greedy clique separation:\t\t%8.2f\n",
         SCIPgetClockTime(scip, conshdlrdata->timegreedysepa));
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Time for pool clique separation:\t\t%8.2f\n",
         SCIPgetClockTime(scip, conshdlrdata->timepoolsepa));
   }

   if ( conshdlrdata->timeexactsepa != nullptr )
   {
      SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->timeexactsepa) );
      SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->timegreedysepa) );
      SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->timepoolsepa) );
   }

   // free constraint halder data
   SCIPfreeBlockMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, nullptr);

   return SCIP_OKAY;
}


//! frees constraint data
static
SCIP_DECL_CONSDELETE(BACSdeleteClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( consdata != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIPdebugMsg(scip, "Deleting clique constraint <%s>.\n", SCIPconsGetName(cons));

   unsigned int n = (*consdata)->n;
   vectorpool* pool = (*consdata)->pool;
   if ( pool != nullptr && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_MINIMAL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\nTotal number of cliques stored in pool:\t\t\t%8zu\n", pool->size());

      if ( pool->size() > 0 )
      {
         // get size statistics
         unsigned int* nSize = nullptr;
         SCIP_CALL( SCIPallocBufferArray(scip, &nSize, n) );
         for (unsigned int i = 0; i < n; ++i)
            nSize[i] = 0;
         unsigned int maxSize = 0;
         for (unsigned int i = 0; i < pool->size(); ++i)
         {
            unsigned int size = (unsigned int) (*pool)[i]->size();
            ++nSize[size];
            if ( size > maxSize )
               maxSize = size;
         }

         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Clique sizes:\n");
         for (unsigned int i = 2; i <= maxSize; ++i)
         {
            int w = (int) ceil(log10(MAX(i, nSize[i]) + 0.1));
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "%*u ", w, i);
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\n");
         for (unsigned int i = 2; i <= maxSize; ++i)
         {
            int w = (int) ceil(log10(MAX(i, nSize[i]) + 0.1));
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "%*u ", w, nSize[i]);
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\n");
         SCIPfreeBufferArray(scip, &nSize);
      }
   }
   if ( pool != nullptr )
      delete (*consdata)->pool;

   assert( (*consdata)->vars != nullptr );
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, n);

   if ( (*consdata)->tmpvars != nullptr )
   {
      assert( (*consdata)->vals != nullptr );
      assert( (*consdata)->cliquenodes != nullptr );
      assert( (*consdata)->varschanged != nullptr );

      SCIPfreeBlockMemoryArray(scip, &(*consdata)->tmpvars, n);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->vals, n);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->cliquenodes, n);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->varschanged, n);
   }

   // free TCG if necessary
   if ( (*consdata)->TCG != nullptr )
      tcliqueFree(&(*consdata)->TCG);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


//! transformation of constraints
static
SCIP_DECL_CONSTRANS(BACStransClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != nullptr );
   assert( targetcons != nullptr );

   SCIPdebugMsg(scip, "Transforming clique constraint <%s>.\n", SCIPconsGetName(sourcecons));

   // get data of original constraint
   SCIP_CONSDATA* sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != nullptr );
   assert( sourcedata->G != nullptr );

   // create transformed constraint data (copy data where necessary)
   unsigned int n = sourcedata->n;
   SCIP_CONSDATA* consdata = nullptr;
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->G = sourcedata->G;
   consdata->isolated = sourcedata->isolated;
   consdata->TCG = sourcedata->TCG;
   consdata->n = n;

   // mark all variables as changed in order to also handle fixed bounds from the initialization
   consdata->nvarschanged = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->varschanged, n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vals, n) );
   for (unsigned int i = 0; i < n; ++i)
   {
      assert( sourcedata->vars[i] != nullptr );
      SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->vars[i], &(consdata->vars[i])) );
      consdata->varschanged[consdata->nvarschanged++] = consdata->vars[i];
   }
   assert( consdata->nvarschanged == n );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->tmpvars, n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->cliquenodes, n) );

   consdata->scaleclique = sourcedata->scaleclique;
   consdata->cutoff = FALSE;
   consdata->scip = sourcedata->scip;
   consdata->pool = new vectorpool(1000, 100000);
   consdata->lastclique = 0;
   consdata->trianglefree = FALSE;

   // create transformed constraint
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   // -----------------------------------------------------------------------------------------------
   // catch events
   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != nullptr );
   assert( conshdlrdata->eventhdlr != nullptr );

   for (unsigned int i = 0; i < n; ++i)
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[i], SCIP_EVENTTYPE_LBTIGHTENED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*) consdata, nullptr) );
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[i], SCIP_EVENTTYPE_GLBCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*) consdata, nullptr) );
   }

   return SCIP_OKAY;
}


//! LP initialization method of constraint handler
static
SCIP_DECL_CONSINITLP(BACSinitlpClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( infeasible != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != nullptr );

   if ( conshdlrdata->initialcliques )
   {
      // loop through constraints
      for (int c = 0; c < nconss; ++c)
      {
         SCIP_CONS* cons = conss[c];
         assert( cons != nullptr );
         SCIPdebugMsg(scip, "Initialize LP method for clique constraint <%s>.\n", SCIPconsGetName(cons));

         unsigned int ngen = 0;
         if ( conshdlrdata->initialselection && SCIPgetSubscipDepth(scip) == 0 && SCIPgetNRuns(scip) == 1 )
         {
            SCIP_CALL( generateInitialCliqueCuts(scip, conshdlrdata, cons, ngen, *infeasible) );
         }
         else if ( conshdlrdata->initialpool && SCIPgetSubscipDepth(scip) == 0 && SCIPgetNRuns(scip) == 1 )
         {
            SCIP_CALL( addEdgeExtensionsToPool(scip, conshdlrdata, cons, ngen) );
         }
         else
         {
            // get data of constraint
            SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
            assert( consdata != nullptr );

            // setup initial solution (all 1, except for fixed vars)
            SCIP_CALL( setupInitialSolVals(scip, consdata) );

            // use enforcement method
            SCIP_CALL( enforceLPSol(scip, conshdlrdata, cons, TRUE, ngen, *infeasible) );
         }
         SCIPdebugMsg(scip, "Created %u initial clique cuts.\n", ngen);
         conshdlrdata->ninitialcliques += ngen;
      }
   }

   return SCIP_OKAY;
}


//! Separation method for clique constraints.
static
SCIP_DECL_CONSSEPALP(BACSsepalpClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != nullptr );

   *result = SCIP_DIDNOTRUN;
   unsigned int ngen = 0;
   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != nullptr );

   // loop through constraints
   for (int c = 0; c < nconss; ++c)
   {
      // get data of constraint
      assert( conss[c] != nullptr );
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]);
      assert( consdata != nullptr );
      assert( consdata->vals != nullptr );

      SCIPdebugMsg(scip, "Separating clique constraint <%s> ...\n", SCIPconsGetName(conss[c]));
      *result = SCIP_DIDNOTFIND;

      // get LP solution
      SCIP_CALL( getSolValues(scip, nullptr, consdata) );

      unsigned int ngenpool = 0;
      SCIP_Bool cutoff;
      SCIP_CALL( separatePoolCliques(scip, conshdlrdata, conss[c], MAX_CLIQUEPOOL_GEN, ngenpool, cutoff) );
      if ( cutoff )
      {
         SCIPdebugMsg(scip, "Clique separation from the pool detected a cutoff.\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      SCIPdebugMsg(scip, "Separated clique cuts from the pool: %u [%lu].\n", ngenpool, consdata->pool->size());
      conshdlrdata->nsepacliques += ngenpool;
      ngen += ngenpool;

      // separate all cliques if we have not found enough cuts or at the beginning of the root node
      int depth = SCIPgetDepth(scip);
      if ( ngenpool <= 10 || (depth == 0 && SCIPgetNSepaRounds(scip) <= 10) )
      {
         unsigned int ngennew = 0;
         SCIP_CALL( separateCliques(scip, conshdlrdata, conss[c], ngennew, cutoff) );
         if ( cutoff )
         {
            SCIPdebugMsg(scip, "Clique separation from the pool detected a cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         SCIPdebugMsg(scip, "Separated clique inequalities: %u\n", ngennew);
         ngen += ngennew;
      }
   }
   if ( ngen > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


//! solution separation method for clique constraints
static
SCIP_DECL_CONSSEPASOL(BACSsepasolClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != nullptr );

   *result = SCIP_DIDNOTRUN;
   unsigned int ngen = 0;
   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != nullptr );

   // loop through constraints
   for (int c = 0; c < nconss; ++c)
   {
      // get data of constraint
      assert( conss[c] != nullptr );
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]);
      assert( consdata != nullptr );
      assert( consdata->vals != nullptr );

      SCIPdebugMsg(scip, "Separating clique constraint <%s> ...\n", SCIPconsGetName(conss[c]));
      *result = SCIP_DIDNOTFIND;

      // get solution
      SCIP_CALL( getSolValues(scip, sol, consdata) );

      unsigned int ngenpool = 0;
      SCIP_Bool cutoff;
      SCIP_CALL( separatePoolCliques(scip, conshdlrdata, conss[c], MAX_CLIQUEPOOL_GEN, ngenpool, cutoff) );
      if ( cutoff )
      {
         SCIPdebugMsg(scip, "Clique separation from the pool detected a cutoff.\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      SCIPdebugMsg(scip, "Separated clique cuts from the pool: %u [%lu].\n", ngenpool, consdata->pool->size());
      conshdlrdata->nsepacliques += ngenpool;
      ngen += ngenpool;

      // separate all cliques if we have not found enough cuts
      if ( ngenpool <= 10 )
      {
         unsigned int ngennew = 0;
         SCIP_CALL( separateCliques(scip, conshdlrdata, conss[c], ngennew, cutoff) );
         if ( cutoff )
         {
            SCIPdebugMsg(scip, "Clique separation from the pool detected a cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         SCIPdebugMsg(scip, "Separated clique inequalities: %u\n", ngennew);
         ngen += ngennew;
      }
   }
   if ( ngen > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** Enforcing method for clique constraints: check feasibility of an integral solution.
 *
 *  We only enforce edge inequalities since this implies enforcing of
 *  the clique inequalities if the solution is integral.
 *
 *  @pre It is assumed that the solution is integral (this can be ensured by appropriate priorities).
 */
static
SCIP_DECL_CONSENFOLP(BACSenfolpClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != nullptr );

   *result = SCIP_FEASIBLE;
   unsigned int ngen = 0;
   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != nullptr );

   // we have a negative priority, so we should come after the integrality conshdlr.
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   // loop through constraints
   for (int c = 0; c < nconss; ++c)
   {
      SCIP_CONS* cons = conss[c];
      assert( cons != nullptr );

      // get data of constraint
      SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
      assert( consdata != nullptr );
      assert( consdata->vals != nullptr );

      SCIPdebugMsg(scip, "Enforcing clique constraint <%s>.\n", SCIPconsGetName(cons));

      // get LP solution
      SCIP_CALL( getSolValues(scip, nullptr, consdata) );

      if ( conshdlrdata->edgeenforce )
      {
         SCIP_Bool cutoff;
         SCIP_CALL( enforceLPSol(scip, conshdlrdata, cons, FALSE, ngen, cutoff) );
         if ( cutoff )
         {
            SCIPdebugMsg(scip, "Clique enforcing detected a cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         if ( ngen > 0 )
            *result = SCIP_SEPARATED;
      }
      else
      {
         SCIP_VAR** vars = consdata->vars;
         SCIP_Real* vals = consdata->vals;
         SCIP_VAR** tmpvars = consdata->tmpvars;
         const Graph* G = consdata->G;
         assert( vars != nullptr );
         assert( vals != nullptr );
         assert( tmpvars != nullptr );
         assert( G != nullptr );

         SCIPdebugMsg(scip, "Separating pool of clique constraint <%s> ...\n", SCIPconsGetName(cons));
         unsigned int ngenpool = 0;
         SCIP_Bool cutoff;
         SCIP_CALL( separatePoolCliques(scip, conshdlrdata, cons, MAX_CLIQUEPOOL_GEN, ngenpool,  cutoff) );
         if ( cutoff )
         {
            SCIPdebugMsg(scip, "Clique separation from the pool detected a cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         SCIPdebugMsg(scip, "Separated clique cuts from the pool: %u [%lu].\n", ngenpool, consdata->pool->size());
         conshdlrdata->nsepacliques += ngenpool;
         ngen += ngenpool;

         // only check edge inequalities if no cuts were separated from the pool
         if ( ngenpool > 0 )
            *result = SCIP_SEPARATED;
         else
         {
            // separate cliques by branch-and-bound if we found no cuts in the pool
            unsigned int ngennew = 0;
            SCIP_CALL( separateCliques(scip, conshdlrdata, cons, ngennew, cutoff) );
            if ( cutoff )
            {
               SCIPdebugMsg(scip, "Clique separation detected a cutoff.\n");
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }

            SCIPdebugMsg(scip, "Separated clique inequalities: %u\n", ngennew);
            ngen += ngennew;

            // loop through edges if no cliques where found (to ensure correct behavior)
            if ( ngennew == 0 )
            {
               // loop through all edges
               EdgeIterator eit, eend;
               for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
               {
                  Edge e = *eit;
                  Vertex source = boost::source(e, *G);
                  Vertex target = boost::target(e, *G);

                  // if the weight of the edge is larger than 1.0
                  if ( SCIPisEfficacious(scip, vals[source] + vals[target] - 1.0) )
                  {
                     tmpvars[0] = vars[source];
                     tmpvars[1] = vars[target];

                     SCIP_ROW* row;
#ifdef NDEBUG
                     char name[SCIP_MAXSTRLEN];
                     (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "edge%d#%d", source, target);
                     SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, name, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
#else
                     SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "edge", -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
#endif
                     SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, row, 2, tmpvars, 1.0) );
#ifdef SHOW_CLIQUES
                     SCIP_CALL( SCIPprintRow(scip, row, nullptr) );
#endif
                     SCIP_CALL( SCIPaddRow(scip, row, FALSE, &cutoff) );
                     if ( cutoff )
                     {
                        SCIP_CALL( SCIPreleaseRow(scip, &row) );
                        *result = SCIP_CUTOFF;
                        return SCIP_OKAY;
                     }

                     // cut should be violated:
                     assert( SCIPisFeasNegative(scip, SCIPgetRowLPFeasibility(scip, row)) );

                     SCIP_CALL( SCIPreleaseRow(scip, &row) );
                     SCIPdebugMsg(scip, "Found violated edge inequality.\n");
                     *result = SCIP_SEPARATED;
                     ++ngen;
                  }
               }
            }
            else
               *result = SCIP_SEPARATED;
         }
      }
      SCIPdebugMsg(scip, "Enforced %u cuts from clique constraint <%s>.\n", ngen, SCIPconsGetName(conss[c]));
   }
   return SCIP_OKAY;
}


//! enforcing method for pseudo solutions
static
SCIP_DECL_CONSENFOPS(BACSenfopsClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != nullptr );

   *result = SCIP_DIDNOTRUN;
   if ( objinfeasible )
      return SCIP_OKAY;

   *result = SCIP_FEASIBLE;
   if ( solinfeasible )
      return SCIP_OKAY;

   // loop through constraints
   for (int c = 0; c < nconss; ++c)
   {
      // get data of constraint
      assert( conss[c] != nullptr );
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]);
      assert( consdata != nullptr );
      assert( consdata->vars != nullptr );
      assert( consdata->G != nullptr );

      SCIP_VAR** vars = consdata->vars;
      const Graph* G = consdata->G;

      // loop through all edges
      EdgeIterator eit, eend;
      for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
      {
         Edge e = *eit;
         Vertex source = boost::source(e, *G);
         Vertex target = boost::target(e, *G);

         SCIP_Real valss = SCIPgetSolVal(scip, nullptr, vars[source]);
         SCIP_Real valst = SCIPgetSolVal(scip, nullptr, vars[target]);

         assert( SCIPisIntegral(scip, valss) );
         assert( SCIPisIntegral(scip, valst) );

         if ( SCIPisGT(scip, valss + valst, 1.0) )
         {
            SCIPdebugMsg(scip, "Solution is infeasible.\n");
            *result = SCIP_INFEASIBLE;
            return SCIP_OKAY;
         }
      }
   }
   SCIPdebugMsg(scip, "Solution is feasible.\n");
   return SCIP_OKAY;
}


//! check feasibility of an integral soltuion
static
SCIP_DECL_CONSCHECK(BACScheckClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != nullptr );

   SCIPdebugMsg(scip, "Check method for clique constraint handler. Checking edge inequalities.\n");
   *result = SCIP_FEASIBLE;

   // loop through constraints
   for (int c = 0; c < nconss; ++c)
   {
      // get data of constraint
      assert( conss[c] != nullptr );
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]);
      assert( consdata != nullptr );
      assert( consdata->vars != nullptr );
      assert( consdata->G != nullptr );

      SCIP_VAR** vars = consdata->vars;
      const Graph* G = consdata->G;

      // possibly temporarily allocate space
      SCIP_Bool allocvals = FALSE;
      if ( consdata->vals == nullptr )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &consdata->vals, consdata->n) );
         allocvals = TRUE;
      }
      // get solution values
      SCIP_CALL( getSolValues(scip, sol, consdata) );
      SCIP_Real* values = consdata->vals;
      assert( values != nullptr );

      // loop through all edges
      EdgeIterator eit, eend;
      for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
      {
         Edge e = *eit;
         Vertex source = boost::source(e, *G);
         Vertex target = boost::target(e, *G);

         SCIP_Real vals = values[source];
         SCIP_Real valt = values[target];

         assert( SCIPisFeasIntegral(scip, vals) );
         assert( SCIPisFeasIntegral(scip, valt) );

         if ( SCIPisFeasGT(scip, vals + valt, 1.0) )
         {
            if ( printreason )
            {
               SCIPinfoMessage(scip, nullptr, "Edge constraint %s + %s <= 1.0 violated (%f + %f > 1.0).\n",
                  SCIPvarGetName(vars[source]), SCIPvarGetName(vars[target]), vals, valt);
            }
            SCIPdebugMsg(scip, "Solution is infeasible.\n");
            *result = SCIP_INFEASIBLE;
            if ( allocvals )
               SCIPfreeBufferArray(scip, &consdata->vals);
            return SCIP_OKAY;
         }
      }
      if ( allocvals )
         SCIPfreeBufferArray(scip, &consdata->vals);
   }
   SCIPdebugMsg(scip, "Solution is feasible.\n");
   return SCIP_OKAY;
}


//! Domain propagation
static
SCIP_DECL_CONSPROP(BACSpropClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != nullptr );

   *result = SCIP_DIDNOTFIND;

   // loop through constraints
   for (int c = 0; c < nconss; ++c)
   {
      SCIP_Bool cutoff = FALSE;
      unsigned int ngen = 0;

      assert( conss[c] != nullptr );
      SCIP_CALL( propMarkedVariables(scip, conss[c], cutoff, ngen) );

      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if ( ngen > 0 )
         *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}


/** presolving method of constraint handler
 *
 *  We perform a simple propagation round
 */
static
SCIP_DECL_CONSPRESOL(BACSpresolClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != nullptr );

   SCIPdebugMsg(scip, "Presolving method of clique constraint handler. Propagating edge inequalities.\n");
   *result = SCIP_DIDNOTRUN;

   // loop through constraints
   for (int c = 0; c < nconss; ++c)
   {
      SCIP_Bool cutoff = FALSE;
      unsigned int ngen = 0;

      *result = SCIP_DIDNOTFIND;

      assert( conss[c] != nullptr );
      SCIP_CALL( propMarkedVariables(scip, conss[c], cutoff, ngen) );

      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if ( ngen > 0 )
      {
         *nfixedvars += (int) ngen;
         *result = SCIP_SUCCESS;
      }
   }

   return SCIP_OKAY;
}


/** Propagation resolution for conflict analysis
 *
 *  The inferinfo integer is set by scip_prop() as follows:
 *  If a variable is set to 0, because of a variable for node v set to 1, inferinfo is @a v.
 */
SCIP_DECL_CONSRESPROP(BACSrespropClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != nullptr );
   assert( infervar != nullptr );
   assert( bdchgidx != nullptr );
   assert( result != nullptr );

   SCIPdebugMsg(scip, "Propagation resolution method of clique constraint handler.\n");
   *result = SCIP_DIDNOTFIND;

   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   assert( consdata != nullptr );
   assert( consdata->n > 0 );
   assert( consdata->vars != nullptr );
   assert( consdata->G != nullptr );

   assert( 0 <= inferinfo && inferinfo < (int) consdata->n );

   SCIP_VAR** vars = consdata->vars;

   // if the variable was fixed to 0
   if ( SCIPvarGetUbAtIndex(infervar, bdchgidx, FALSE) > 0.5 && SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) < 0.5 )
   {
      SCIPdebugMsg(scip, " -> reason for setting <%s> = 0 was <%s> = 1.\n", SCIPvarGetName(infervar), SCIPvarGetName(vars[inferinfo]));
      SCIP_CALL( SCIPaddConflictLb(scip, vars[inferinfo], bdchgidx) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


//! Lock variables
static
SCIP_DECL_CONSLOCK(BACSlockClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != nullptr );

   SCIPdebugMsg(scip, "Locking method for clique constraint <%s>.\n", SCIPconsGetName(cons));

   // get data of original constraint
   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   assert( consdata != nullptr );
   assert( consdata->vars != nullptr );

   // we lock all variables corresponding to non-isolated nodes
   unsigned int n = consdata->n;
   SCIP_VAR** vars = consdata->vars;

   if ( consdata->isolated != nullptr )
   {
      const std::vector<bool> isolated = *consdata->isolated;

      for (unsigned int i = 0; i < n; ++i)
      {
         if ( ! isolated[i] )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, vars[i], nlocksneg, nlockspos) );
         }
      }
   }
   else
   {
      for (unsigned int i = 0; i < n; ++i)
      {
         SCIP_CALL( SCIPaddVarLocks(scip, vars[i], nlocksneg, nlockspos) );
      }
   }

   return SCIP_OKAY;
}


//! printing method of constraint handler
static
SCIP_DECL_CONSPRINT(BACSprintClique)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != nullptr );

   SCIPdebugMsg(scip, "Printing method for clique constraint <%s>.\n", SCIPconsGetName(cons));

   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   assert( consdata != nullptr );
   assert( consdata->G != nullptr );
   assert( consdata->vars != nullptr );

   SCIP_VAR** vars = consdata->vars;
   const Graph* G = consdata->G;

   SCIPinfoMessage(scip, file, "clique(");

   // loop through all edges
   unsigned int cnt = 0;
   EdgeIterator eit, eend;
   for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
   {
      Edge e = *eit;
      Vertex source = boost::source(e, *G);
      Vertex target = boost::target(e, *G);

      if ( cnt > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "<%s> -- <%s>", SCIPvarGetName(vars[source]), SCIPvarGetName(vars[target]));
      ++cnt;
   }
   SCIPinfoMessage(scip, file, ")");

   return SCIP_OKAY;
}


//! copying method of constraint handler
static
SCIP_DECL_CONSHDLRCOPY(BACSconshdlrCopyCliuqe)
{
   assert( scip != nullptr );
   assert( conshdlr != nullptr );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( valid != nullptr );

   // call inclusion method of constraint handler
   SCIP_CALL( BACSincludeConshdlrClique(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}


//! copying method of constraint
static
SCIP_DECL_CONSCOPY(BACSconsCopyClique)
{
   assert( scip != nullptr );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) == 0 );
   assert( cons != nullptr );
   assert( sourcescip != nullptr );
   assert( sourcecons != nullptr );
   assert( varmap != nullptr );
   assert( valid != nullptr );

   *valid = TRUE;

   SCIPdebugMsg(scip, "Copying method for clique constraint <%s>.\n", SCIPconsGetName(sourcecons));

   // get data of original constraint
   SCIP_CONSDATA* sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != nullptr );
   assert( sourcedata->G != nullptr );
   assert( sourcedata->n > 0 );
   assert( sourcedata->vars != nullptr );

   // copy constraint data
   unsigned int n = sourcedata->n;
   SCIP_VAR** sourcevars = sourcedata->vars;

   // create new variables
   SCIP_VAR** vars = nullptr;

   // separately allocate space to account for unsuccessful copying
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, n) );

   // get copies
   for (unsigned int i = 0; i < n; ++i)
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[i], &vars[i], varmap, consmap, global, valid) );
      assert( *valid );
      assert( vars[i] != nullptr );

      // make sure for older SCIP versions that the data is correct
      SCIPvarSetData(vars[i], (SCIP_VARDATA*) (size_t) i);
   }

   // create copied constraint
   if ( name == nullptr )
      name = SCIPconsGetName(sourcecons);

   SCIP_CALL( BACScreateConsClique(scip, cons, name,
         sourcedata->G, sourcedata->isolated, n, vars,
         initial, separate, enforce, check, propagate) );

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

#if SCIP_VERSION >= 900
//! construct symmetry detection graph
static
SCIP_DECL_CONSGETPERMSYMGRAPH(BACSconsGetPermsymClique)
{
   assert( scip != nullptr );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0 );
   assert( cons != nullptr );
   assert( success != nullptr );

   SCIPdebugMsg(scip, "Creating symmetry detection graph for clique constraint <%s>.\n", SCIPconsGetName(cons));

   // get data
   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   assert( consdata != nullptr );
   assert( consdata->G != nullptr );

   unsigned int n = consdata->n;
   if ( n > 0 )
   {
      const Graph& G = *consdata->G;
      SCIP_VAR** vars = consdata->vars;

      // add edges as set packing constraints
      SCIP_VAR** activevars;
      SCIP_Real* activevals;

      SCIP_CALL( SCIPallocBufferArray(scip, &activevars, 2) );
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, 2) );

      // loop through graph
      EdgeIterator eit, eend;
      for (boost::tie(eit, eend) = boost::edges(G); eit != eend; ++eit)
      {
         Edge e = *eit;
         Vertex s = boost::source(e, G);
         Vertex t = boost::target(e, G);

         if ( ( SCIPvarGetLbLocal(vars[s]) < 0.5 && SCIPvarGetUbLocal(vars[s]) > 0.5 )
            && ( SCIPvarGetLbLocal(vars[t]) < 0.5 && SCIPvarGetUbLocal(vars[t]) > 0.5 ) )
         {
            activevars[0] = vars[s];
            activevars[1] = vars[t];
            activevals[0] = 1.0;
            activevals[1] = 1.0;

            int nlocvars = 2;
            SCIP_Real constant = 0.0;
            SCIP_CALL( SCIPgetSymActiveVariables(scip, SYM_SYMTYPE_PERM, &activevars, &activevals, &nlocvars, &constant, TRUE) );

            SCIP_Real lhs = - SCIPinfinity(scip);
            SCIP_Real rhs = 1.0 - constant;

            SCIP_CALL( SCIPextendPermsymDetectionGraphLinear(scip, graph, activevars, activevals, nlocvars, cons, lhs, rhs, success) );
            assert( *success );
         }
      }

      SCIPfreeBufferArray(scip, &activevals);
      SCIPfreeBufferArray(scip, &activevars);
   }

   *success = TRUE;

   return SCIP_OKAY;
}
#endif


//--------------------------------------------------------------------------------------------
//------------------------------ constraint (handler) generation -----------------------------
//--------------------------------------------------------------------------------------------

//! creates the handler for clique constraints and includes it in SCIP
SCIP_RETCODE BACSincludeConshdlrClique(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = nullptr;
   SCIP_CONSHDLR* conshdlr;

   // create constraint handler data
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   // initialize data
   conshdlrdata->nsepacliques = 0;
   conshdlrdata->ninitialcliques = 0;
   conshdlrdata->nsepapoolcliques = 0;
   conshdlrdata->nsepaexactcliques = 0;
   conshdlrdata->nsepagreedycliques = 0;
   conshdlrdata->nsepabbcliquerounds = 0;
   conshdlrdata->nbbnodelimit = 0;
   conshdlrdata->eventhdlr = nullptr;

   SCIP_CALL( SCIPcreateCPUClock(scip, &conshdlrdata->timeexactsepa) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &conshdlrdata->timegreedysepa) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &conshdlrdata->timepoolsepa) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         BACSenfolpClique, BACSenfopsClique, BACScheckClique, BACSlockClique,
         conshdlrdata) );
   assert( conshdlr != nullptr );

   //  set non-fundamental callbacks via specific setter functions
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, BACSconshdlrCopyCliuqe, BACSconsCopyClique) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, BACSconshdlrFreeClique) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, BACSdeleteClique) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, BACSpresolClique, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, BACSprintClique) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, BACSpropClique, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROPTIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, BACSrespropClique) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, BACSinitlpClique) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, BACSsepalpClique, BACSsepasolClique, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, BACStransClique) );
#if SCIP_VERSION >= 900
   SCIP_CALL( SCIPsetConshdlrGetPermsymGraph(scip, conshdlr, BACSconsGetPermsymClique) );
#endif

   // include event handler
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, BACSeventExecClique, (SCIP_EVENTHDLRDATA*) conshdlrdata) );

   // get event handler
   conshdlrdata->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if ( conshdlrdata->eventhdlr == nullptr )
   {
      SCIPerrorMessage("event handler for clique constraints not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/clique/maxnsepa",
         "maximal number of cuts separated in clique constraints",
         &conshdlrdata->maxnsepa, TRUE, DEFAULT_MAXNSEPA, 0, INT_MAX, nullptr, nullptr) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/clique/maxninitial",
         "maximal number of clique inequalities added initially",
         &conshdlrdata->maxninitial, TRUE, DEFAULT_MAXNINITIAL, 0, INT_MAX, nullptr, nullptr) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/clique/maxcliquebbnodes",
         "maximal number of nodes in clique B&B code",
         &conshdlrdata->maxcliquebbnodes, TRUE, DEFAULT_MAXCLIQUEBBNODES, 0, INT_MAX, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/clique/initialcliques",
         "Should initial cliques be generated?",
         &(conshdlrdata)->initialcliques, TRUE, DEFAULT_INITIALCLIQUES, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/clique/initialselection",
         "Should initial cliques be selected from an edge covering?",
         &(conshdlrdata)->initialselection, TRUE, DEFAULT_INITIALSELECTION, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/clique/initialpool",
         "Should initial cliques be added to the pool?",
         &(conshdlrdata)->initialpool, TRUE, DEFAULT_INITIALPOOL, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/clique/exactcliques",
         "Should exact clique generation be called if no other cuts have been found?",
         &(conshdlrdata)->exactcliques, TRUE, DEFAULT_EXACTCLIQUES, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/clique/greedycliques",
         "Should greedy clique generation be called efore exact code?",
         &(conshdlrdata)->greedycliques, TRUE, DEFAULT_GREEDYCLIQUES, nullptr, nullptr) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/clique/edgeenforce",
         "Run enforcement based on edges with clique extension?",
         &conshdlrdata->edgeenforce, TRUE, DEFAULT_EDGEENFORCE, nullptr, nullptr) );

   return SCIP_OKAY;
}


//! creates a clique constraint
SCIP_RETCODE BACScreateConsClique(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS**           cons,               /**< pointer to constraint (output) */
   const char*           name,               /**< name of constraint  */
   const Graph*          G,                  /**< underlying graph */
   const std::vector<bool>* isolated,        /**< marks whether a node is isolated */
   unsigned int          n,                  /**< number of nodes */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_Bool             initial,            /**< whether constraint is initial */
   SCIP_Bool             separate,           /**< whether constraint should be separated */
   SCIP_Bool             enforce,            /**< whether constraint should be enforced */
   SCIP_Bool             check,              /**< whether constraint should be checked */
   SCIP_Bool             propagate           /**< whether constraint should be propagated */
   )
{
   assert( scip != nullptr );
   assert( cons != nullptr );
   assert( G != nullptr );
   assert( vars != nullptr );
   assert( isolated == nullptr || isolated->size() == n );

   // find constraint handler
   SCIP_CONSHDLR* conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == nullptr )
   {
      SCIPerrorMessage("clique constraint handler not found\n");
      return SCIP_INVALIDCALL;
   }

   // init constraint data
   SCIP_CONSDATA* consdata = nullptr;
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   consdata->G = G;
   consdata->isolated = isolated;
   consdata->TCG = nullptr;
   consdata->n = n;

   // copy variables
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->vars, vars, n) );

   // only needed in transformed constraint:
   consdata->vals = nullptr;
   consdata->varschanged = nullptr;
   consdata->tmpvars = nullptr;

   // set default values:
   consdata->nvarschanged = 0;
   consdata->cliquenodes = nullptr;
   consdata->scaleclique = 1000.0;
   consdata->cutoff = FALSE;
   consdata->scip = scip;
   consdata->pool = nullptr;
   consdata->lastclique = 0;
   consdata->trianglefree = FALSE;

   // generate constraint
   // Bools: initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata,
         initial, separate, enforce, check, propagate, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
