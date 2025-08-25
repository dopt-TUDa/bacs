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

/**@file   heur_tabu.cpp
 * @ingroup DEFPLUGINS_HEUR
 * @brief  tabu search primal heuristic
 * @author Marc Pfetsch
 * @author Jonas Alker
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <set>

#include "heur_tabu.h"
#include "graph.h"
#include "struct_probdata.h"
#include "scip/scip_sol.h"
#include "scip/sol.h"

#define HEUR_NAME            "tabu"
#define HEUR_DESC            "tabu search for stable set solutions"
#define HEUR_DISPCHAR        'T'
#define HEUR_PRIORITY        10000
#define HEUR_FREQ            1
#define HEUR_FREQOFS         0
#define HEUR_MAXDEPTH        -1
#define HEUR_TIMING          ( SCIP_HEURTIMING_BEFOREPRESOL | SCIP_HEURTIMING_BEFORENODE )
#define HEUR_USESSUBSCIP     FALSE /**< does the heuristic use a secondary SCIP instance? */

#define TABU_L                50        //!< tabu interval
#define TABU_maxIter          30000     //!< maximal number of iterations
#define TABU_fixed            100000    //!< fixed nodes get huge iter number, assert( TABU_fixed > TABU_maxIter + TABU_L )
#define INIT_INDEX_SIZE       1000
#ifdef SCIP_DEBUG
#define TABU_freq             100       //!< frequency for output
#endif
#define DEFAULT_GLOBALBOUNDS TRUE       //!< use global variable bounds for initialization
#define DEFAULT_RANDOMTABU   TRUE       //!< use random tabu length
#define DEFAULT_TABUOUTNODE  TRUE       //!< apply tabu additionally on solution leaving node

/** primal heuristic data */
struct SCIP_HeurData
{
   int*                  knownindex;          //!< array of known indices
   int                   nknownsols;          //!< number of stored solutions
   int                   maxindex;            //!< maximum index seen by heuristic
   int                   knownindexsize;      //!< size of knownindex array, used for freeing data
   SCIP_RANDNUMGEN*      randnumgen;          //!< randon number generator
   SCIP_Bool             globalbounds;        //!< use global variable bounds for initialization
   SCIP_Bool             randomtabu;          //!< use random tabu length
   SCIP_Bool             tabuoutnode;         //!< apply tabu additionally on solution leaving node
};


#ifdef SCIP_DEBUG
//! compute objective value of a given solution
static
SCIP_Real computeSolObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   SCIP_SOL*             sol                 /**< sol */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( sol != nullptr );

   SCIP_Real* vals;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vals, probdata->n) );

   SCIPgetSolVals(scip, sol, probdata->n, probdata->vars, vals);

   SCIP_Real obj = 0.0;
   for (size_t i = 0; i < probdata->n; ++i)
   {
      obj += SCIPvarGetObj(probdata->vars[i]) * vals[i];
   }

   SCIPfreeBlockMemoryArray(scip, &vals, probdata->n);
   return obj;
}
#endif

//! check wether index of a given solution is already known
static
SCIP_Bool isIndexKnown(
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   index               /**< index of solution */
   )
{
   assert( heurdata != nullptr );

   // loop through knownsols array and check wether solution is known to heuristic
   // TODO: make knownindex array sorted, to improve running time
   for (int i = 0; i < heurdata->nknownsols; ++i)
   {
      if ( index == heurdata->knownindex[i] )
         return TRUE;
   }
   return FALSE;
}


/** calculation of tabu */
static
size_t calcTabu(
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   size_t                size,               /**< size of current S */
   int                   degree              /**< degree of node */
   )
{
   assert( heurdata != nullptr );
   assert( heurdata->randnumgen != nullptr );
   assert( degree >= 0 );

   return (size_t) MAX( (int) size, degree * SCIPrandomGetInt(heurdata->randnumgen, 0, 10) );
}


//! run tabu seach algorithm on given initial solution
static
SCIP_RETCODE tabuSearch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   SCIP_HEUR*            heur,               /**< current heuristic */
   SCIP_SOL*             init_sol,           /**< initial solution to start tabu search with */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( heur != nullptr );
   assert( init_sol != nullptr );
   assert( result != nullptr );

#ifdef SCIP_DEBUG
   SCIP_Real sol_obj = computeSolObj(scip, probdata, init_sol);
#endif

   size_t n = probdata->n;
   const Graph* G = probdata->G;
   assert( G != nullptr );

   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);

   // init matrices
   std::size_t* T;   // stores iteration at which tabu node will expire to be tabu
   SCIP_Real* A;     // weights of adjacent nodes in the solution
   SCIP_Bool* S;     // stores solution (initialized to 0 sol)
   SCIP_Bool* B;     // stores best solution
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &T, n) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &A, n) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &S, n) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &B, n) );

   // perform predefined number of iterations
   SCIP_Real bestobj = 0.0;
   SCIP_Real obj = 0.0;
   size_t solsize = 0;

   // add local fixings
   // - fixed to zero: set the tabu value for this vertex to TABU_fixed,
   //                  the vertex will never be added to a solution
   // - fixed to one: add vertex to current and best solution and fix all
   //                 neighbors to zero
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;

      SCIP_Real lb;
      SCIP_Real ub;
      if ( heurdata->globalbounds )
      {
         lb = SCIPvarGetLbGlobal(probdata->vars[v]) ;
         ub = SCIPvarGetUbGlobal(probdata->vars[v]) ;
      }
      else
      {
         lb = SCIPvarGetLbLocal(probdata->vars[v]) ;
         ub = SCIPvarGetUbLocal(probdata->vars[v]) ;
      }

      if ( ub < 0.5 )
         T[v] = TABU_fixed;
      else if ( lb > 0.5 )
      {
         S[v] = TRUE;
         B[v] = TRUE;

         SCIP_Real weight = boost::get(vertex_weight_t(), *G, v);
         bestobj += weight;
         obj += weight;

         ++solsize;
         assert( solsize < n );

         assert( SCIPisEQ(scip, A[v], 0.0) );
         AdjacencyIterator ait, aend;
         for (boost::tie(ait,aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            Vertex w = *ait;
            A[w] += weight;
            T[w] = TABU_fixed;
         }
      }
   }

   // add given solution
   // we want to start the tabu search with the currently best known solution sol, so we try to add all vertices included in init_sol to S.
   // Let v be a vertex in sol
   // - v is locally fixed to one:  v is alredy included in S by the previous step (add local fixings)
   // - v is locally fixed to zero: v cannot be added to S, we store the difference between the value of init_sol and S in delta_sol_obj
   // - else:                        add v to S
#ifdef SCIP_DEBUG
   SCIP_Real delta_sol_obj = 0.0;
#endif
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      if ( SCIPgetSolVal(scip, init_sol, probdata->vars[v]) > 0.5 )
      {
         if ( S[v] )                // v is locally fixed to 1.0
            continue;

         SCIP_Real weight = boost::get(vertex_weight_t(), *G, v);
         if ( T[v] > TABU_maxIter ) // v is locally fixed to 0.0
         {
#ifdef SCIP_DEBUG
            delta_sol_obj -= weight;
#endif
            continue;
         }
         assert( ! S[v] );
         assert( T[v] == 0);

         S[v] = TRUE;
         B[v] = TRUE;
         A[v] = 0.0;
         bestobj += weight;
         obj += weight;

         ++solsize;
         assert( solsize <= n );

         AdjacencyIterator ait, aend;
         for (boost::tie(ait,aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            Vertex w = *ait;
            A[w] += weight;
         }
      }
   }
#ifdef SCIP_DEBUG
   assert( SCIPisFeasEQ(scip, bestobj, sol_obj + delta_sol_obj) );
   assert( SCIPisFeasEQ(scip, obj, sol_obj + delta_sol_obj) );
#endif

   // main tabu loop with fixed number of iterations
   for (size_t iter = 1; iter <= TABU_maxIter; ++iter)
   {
      SCIP_Real bestdelta = - SCIP_REAL_MAX;
      SCIP_Real bestweight = SCIP_INVALID;
      Vertex bestnode = 0;

      // try all non-tabu nodes that are not in the solution and add the best
#ifdef SCIP_DEBUG
      SCIP_Bool aspiration = FALSE;
#endif
      for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
      {
         Vertex v = *vit;

         // skip nodes in the solution
         if ( S[v] )
            continue;

         // compute cost for adding v
         SCIP_Real weight = boost::get(vertex_weight_t(), *G, v);
         SCIP_Real objdelta = weight - A[v];

         // aspiration - directly take node if it improves objective and v is not locally fixed to zero
         if ( T[v] < TABU_fixed && objdelta > 0.0 )
         {
            bestdelta = objdelta;
            bestweight = weight;
            bestnode = v;
#ifdef SCIP_DEBUG
            aspiration = TRUE;
#endif
            break;
         }

         // take best non-tabu node
         if ( T[v] <= iter && objdelta > bestdelta )
         {
#ifdef NDEBUG
            if ( heurdata->globalbounds )
            {
               assert( SCIPvarGetUbGlobal(probdata->vars[v]) > 0.5 );
               assert( SCIPvarGetLbGlobal(probdata->vars[v]) < 0.5 );
            }
            else
            {
               assert( SCIPvarGetUbLocal(probdata->vars[v]) > 0.5 );
               assert( SCIPvarGetLbLocal(probdata->vars[v]) < 0.5 );
            }
#endif
            bestdelta = objdelta;
            bestweight = weight;
            bestnode = v;
         }
      }
      // if no node was selected we break the tabu loop, this happens on instances where #(unfixed variables) < TABU_L
      if ( bestweight == SCIP_INVALID )  /*lint !e777*/
         break;

      assert( ! S[bestnode] );
      assert( SCIPisEQ(scip, bestdelta, bestweight - A[bestnode]) );

      // put best node into solution
      S[bestnode] = TRUE;
      A[bestnode] = 0.0;

      ++solsize;
      assert( solsize <= n );

      // update tabu list
      if ( heurdata->randomtabu )
      {
         if ( heurdata->globalbounds )
            T[bestnode] = iter + calcTabu(heurdata, solsize, probdata->degrees[bestnode]);
         else
            T[bestnode] = iter + calcTabu(heurdata, solsize, probdata->localdegrees[bestnode]);
      }
      else
         T[bestnode] = iter + TABU_L;

      // remove adjacent nodes from solution and update data
      SCIP_Real newobj = obj + bestweight;
      AdjacencyIterator ait, aend;
      for (boost::tie(ait,aend) = boost::adjacent_vertices(bestnode, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;

         // if neighbor is in solution then remove it
         if ( S[w] )
         {
            SCIP_Real weight = boost::get(vertex_weight_t(), *G, w);
            newobj -= weight;
            S[w] = FALSE;

            assert( solsize > 0 );
            --solsize;

            if ( heurdata->tabuoutnode )
            {
               if ( heurdata->randomtabu )
               {
                  if ( heurdata->globalbounds )
                     T[w] = iter + calcTabu(heurdata, solsize, probdata->degrees[w]);
                  else
                     T[w] = iter + calcTabu(heurdata, solsize, probdata->localdegrees[w]);
               }
               else
                  T[w] = iter + TABU_L;
            }

            // correct weights of neighbors of w
            AdjacencyIterator ait2, aend2;
            for (boost::tie(ait2,aend2) = boost::adjacent_vertices(w, *G); ait2 != aend2; ++ait2)
            {
               Vertex x = *ait2;
               if ( ! S[x] )
                  A[x] -= weight;
            }
         }
         assert( ! S[w] );

         // bestnode entered solution, so correct weight for neighbor w
         A[w] += bestweight;
      }
      // update objective
      obj += bestdelta;
      assert( SCIPisEQ(scip, obj, newobj) );

      // store best solution
      if ( obj > bestobj )
      {
         bestobj = obj;
         BMScopyMemoryArray(B, S, n);
      }

#ifdef SCIP_DEBUG
      // possibly output
      if ( ((iter-1) % TABU_freq == 0 || obj > bestobj) )
      {
         if ( probdata->unweighted )
            SCIPinfoMessage(scip, nullptr, "Iter: %7lu\tobj: %4d\tbest: %4d\tnode: %4lu\tdelta: %d\t%c\n",
               iter, (int) obj, (int) bestobj, bestnode, (int) bestdelta, aspiration ? '*' : ' ');
         else
            SCIPinfoMessage(scip, nullptr, "Iter: %7lu\tobj: %4f\tbest: %4f\tnode: %4lu\tdelta: %f\t%c\n",
               iter, obj, bestobj, bestnode, bestdelta, aspiration ? '*' : ' ');
      }
#endif
   }

   // pass solution to SCIP
   SCIP_SOL* sol;
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
   for (size_t i = 0; i < n; ++i)
   {
      if ( B[i] )
         SCIP_CALL( SCIPsetSolVal(scip, sol, probdata->vars[i], 1.0) );
   }
   assert( SCIPisSumEQ(scip, bestobj + SCIPgetOrigObjoffset(scip), SCIPgetSolOrigObj(scip, sol)) );

   SCIP_Bool stored;
   SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, TRUE, FALSE, &stored) );

   if ( stored )
   {
      *result = SCIP_FOUNDSOL;
      SCIPdebugMsg(scip, "Tabu search found solution of value %f.\n", bestobj);
   }
   else
      SCIPdebugMsg(scip, "Tabu search did not find new solution.\n");

   SCIPfreeBlockMemoryArray(scip, &B, n);
   SCIPfreeBlockMemoryArray(scip, &S, n);
   SCIPfreeBlockMemoryArray(scip, &A, n);
   SCIPfreeBlockMemoryArray(scip, &T, n);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTabu)
{
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( result != nullptr );
   static_assert(TABU_fixed > TABU_maxIter + TABU_L, "TABU_fixed is not large enough.");

   *result = SCIP_DIDNOTRUN;

   SCIP_HeurData* heurdata = SCIPheurGetData(heur);

   // get nsols in solution storage, we cannot check whether nsols > heurdata->nknownsols as solutions from storage get removed frequently
   int nsols = SCIPgetNSols(scip);

   SCIP_PROBDATA* probdata;
   probdata = SCIPgetProbData(scip);
   assert( probdata != nullptr );

   size_t n = probdata->n;
   if ( n == 0 )
      return SCIP_OKAY;

   SCIP_SOL** init_sols = SCIPgetSols(scip);
   for ( int s = 0; s < nsols; ++s )
   {
      SCIP_SOL* init_sol = init_sols[s];
      int index = SCIPsolGetIndex(init_sol);

      if ( index > heurdata->maxindex )
         heurdata->maxindex = index;

      // check whether tabu search was executed on init_sol before
      if ( isIndexKnown(heurdata, index) )
         continue;

      // if tabuSearch finds solution, we don't want the result pointer to be set to DIDNOTFIND again
      if ( *result == SCIP_DIDNOTRUN )
      {
         SCIPdebugMsg(scip, "Running tabu search heuristic ...\n");
         *result = SCIP_DIDNOTFIND;
      }

      // store solution to known solutions, and reallocate memory if needed
      if ( heurdata->nknownsols == heurdata->knownindexsize )
      {
         int newsize = 2 * heurdata->knownindexsize;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &heurdata->knownindex, heurdata->knownindexsize, newsize) );
         heurdata->knownindexsize = newsize;
      }
      heurdata->knownindex[heurdata->nknownsols] = index;
      ++heurdata->nknownsols;

      SCIP_CALL( tabuSearch(scip, probdata, heur, init_sol, result) );
   }
   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeTabu)
{
   assert( scip != nullptr );
   assert( heur != nullptr );

   SCIP_HEURDATA* heurdata;
   heurdata = SCIPheurGetData(heur);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\nNumber of solutions seen by tabu heuristic:\t %d\n", heurdata->nknownsols);

   if ( heurdata != nullptr )
   {
      SCIPfreeRandom(scip, &(heurdata->randnumgen));
      SCIPfreeBlockMemoryArray(scip, &heurdata->knownindex, heurdata->knownindexsize);
      SCIPfreeBlockMemory(scip, &heurdata);
   }

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates tabu heuristic and includes it in SCIP */
SCIP_RETCODE BACSincludeHeurTabu(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata = nullptr;
   SCIP_HEUR* heur = nullptr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   heurdata->maxindex = -1;
   heurdata->nknownsols = 0;
   heurdata->knownindex = nullptr;
   heurdata->knownindexsize = INIT_INDEX_SIZE;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->knownindex, heurdata->knownindexsize) );
   SCIP_CALL( SCIPcreateRandom(scip, &(heurdata->randnumgen), 0, TRUE) );

   SCIP_CALL(SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecTabu, heurdata));

   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeTabu) );

   assert( heur != nullptr );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "heuristics/" HEUR_NAME "/globalbounds",
      "use global variable bounds",
      &heurdata->globalbounds, FALSE, DEFAULT_GLOBALBOUNDS, nullptr, nullptr) );

      SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/randomtabu",
         "use global variable bounds",
         &heurdata->randomtabu, FALSE, DEFAULT_RANDOMTABU, nullptr, nullptr) );

      SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/tabuoutnode",
         "use global variable bounds",
         &heurdata->tabuoutnode, FALSE, DEFAULT_TABUOUTNODE, nullptr, nullptr) );

   return SCIP_OKAY;
}
