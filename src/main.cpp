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

/**@file   main.cpp
 * @brief  main file for BACS - branch-and-cut for the stable set problem
 * @author Marc Pfetsch
 */

#include "probdata_bacs.h"
#include <scip/scipdefplugins.h>
#include "changesoltime.h"
#include "cons_clique.h"
#include "sepa_neigh.h"
#include "presol_persistence.h"
#include "presol_stableprobing.h"
#include "prop_comp.h"
#include "prop_neighborhoods.h"
#include "prop_dominance.h"
#include "prop_simplicial.h"
#include "prop_lowdegreenodes.h"
#include "prop_cliquefixing.h"
#include "branch_degree.h"
#include "branch_maxlpneigh.h"
#include "branch_cliquepartition.h"
#include "heur_greedylp.h"
#include "heur_greedydeg.h"
#include "heur_greedyrounding.h"
#include "heur_tabu.h"
#include "heur_dynamicdeg.h"
#include "heur_tabuseq.h"

#include <sys/times.h>

//! macro to check for a SCIP error and possibly exit
#define SCIP_CALL_ERROR(x)   do                                                                               \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             SCIPprintError(_restat_);                                                        \
                             return -1;                                                                       \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )


//! read comand line arguments
static
SCIP_Bool readArguments(
   SCIP*                 scip,               //!< SCIP pointer
   int                   argc,               //!< number of shell parameters
   const char**          argv,               //!< array with shell parameters
   const char**          filename,           //!< file name for network from arguments
   const char**          settingsname,       //!< name of settings file
   const char**          solutionfile,       //!< name of solution file
   SCIP_Real&            timelimit,          //!< time limit read from arguments
   SCIP_Real&            memlimit,           //!< memory limit read from arguments
   SCIP_Longint&         nodelimit,          //!< node limit read from arguments
   int&                  dispfreq,           //!< display frequency
   int&                  verblvl,            //!< display verbosity
   SCIP_Real&            cutoffvalue         //!< cutoff value
   )
{
   char usage[5000];
   int i;

   assert( scip != nullptr );
   assert( argv != nullptr );
   assert( filename != nullptr );
   assert( settingsname != nullptr );
   assert( solutionfile != nullptr );

   // init usage text
   (void) snprintf(usage, (unsigned long)5000,
      "usage: <filename> [-s <setting file>] [-l <solution file>] [-t <time limit>] [-m <mem limit>] [-n <node limit>] [-d <display frequency>] [-v <display verbosity>] [-c cutoff]");

   assert( strlen(usage) < 5000 );

   // init arguments
   *filename  = nullptr;
   *settingsname = nullptr;
   *solutionfile = nullptr;
   timelimit = SCIPinfinity(scip);
   memlimit = SCIPinfinity(scip);
   nodelimit = SCIP_LONGINT_MAX;
   dispfreq = -1;
   verblvl = -1;
   cutoffvalue = -SCIPinfinity(scip);

   if ( argc <= 1 )
   {
      SCIPinfoMessage(scip, nullptr, "%s\n", usage);
      return FALSE;
   }

   if ( *(argv[1]) == '-' )
   {
      SCIPinfoMessage(scip, nullptr, "Error: Filename should not start with '-'.\n\n");
      SCIPinfoMessage(scip, nullptr, "%s\n", usage);
      return FALSE;
   }

   *filename = argv[1];

   // check all remaining arguments
   for (i = 2; i < argc; ++i)
   {
      // check for settings file
      if ( strcmp(argv[i], "-s") == 0 )
      {
         ++i;  /*lint !e850*/
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "Error: No setting file name supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         *settingsname = argv[i];
      }
      // check for solution file
      else if ( strcmp(argv[i], "-l") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "Error: No solution file name supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         *solutionfile = argv[i];
      }
      // check for time limit
      else if ( strcmp(argv[i], "-t") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "No time limit supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         timelimit = atof(argv[i]);
      }
      // check for memory limit
      else if ( strcmp(argv[i], "-m") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "No memory limit supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         memlimit = atof(argv[i]);
      }
      // check for node limit
      else if ( strcmp(argv[i], "-n") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "No node limit supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         nodelimit = atol(argv[i]);
      }
      // check for display frequency
      else if ( strcmp(argv[i], "-d") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "No display frequency supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         dispfreq = atoi(argv[i]);
      }
      // check for display verbosity
      else if ( strcmp(argv[i], "-v") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "No verbosity level supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         verblvl = atoi(argv[i]);
      }
      // check for cutoff value
      else if ( strcmp(argv[i], "-c") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == nullptr )
         {
            SCIPinfoMessage(scip, nullptr, "No cutoff value supplied.\n\n");
            SCIPinfoMessage(scip, nullptr, "%s\n", usage);
            return FALSE;
         }
         cutoffvalue = atof(argv[i]);
      }
      else
      {
	SCIPinfoMessage(scip, nullptr, "Unknown argument %s.\n\n", argv[i]);
	SCIPinfoMessage(scip, nullptr, "%s\n", usage);
        return FALSE;
      }
   }

   return TRUE;
}


//! main function for solving the stable set problem
int main(
   int                   argc,               //!< number of command line arguments
   const char**          argv                //!< array of command line arguments
   )
{
   // initialize SCIP
   SCIP* scip;
   SCIP_CALL_ERROR( SCIPcreate(&scip) );

   // output BACS banner
   SCIPinfoMessage(scip, nullptr, "BACS: branch-and-cut for the stable set problem\n");
   const char* filename;
   const char* settingsname;
   const char* solutionfile;
   SCIP_Real timelimit;
   SCIP_Real memlimit;
   SCIP_Longint nodelimit;
   int dispfreq;
   int verblvl;
   SCIP_Real cutoffvalue;

   if ( ! readArguments(scip, argc, argv, &filename, &settingsname, &solutionfile, timelimit, memlimit, nodelimit, dispfreq, verblvl, cutoffvalue) )
   {
      SCIP_CALL_ERROR( SCIPfree(&scip) );
      BMScheckEmptyMemory();
      return -1;
   }

   // include basic plugins
   SCIP_CALL_ERROR( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL_ERROR( BACSincludeConshdlrClique(scip) );
   SCIP_CALL_ERROR( BACSincludeSepaNeigh(scip) );
   SCIP_CALL_ERROR( BACSincludePresolPersistence(scip) );
   SCIP_CALL_ERROR( BACSincludePresolStableprobing(scip) );
   SCIP_CALL_ERROR( BACSincludePropComp(scip) );
   SCIP_CALL_ERROR( BACSincludePropNeighborhoods(scip) );
   SCIP_CALL_ERROR( BACSincludePropDominance(scip) );
   SCIP_CALL_ERROR( BACSincludePropSimplicial(scip) );
   SCIP_CALL_ERROR( BACSincludePropLowdegreenodes(scip) );
   SCIP_CALL_ERROR( BACSincludePropCliquefixing(scip) );
   SCIP_CALL_ERROR( BACSincludeBranchruleDegree(scip) );
   SCIP_CALL_ERROR( BACSincludeBranchruleMaxlpneigh(scip) );
   SCIP_CALL_ERROR( BACSincludeBranchruleCliquePartition(scip) );
   SCIP_CALL_ERROR( BACSincludeHeurGreedyLP(scip) );
   SCIP_CALL_ERROR( BACSincludeHeurGreedyDeg(scip) );
   SCIP_CALL_ERROR( BACSincludeHeurGreedyRounding(scip) );
   SCIP_CALL_ERROR( BACSincludeHeurTabu(scip) );
   SCIP_CALL_ERROR( BACSincludeHeurDynamicdeg(scip) );
   SCIP_CALL_ERROR( BACSincludeHeurTabuseq(scip) );

   // output SCIP banner
   SCIPprintVersion(scip, nullptr);
   SCIPinfoMessage(scip, nullptr, "\n");
   SCIPprintExternalCodes(scip, nullptr);
   SCIPinfoMessage(scip, nullptr, "\n");

   // init problem data (allocate memory)
   SCIP_PROBDATA* probdata;
   SCIP_CALL_ERROR( BACSinitProblem(scip, &probdata) );
   assert( probdata != nullptr );

   // check for parameters
   if ( settingsname != nullptr )
   {
      if ( SCIPfileExists(settingsname) )
      {
         SCIPinfoMessage(scip, nullptr, "Reading parameter file <%s> ...\n\n", settingsname);
         SCIP_CALL_ERROR( SCIPreadParams(scip, settingsname) );
      }
      else
      {
         SCIPerrorMessage("Parameter file <%s> not found.\n", settingsname);
         return SCIP_NOFILE;
      }
   }

   // print changed paramters
   SCIPinfoMessage(scip, nullptr, "\nChanged settings:\n");
   SCIP_CALL_ERROR( SCIPwriteParams(scip, nullptr, FALSE, TRUE) );
   SCIPinfoMessage(scip, nullptr, "\n");

   // read graph
   struct tms timer_beg;
   struct tms timer_end;

   (void) times(&timer_beg);
   SCIP_CALL_ERROR( BACSreadGraph(scip, probdata, argv[1]) );
   (void) times(&timer_end);

   double t = (timer_end.tms_utime - timer_beg.tms_utime) / (double)sysconf(_SC_CLK_TCK);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Reading time: %.2f\n", t);

   // perform graph presolving
   (void) times(&timer_beg);
   SCIP_CALL_ERROR( BACSgraphPresolving(scip, probdata, timelimit) );
   (void) times(&timer_end);
   double presolvingtime = (timer_end.tms_utime - timer_beg.tms_utime) / (double)sysconf(_SC_CLK_TCK);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Presolving time: %.2f\n", presolvingtime);

   if ( probdata->onlypresol )
      return 0;

   // initialize data structures (tclique graph etc.)
   SCIP_CALL_ERROR( BACSinitDatastructures(scip, probdata) );

   // setup problem (create variables and constraints)
   SCIP_CALL_ERROR( BACSsetupProblem(scip, "BACS", probdata) );

   // set time, node, and memory limit
   if ( ! SCIPisInfinity(scip, timelimit) )
   {
      t = MAX(0.0, timelimit - presolvingtime);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Reducing time limit to: %.2f\n", t);
      SCIP_CALL_ERROR( SCIPsetRealParam(scip, "limits/time", t) );
   }

   if ( ! SCIPisInfinity(scip, memlimit) )
   {
      SCIP_CALL_ERROR( SCIPsetRealParam(scip, "limits/memory", memlimit) );
   }

   if ( nodelimit < SCIP_LONGINT_MAX )
   {
      SCIP_CALL_ERROR( SCIPsetLongintParam(scip, "limits/nodes", nodelimit) );
   }

   if ( dispfreq >= 0 )
   {
      SCIP_CALL_ERROR( SCIPsetIntParam(scip, "display/freq", dispfreq) );
   }

   if ( verblvl >= 0 && verblvl <= 5 )
   {
      SCIP_CALL_ERROR( SCIPsetIntParam(scip, "display/verblevel", verblvl) );
   }

   // set cutoff and turn heuristics off
   if ( ! SCIPisInfinity(scip, -cutoffvalue) )
   {
      SCIPinfoMessage(scip, nullptr, "set limits objective %g\n", cutoffvalue);
      SCIP_CALL_ERROR( SCIPsetObjlimit(scip, cutoffvalue) );
      SCIPinfoMessage(scip, nullptr, "Turn all heuristics off!\n");
      SCIP_CALL_ERROR( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   }

   // some parameter changes
   SCIP_CALL_ERROR( SCIPsetIntParam(scip, "display/verblevel", 5) );

   // deactivate diving heuristics
   SCIP_HEUR** heurs = SCIPgetHeurs(scip);
   int nheurs = SCIPgetNHeurs(scip);
   for (int i = 0; i < nheurs; ++i)
   {
      const char* heurname = SCIPheurGetName(heurs[i]);
      if ( strstr(heurname, "diving") != nullptr )
         SCIPheurSetFreq(heurs[i], -1);
   }

   // deactivate probing, because it will never find something
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/probing/maxprerounds", 0) );

   // possibly read solution
   if ( solutionfile != nullptr )
   {
      SCIPinfoMessage(scip, nullptr, "Reading solution file <%s> ...\n", solutionfile);
      SCIP_CALL_ERROR( SCIPreadSol(scip, solutionfile) );
   }
   else
   {
      SCIPinfoMessage(scip, nullptr, "\nOriginal problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints.\n",
         SCIPgetNOrigVars(scip), SCIPgetNOrigBinVars(scip), SCIPgetNOrigIntVars(scip), SCIPgetNOrigImplVars(scip), SCIPgetNOrigContVars(scip), SCIPgetNOrigConss(scip));
   }

   // enable debugging
#ifdef WITH_DEBUG_SOLUTION
   SCIPenableDebugSol(scip);
#endif

   // correct solving time; do not include time for reading graph
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "\nSetting solving time to %.2f.\n\n", presolvingtime);
   changeSolvingTime(scip, presolvingtime);

   // solve problem
   SCIP_CALL_ERROR( SCIPsolve(scip) );
   SCIP_CALL_ERROR( SCIPprintStatistics(scip, nullptr) );

   // check whether final solution is feasible
   SCIPinfoMessage(scip, nullptr, "\n");
   SCIP_CALL_ERROR( BACScheckBestSol(scip, probdata) );
   SCIP_CALL_ERROR( BACScheckBestSolOrig(scip, probdata) );

   // free SCIP
   // SCIPprintMemoryDiagnostic(scip);
   SCIP_CALL_ERROR( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return 0;
}
