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

/**@file   probdata_bacs.h
 * @brief  Problem class
 * @author Marc Pfetsch
 */

#ifndef PROBDATA_BACS_H
#define PROBDATA_BACS_H

#include <scip/scip.h>
#include "struct_probdata.h"

//! init problem data (allocate memory, init variables)
SCIP_EXPORT
SCIP_RETCODE BACSinitProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< problem data structure */
   );

//! free problem data
SCIP_EXPORT
SCIP_RETCODE BACSfreeProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< problem data structure */
   );

//! read graph
SCIP_EXPORT
SCIP_RETCODE BACSreadGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   const char*           filename            /**< name of graph file to read */
   );

//! setup problem (create variables and constraints)
SCIP_EXPORT
SCIP_RETCODE BACSsetupProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   );

//! perform graph presolving
SCIP_EXPORT
SCIP_RETCODE BACSgraphPresolving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   SCIP_Real             timelimit           /**< time limit for presolving */
   );

//! init data structures
SCIP_EXPORT
SCIP_RETCODE BACSinitDatastructures(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   );

//! create complemented graph
SCIP_EXPORT
SCIP_RETCODE BACScreateComplementedGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   );

//! debug routine to check whether the local degrees are correct
SCIP_EXPORT
SCIP_RETCODE BACScheckLocalDegrees(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   );

//! check whether best solution of SCIP is a stable set
SCIP_EXPORT
SCIP_RETCODE BACScheckBestSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   );

//! check whether best solution of SCIP corresponds to a stable set in the original graph
SCIP_EXPORT
SCIP_RETCODE BACScheckBestSolOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   );

#endif
