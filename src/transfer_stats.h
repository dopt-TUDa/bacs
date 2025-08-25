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

/**@file   transfer_stats.h
 * @brief  code to transfer statistics from one SCIP to another one
 * @author Marc Pfetsch
 */

#ifndef TRANSFER_STATS_H
#define TRANSFER_STATS_H

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** transfers statistics from sourcescip to targetscip */
SCIP_EXPORT
SCIP_RETCODE SCIPtransferStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
