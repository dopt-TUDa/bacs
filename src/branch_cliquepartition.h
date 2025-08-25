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

/**@file   branch_cliquepartition.h
 * @brief  branching on a clique partition (a.k.a. Balas-Yu) rule
 * @author Annika Jaeger
 * @author Erik Jansen
 * @author Jonas Alker
 * @author Marc Pfetsch
 */

#ifndef __SCIP_BRANCH_CLIQUEPARTITION_H__
#define __SCIP_BRANCH_CLIQUEPARTITION_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the clique partition branching rule and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE BACSincludeBranchruleCliquePartition(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
