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

/**@file   heur_greedydeg.h
 * @ingroup PRIMALHEURISTICS
 * @brief  greedy degree primal heuristic with smart degree computation
 * @author Erik Jansen
 * @author Jonas Alker
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_GREEDY_DEG_H__
#define __SCIP_HEUR_GREEDY_DEG_H__

#include <scip/scip.h>

#ifdef __cplusplus
extern "C"{
#endif

//! creates the greedy degree primal heuristic and includes it in SCIP
SCIP_EXPORT
SCIP_RETCODE BACSincludeHeurGreedyDeg(
   SCIP*                 scip                //! SCIP data structure
   );

#ifdef __cplusplus
}
#endif

#endif
