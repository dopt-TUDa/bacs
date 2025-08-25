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
/*    Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                    */
/*                                                                           */
/*    Both are licensed under the Apache License, Version 2.0.               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_stableprobing.h
 * @ingroup PRESOLVERS
 * @brief  presolver to perform probing
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef PRESOL_STABLEPROBING_H
#define PRESOL_STABLEPROBING_H

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

//! creates the stable probing presolver and includes it in SCIP
SCIP_EXPORT
SCIP_RETCODE BACSincludePresolStableprobing(
   SCIP*                 scip                //! SCIP data structure
   );

#ifdef __cplusplus
}
#endif

#endif
