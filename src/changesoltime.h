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

/**@file   changesoltime.h
 * @brief  Allow to change the solving time
 * @author Marc Pfetsch
 */

#ifndef CHANGESOLTIME_H
#define CHANGESOLTIME_H

#include "scip/def.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** change solving time */
void changeSolvingTime(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_Real             newtime             /**< new solving time */
   );

#ifdef __cplusplus
}
#endif

#endif
