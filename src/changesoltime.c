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

/**@file   changesoltime.c
 * @brief  Allow to change the solving time
 * @author Marc Pfetsch
 */

#include "changesoltime.h"

#include "scip/scip_timing.h"
#include "scip/clock.h"
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/struct_clock.h"


/** change solving time */
void changeSolvingTime(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_Real             newtime             /**< new solving time */
   )
{
   SCIPclockSetTime(scip->stat->solvingtime, newtime);
}
