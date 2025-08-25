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

/**@file   sepa_neigh.h
 * @brief  Separator for neighborhood inequalities
 * @author Marc Pfetsch
 */

#ifndef SEPA_NEIGH_H
#define SEPA_NEIGH_H

#include <scip/scip.h>

//! creates the separator for neighborhood constraints and includes it in SCIP
SCIP_EXPORT
SCIP_RETCODE BACSincludeSepaNeigh(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
