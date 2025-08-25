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

/**@file   cons_clique.h
 * @brief  Constraint handler class for clique inequalities
 * @author Marc Pfetsch
 */

#ifndef CONS_CLIQUE_H
#define CONS_CLIQUE_H

#include <vector>
#include <queue>
#include <scip/scip.h>
#include <tclique/tclique.h>   // definition of clique data structures
#include "graph.h"

//! creates a clique constraint
SCIP_EXPORT
SCIP_RETCODE BACScreateConsClique(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS**           cons,               /**< pointer to constraint (output) */
   const char*           name,               /**< name of constraint  */
   const Graph*          G,                  /**< underlying graph */
   const std::vector<bool>* isolated,        /**< marks whether a node is isolated */
   unsigned int          n,                  /**< number of nodes */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_Bool             initial,            /**< whether constraint is initial */
   SCIP_Bool             separate,           /**< whether constraint should be separated */
   SCIP_Bool             enforce,            /**< whether constraint should be enforced */
   SCIP_Bool             check,              /**< whether constraint should be checked */
   SCIP_Bool             propagate           /**< whether constraint should be propagated */
   );

//! creates the handler for clique constraints and includes it in SCIP
SCIP_EXPORT
SCIP_RETCODE BACSincludeConshdlrClique(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
