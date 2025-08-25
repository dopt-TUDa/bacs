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

/**@file   type_graphpresolve.h
 * @brief  definitions for graph presolving
 * @author Marc Pfetsch
 */

#ifndef TYPE_GRAPHPRESOLVE_H
#define TYPE_GRAPHPRESOLVE_H

//! enum to encode node fixings
enum BACS_Nodefixing
{
   BACSfixedzero  = -1,  /**< node is fixed to 0 (not part of the optimal stable set) */
   BACSunfixed    = 0,   /**< node is unfixed */
   BACSfixedone   = 1    /**< node is fixed to 1 (part of the optimal stable set) */
};
typedef enum BACS_Nodefixing BACS_NODEFIXING;

#endif
