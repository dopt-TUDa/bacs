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

/**@file   type_symmetry.h
 * @brief  type definitions for symmetry
 * @author Marc Pfetsch
 */

#ifndef TYPE_SYMMETRY_H
#define TYPE_SYMMETRY_H

/** different types of orbit selection rule */
enum ORBIT_Rule
{
   ORBIT_RULE_MAXORBIT         = 0,                 //!< rule to select orbit of maximum size
   ORBIT_RULE_MINORBIT         = 1,                 //!< rule to select orbit of minimum size
   ORBIT_RULE_STRINGENT        = 2,                 //!< rule to select orbits to generate stringent SST cuts
   ORBIT_RULE_FIRST            = 3,                 //!< select the first orbit
   ORBIT_RULE_DEGREE           = 4,                 //!< rule to select orbits with node of highest degree
   ORBIT_RULE_STRINGENT_DEGREE = 5,                 //!< rule to select orbits to generate stringent SST cuts where nodes with maximum degree are prioritized
   ORBIT_RULE_MAXDEL           = 6,                 //!< rule to select leader such that greedily most variables are fixed to 0
   ORBIT_RULE_WSTRINGENT_MAXDEL= 7,                 //!< rule to select orbits to generate weakly stringent SST cuts where nodes that cause maximal variable fixing are prioritized
   ORBIT_RULE_STRINGENT_MAXDEL = 8,                 //!< rule to select orbits to generate stringent SST cuts where nodes that cause maximal variable fixing are prioritized
   ORBIT_RULE_WSTRINGENT       = 9                  //!< rule to select orbits to generate weakly stringent SST cuts
};
typedef enum ORBIT_Rule ORBIT_RULE;

#endif
