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

/**@file   symmetry.h
 * @brief  Methods to deal with symmetries
 * @author Marc Pfetsch
 *
 * This code is based on the implementation in maxkcol by Christopher Hojny.
 */

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "graph.h"
#include "scip/scip.h"
#include "type_graphpresolve.h"
#include "type_symmetry.h"

// typedefs
typedef std::vector<std::vector<Vertex> >         VecVecVertex;
typedef std::vector<std::vector<unsigned int> >   VecVecUInt;
typedef std::vector<Vertex>                       VecVertex;
typedef std::list<std::pair<Vertex,Vertex> >      VertexPairList;

//! Compute leaders and followers for SST cuts heuristically by filtering generators
SCIP_RETCODE computeLeaderFollowersSSTgensFilter(
   SCIP*                 scip,               //!< SCIP instance
   const Graph*          G,                  //!< graph
   size_t                nelems,             //!< number of elements the permutations are acting on
   const VecVecUInt&     gens,               //!< group generators
   ORBIT_RULE            orbitrule,          //!< rule to choose next orbit
   BACS_NODEFIXING*      fixed,              //!< array for fixings of nodes
   VecVertex&            leaders,            //!< returns leaders
   VecVecVertex&         followers,          //!< returns followers
   size_t*               degrees             //!< degrees
   );

//! Compute generators of automorphism group of a graph
SCIP_RETCODE computeAutomorphismsGraph(
   SCIP*                 scip,               //!< SCIP instance
   const Graph*          G,                  //!< graph
   const Graph*          E,                  //!< additional graph
   const BACS_NODEFIXING* fixed,             //!< array to mark fixed nodes
   const SCIP_Real*      weights,            //!< current weights
   SCIP_Bool             silent,             //!< whether no output should be produced
   long double&          groupsize,          //!< group size
   VecVecUInt&           generators          //!< returns generators
   );

#endif
