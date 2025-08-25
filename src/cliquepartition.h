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

/**@file   cliquepartition.h
 * @brief  functions to compute cliquepartition on local graph
 * @author Jonas Alker
 */


#ifndef CLIQUEPARTITION_H
#define CLIQUEPARTITION_H

#include <scip/scip.h>
#include "graph.h"


/** different types of partitioning methods */
enum PARTITIONING_Method
{
   PARTITIONING_COLOR    = 0,                //!< use coloring heuristic
   PARTITIONING_LINEAR   = 1,                //!< use linear clique partition
   PARTITIONING_GREEDY   = 2,                //!< use greedy clique partition
};
typedef enum PARTITIONING_Method PARTITIONING_METHOD;

#ifndef NDEBUG
SCIP_EXPORT
void BACScheckCliquesizes(
   SCIP*                 scip,               /**< SCIP pointer */
   int*                  cand,               /**< clique partition, representant for every vertex */
   size_t*               cliquesizes,        /**< sizes of cliques */
   size_t                ngen,               /**< number of initially generated cliques */
   size_t                n                   /**< number of vertices */
   );
#endif

/** compute lowerbound on current locally valid graph respecting the current primalbound, the objective offset and locally fixed variables */
SCIP_EXPORT
SCIP_RETCODE BACScomputeLocalLowerbound(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   size_t&               lowerbound          /**< lowerbound regarding free nodes */
   );

/** try to improve given partition by removing cliques from the partition and moving contained vertices to other cliques
 *
 *  if successful, the method runs again until no further improvement can be found or lowerbound cliques have been reached
 */
SCIP_EXPORT
SCIP_RETCODE improvePartition(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int*                  cand,               /**< clique partition, representant for every vertex */
   size_t*               cliquesizes,        /**< sizes of cliques */
   size_t                ngen,               /**< number of initially generated cliques */
   size_t                lowerbound,         /**< lowerbound regarding free nodes */
   size_t&               ncliques,           /**< current number of cliques */
   SCIP_Bool&            cutoff              /**< whether a cutoff was detected */
   );

/** compute clique partition in O(n + m)
 *
 *  this partitioning algorithm is inspired by a simple coloring heuristic
 *
 *  we start with an empty partitioning
 *
 *  for all vertices v in G we check whether v can be added to an existing clique in the partitioning by scanning the neighborhood of v
 *  and checking whether a clique in the current partitioning is covered by the neighborhood of v, then we add v oder create a new clique
 *
 *  we scan every vertex exactly once and all adjacent edges at most once
 */
SCIP_EXPORT
SCIP_RETCODE computeCliquePartitionColoring(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int*                  cand,               /**< clique partition, representant for every vertex */
   size_t*               cliquesizes,        /**< sizes of cliques */
   size_t&               ngen                /**< number of initially generated cliques */
   );

/** compute clique partition in O(m)
 *
 *  start with a trivial partition of the graph G where all nodes are in the same subset
 *
 *  for all vertices v in G we split the subset containing v such that all non-adjacent vertices of v
 *  are in a different subset that v itself
 *
 *  after iterating over all vertices every vertex is connected with every other vertex in the same subset - a clique
 */
SCIP_EXPORT
SCIP_RETCODE computeCliquePartitionLinear(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int*                  cand,               /**< clique partition, representant for every vertex */
   size_t*               cliquesizes,        /**< sizes of cliques */
   size_t&               ngen                /**< number of initially generated cliques */
   );

/** compute greedy clique partition in O(n*deg^2) */
SCIP_EXPORT
SCIP_RETCODE computeCliquePartitionGreedy(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int*                  cand,               /**< clique partition, representant for every vertex */
   size_t*               cliquesizes,        /**< sizes of cliques */
   size_t&               ngen                /**< number of initially generated cliques */
   );

#endif
