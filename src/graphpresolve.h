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

/**@file   graphpresolve.h
 * @brief  functions to preprocess the graph
 * @author Marc Pfetsch
 */

#ifndef GRAPHPRESOLVE_H
#define GRAPHPRESOLVE_H

#include <scip/scip.h>
#include "graph.h"
#include "type_graphpresolve.h"
#include "type_symmetry.h"


//! determine connected components through BFS, ignoring nodes whose variable is fixed
SCIP_EXPORT
SCIP_RETCODE BFSExtended(
   SCIP*                 scip,               //!< SCIP main data structure
   const Graph*          G,                  //!< graph
   const Graph*          E,                  //!< additional edges to be processed
   BACS_NODEFIXING*      fixed,              //!< array for fixings of nodes
   int*                  components,         //!< array to store the component index for each variable
   size_t&               ncomponents         //!< number of components in graph
   );

/** compute degrees */
SCIP_EXPORT
void BACScomputeDegrees(
   const Graph*          G,                  /**< graph to be processed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees to be computed */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo          /**< number of degree 2 nodes */
   );

/** neighborhood graph presolving - fast version relying on sorted adjacency lists in G */
SCIP_EXPORT
void BACSgraphNeighborhoodPresolving(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   size_t&               nfixedzero          /**< return number of nodes fixed to 0 */
   );

/** neighborhood graph presolving - version that takes additional edges in E into account */
SCIP_EXPORT
SCIP_RETCODE BACSgraphNeighborhoodPresolvingExtended(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   const Graph*          E,                  /**< additional edges to be processed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   size_t&               nfixedzero          /**< return number of nodes fixed to 0 */
   );

/** merge operation - fast version relying on sorted adjacency lists in G */
SCIP_EXPORT
SCIP_RETCODE BACSgraphMergePresolving(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   SCIP_Bool&            isweighted,         /**< is graph weighted? */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   SCIP_DISJOINTSET*     mergemapping,       /**< disjoint set structure such that representatives are the nodes in the presolved graph */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   size_t&               nmerged             /**< return number of pairs of nodes have have been merged */
   );

/** merge presolving - version that takes additional edges in E into account */
SCIP_EXPORT
SCIP_RETCODE BACSgraphMergePresolvingExtended(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   const Graph*          E,                  /**< additional edges to be processed */
   SCIP_Bool&            isweighted,         /**< is graph weighted? */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   SCIP_DISJOINTSET*     mergemapping,       /**< disjoint set structure such tha representatives are the nodes in the presolved graph */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   size_t&               nmerged             /**< return number of pairs of nodes have have been merged */
   );

/** perform simplicial presolving w.r.t. degrees 0/1 */
SCIP_EXPORT
void BACSsimplicialPresolving(
   SCIP*                 scip,               /**< SCIP data structure */
   const Graph*          G,                  /**< graph */
   const Graph*          E,                  /**< extended graph to be processed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   size_t&               nfixedzero,         /**< return number of found 0 fixings */
   size_t&               nfixedone           /**< return number of found 1 fixings */
   );

/** perform cycle presolving */
SCIP_EXPORT
SCIP_RETCODE BACScyclePresolving(
   SCIP*                 scip,               /**< SCIP data structure */
   const Graph*          G,                  /**< graph */
   const Graph*          E,                  /**< extended graph to be processed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   size_t&               nfixedzero,         /**< return number of found 0 fixings */
   size_t&               nfixedone           /**< return number of found 1 fixings */
   );

/** perform one complete round of presolving */
SCIP_EXPORT
SCIP_RETCODE BACSpresolvingRound(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   SCIP_Bool&            isweighted,         /**< is graph weighted? */
   SCIP_Bool             neighborhoodpresol, /**< whether neighborhood presolving should be performed */
   SCIP_Bool             simplicialpresol,   /**< whether simplicial presolving should be performed */
   SCIP_Bool             mergepresol,        /**< whether merge presolving should be performed */
   SCIP_Bool             cyclepresol,        /**< whether cycle presolving should be performed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   SCIP_DISJOINTSET*     mergemapping,       /**< disjoint set structure such that representatives are the nodes in the presolved graph */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   size_t&               nfixedzeroneigh,    /**< return number of found 0 fixings by neighborhood presolving */
   size_t&               nfixedzerosimpl,    /**< return number of found 0 fixings by simplicial presolving with degree 0/1 */
   size_t&               nfixedzerocycle,    /**< return number of found 0 fixings by cycle presolving */
   size_t&               nfixedonesimpl,     /**< return number of found 1 fixings by simplicial presolving with degree 0/1 */
   size_t&               nfixedonecycle,     /**< return number of found 1 fixings by cycle presolving */
   size_t&               nmerged             /**< return number of pairs of nodes have have been merged */
   );

/** perform one complete round of presolving, also using additional edges in E */
SCIP_EXPORT
SCIP_RETCODE BACSpresolvingRoundExtend(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   const Graph*          E,                  /**< extended graph to be processed */
   SCIP_Bool&            isweighted,         /**< is graph weighted? */
   SCIP_Bool             neighborhoodpresol, /**< whether neighborhood presolving should be performed */
   SCIP_Bool             simplicialpresol,   /**< whether simplicial presolving should be performed */
   SCIP_Bool             mergepresol,        /**< whether merge presolving should be performed */
   SCIP_Bool             cyclepresol,        /**< whether cycle presolving should be performed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   SCIP_DISJOINTSET*     mergemapping,       /**< disjoint set structure such that representatives are the nodes in the presolved graph */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   size_t&               nfixedzeroneigh,    /**< return number of found 0 fixings by neighborhood presolving */
   size_t&               nfixedzerosimpl,    /**< return number of found 0 fixings by simplicial presolving with degree 0/1 */
   size_t&               nfixedzerocycle,    /**< return number of found 0 fixings by cycle presolving */
   size_t&               nfixedonesimpl,     /**< return number of found 1 fixings by simplicial presolving with degree 0/1 */
   size_t&               nfixedonecycle,     /**< return number of found 1 fixings by cycle presolving */
   size_t&               nmerged             /**< return number of pairs of nodes have have been merged */
   );

/** perform rounds of presolving until no changes happen anymore, also using additional edges in E */
SCIP_EXPORT
SCIP_RETCODE BACSpresolvingRoundsExtend(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   const Graph*          E,                  /**< extended graph to be processed */
   SCIP_Bool&            isweighted,         /**< is graph weighted? */
   SCIP_Bool             neighborhoodpresol, /**< whether neighborhood presolving should be performed */
   SCIP_Bool             simplicialpresol,   /**< whether simplicial presolving should be performed */
   SCIP_Bool             mergepresol,        /**< whether merge presolving should be performed */
   SCIP_Bool             cyclepresol,        /**< whether cycle presolving should be performed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   SCIP_DISJOINTSET*     mergemapping,       /**< disjoint set structure such that representatives are the nodes in the presolved graph */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   size_t&               nrounds,            /**< counter of rounds to increase and output */
   size_t&               nfixedzeroneigh,    /**< return number of found 0 fixings by neighborhood presolving */
   size_t&               nfixedzerosimpl,    /**< return number of found 0 fixings by simplicial presolving with degree 0/1 */
   size_t&               nfixedzerocycle,    /**< return number of found 0 fixings by cycle presolving */
   size_t&               nfixedonesimpl,     /**< return number of found 1 fixings by simplicial presolving with degree 0/1 */
   size_t&               nfixedonecycle,     /**< return number of found 1 fixings by cycle presolving */
   size_t&               nmerged             /**< return number of pairs of nodes have have been merged */
   );

/** perform in-probing */
SCIP_EXPORT
SCIP_RETCODE BACSgraphInProbing(
   SCIP*                 scip,               /**< SCIP pointer */
   const struct tms&     timer_beg,          /**< timer for the start time */
   SCIP_Real             timelimit,          /**< timelimit that we want to obey */
   const Graph*          G,                  /**< graph to be processed */
   SCIP_Bool&            isweighted,         /**< is graph weighted? */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   SCIP_DISJOINTSET*     mergemapping,       /**< disjoint set structure such that representatives are the nodes in the presolved graph */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   Graph*                E,                  /**< extended graph to be processed */
   size_t&               naddededges         /**< return number of added edges */
   );

/** SST presolving */
SCIP_EXPORT
SCIP_RETCODE BACSpresolvingSST(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_Bool             sstaddedges,        /**< whether we should add edges */
   const Graph*          G,                  /**< graph to be processed */
   SCIP_Bool&            isweighted,         /**< is graph weighted? */
   ORBIT_RULE            orbitrule,          /**< rule to choose next orbit */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   SCIP_DISJOINTSET*     mergemapping,       /**< disjoint set structure such that representatives are the nodes in the presolved graph */
   int*                  components,         /**< array to store the component index for each variable */
   size_t                ncomponents,        /**< number of components in graph */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   Graph*                E,                  /**< additional graph to be processed and extended */
   size_t&               nfixedzero,         /**< return number of found 0 fixings */
   size_t&               naddededges,        /**< return number of added edges */
   size_t&               nmerged             /**< return number of pairs of nodes have have been merged */
   );

#endif
