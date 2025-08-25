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

/**@file   graph.h
 * @brief  BOOST data structure definitions
 * @author Marc Pfetsch
 */


#ifndef GRAPH_H
#define GRAPH_H

// turn of warning for boost (just for this file)
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#pragma GCC diagnostic warning "-Wzero-as-null-pointer-constant"
#pragma GCC diagnostic warning "-Wshadow"



// define properties
enum vertex_weight_t  { vertex_weight };          //!< weight of node

// tell boost about these properties
namespace boost
{
   BOOST_INSTALL_PROPERTY(vertex, weight);
}

typedef boost::property<vertex_weight_t, double> VertexWeightProperty;


/**@class Graph
 * @brief an undirected graph
 */
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexWeightProperty, boost::no_property> Graph;

// useful typedefs
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor   Edge;
typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;
typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIterator;
typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIterator;

#endif
