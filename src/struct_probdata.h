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

/**@file   struct_probdata.h
 * @brief  definition of problem data
 * @author Marc Pfetsch
 */

#ifndef STRUCT_PROBDATA_H
#define STRUCT_PROBDATA_H

#include <scip/def.h>
#include <tclique/tclique.h>   // definition of clique data structures
#include "graph.h"
#include "type_graphpresolve.h"


//! BACS problem data
struct SCIP_ProbData
{
   Graph*                Gorig;              //!< original graph
   Graph*                G;                  //!< presolved graph
   SCIP_Bool             changedpresol;      //!< whether graph was changed in presolving (then G != Gorig)
   BACS_NODEFIXING*      origfixings;        //!< fixings of nodes in the original graph
   int*                  origtopresol;       //!< map from original graph to presolved graph (-1 if unmapped)
   int*                  origtopresolparity; //!< parity of nodes that are contracted to reconstruct solution in original graph
   SCIP_DISJOINTSET*     mergemapping;       //!< disjoint set structure such that representants are the nodes in the presolved graph
   SCIP_Real             objoffset;          //!< objective offset
   size_t                norig;              //!< number of nodes of Gorig
   size_t                morig;              //!< number of edges of Gorig
   size_t                n;                  //!< number of nodes of G
   size_t                m;                  //!< number of edges of G
   SCIP_Bool             unweighted;         //!< if all weights of G are 1
   std::vector<bool>*    isolated;           //!< marks whether a node is isolated
   Graph*                CG;                 //!< complemented graph
   SCIP_VAR**            vars;               //!< binary variables for each node

   SCIP_Bool             printprobstats;     //!< Print problem statistics?
   SCIP_Bool             graphpresolving;    //!< Perform graph presolving?
   SCIP_Bool             graphinprobing;     //!< Perform graph in-probing?

   SCIP_Bool             neighborhoodpresol; //!< Perform graph neighborhood presolving?
   SCIP_Bool             simplicialpresol;   //!< Perform graph simplicial presolving?
   SCIP_Bool             mergepresol;        //!< Perform merge presolving?
   SCIP_Bool             cyclepresol;        //!< Perform graph cycle presolving?
   SCIP_Bool             degtwocontract;     //!< Perform contraction of degree 2 paths?
   SCIP_Bool             sstpresolving;      //!< Perfrom SST presolving?
   SCIP_Bool             sstaddedges;        //!< Add edges during SST presolving?
   SCIP_Bool             sstrepeat;          //!< Repeat SST presolving if successful?
   SCIP_Bool             onlypresol;         //!< Only perform graph presolving?
   int                   sstorbitrule;       //!< SST orbit selection (0: maximum length, 1: minimum length, 2: stringent, 3: first variable, 4: maximum degree, 5: stringent degree, 6: max adjacent followers, 7: weakly stringent)"

   int*                  degrees;            //!< degrees of every vertex for global graph
   int*                  localdegrees;       //!< local degrees of every vertex for global graph considering local bounds
   size_t                nlocalones;         //!< number of variables locally fixed to 1.0
   size_t                nlocalzeros;        //!< number of variables locally fixed to 0.0
   size_t                nlocaledges;        //!< number of edges present in locally unfixed graph
   SCIP_Real             localdensity;       //!< density of graph induced by local bounds
   SCIP_EVENTHDLR*       eventhdlr;          //!< event handler

   int*                  cliqueindex;        //!< clique partition, giving every nonfixed vertex a clique index
   size_t*               cliquesizes;        //!< size of each clique
   size_t                ncliques;           //!< current number of cliques in partition
   size_t                maxcliqueindex;     //!< maximal index of cliques in partition
   SCIP_NODE*            cliquenode;         //!< node, where clique partition was computed
};
// typedef struct SCIP_ProbData SCIP_PROBDATA;

#endif
