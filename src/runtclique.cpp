/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is part of the program colorbitopt                           */
/*                                                                           */
/*    an implementation of a branch-and-cut algorithm to solve the           */
/*    coloring problem by symmetry breaking methods based on orbitopes.      */
/*                                                                           */
/*    Copyright (C) 2005-2018  Marc Pfetsch                                  */
/*                                                                           */
/*                                                                           */
/*    Based on SCIP  --- Solving Constraint Integer Programs                 */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*       mailto: scip@zib.de                                                 */
/*       SCIP is distributed under the terms of the SCIP Academic Licence,   */
/*       see file COPYING in the SCIP distribution.                          */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   runtclique.cpp
 * @brief  compute maximum stable set by calling tclique on complemented graph
 * @author Marc Pfetsch
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include <vector>
#include <string.h>
#include "tclique/tclique.h"   // def. of clique data structures
#include "scip/pub_fileio.h"
#include "graph.h"

//! get next number from string s
static
SCIP_Longint getNextNumber(
   char*&                s                   //!< string pointer (modified)
   )
{
   // skip whitespace
   while ( isspace(*s) )
      ++s;

   // get number
   SCIP_Longint tmp = std::atol(s);

   // skip number
   while ( (*s != 0) && (! isspace(*s)) )
      ++s;

   return tmp;
}


/** Read graph */
void readGraph(
   const char*           filename,           /**< file name */
   TCLIQUE_GRAPH*&       TCG                 /**< graph (on exit) */
   )
{
   // open file
   SCIP_FILE* file = SCIPfopen(filename, "r");
   if ( ! file )
   {
      fprintf(stderr, "Error: could not open %s\n", filename);
      return;
   }

   printf("Reading file %s ...\n", filename);

   // get buffer for reading input lines
   char* buffer = new char [SCIP_MAXSTRLEN];

   // get number of nodes and edges: find line starting with 'p'
   unsigned int line = 1;
   char* str = buffer;
   do
   {
      (void) SCIPfgets(buffer, SCIP_MAXSTRLEN, file);
      ++line;

      // skip whitespace
      str = buffer;
      while ( isspace(*str) )
         ++str;
   }
   while ( ! SCIPfeof(file) && *str != 'p' );

   if ( SCIPfeof(file) )
   {
      fprintf(stderr, "%s: Could not find line starting with 'p'.\n", filename);
      delete [] buffer;
      (void) SCIPfclose(file);
      abort();
   }

   // skip 'p'
   ++str;

   // skip whitespace
   while ( isspace(*str) )
      ++str;

   // check whether next string part is "edge"
   if ( strncmp(str, "edge", 4UL) == 0 )
   {
      // skip 'edge'
      str = str + 5;

      // if line reads 'edges' (non-standard!), instead of 'edge'.
      if ( *str == 's' )
         ++str;
   }
   else
   {
      // check nonstandard alternative 'col' instead of 'edge'
      if ( strncmp(str, "col", 3UL) == 0 )
      {
         // skip 'col'
         str = str + 4;
      }
      else
      {
         fprintf(stderr, "%s: line %u: Line starting with 'p' must continue with 'edge' or 'col'!\n%s\n", filename, line, buffer);
         delete [] buffer;
         (void) SCIPfclose(file);
         abort();
      }
   }

   // read sizes
   SCIP_Longint nn = getNextNumber(str);
   SCIP_Longint mm = getNextNumber(str);

   if ( nn <= 0 )
   {
      fprintf(stderr, "%s: line %u: Number of nodes must be positive!\n%s\n", filename, line, buffer);
      delete [] buffer;
      (void) SCIPfclose(file);
      abort();
   }

   if ( mm < 0 )
   {
      fprintf(stderr, "%s: line %u: Number of edges must be nonnegative!\n%s\n", filename, line, buffer);
      delete [] buffer;
      (void) SCIPfclose(file);
      abort();
   }
   size_t n = (size_t) nn;
   size_t m = (size_t) mm;

   // set to store edges to avoid duplicates
   std::set<std::pair<size_t, size_t> > EdgeSet;

   // start to read nodes and edges into set
   (void) SCIPfgets(buffer, SCIP_MAXSTRLEN, file);
   ++line;
   assert( ! SCIPfeof(file) );

   // detect isolated nodes, loops, and duplicate edges
   std::vector<SCIP_Real> weights(n, 0.0); // default value for the weighted version is 0
   bool foundweights = false;
   size_t nloops = 0;
   size_t nduplicateedges = 0;

   // loop through file
   while ( ! SCIPfeof(file) )
   {
      str = buffer;

      if ( str[0] == 'n' )
      {
         ++str;
         SCIP_Longint node = getNextNumber(str) - 1;
         while ( isspace(*str) )
            ++str;
         SCIP_Real weight = atof(str);

         // check node
         if ( node < 0 || node >= nn )
         {
            fprintf(stderr, "%s: line %u: Node number invalid!\n%s\n", filename, line, buffer);
            delete [] buffer;
            (void) SCIPfclose(file);
            abort();
         }
         weights[node] = weight;
         foundweights = true;
      }
      else if ( str[0] == 'e' )
      {
         ++str;
         SCIP_Longint source = getNextNumber(str) - 1;
         SCIP_Longint target = getNextNumber(str) - 1;

         // check source
         if ( source < 0 || source >= nn || target < 0 || target >= nn )
         {
            fprintf(stderr, "%s: line %u: Node number invalid!\n%s\n", filename, line, buffer);
            delete [] buffer;
            (void) SCIPfclose(file);
            abort();
         }
         size_t s = (size_t) source;
         size_t t = (size_t) target;

         // detect self-loops
         if ( s != t )
         {
            std::pair<std::set<std::pair<size_t, size_t> >::iterator, bool> ret;
            if ( s < t )
               ret = EdgeSet.insert(std::make_pair(s, t));
            else
               ret = EdgeSet.insert(std::make_pair(t, s));

            // if edge is already present just count it
            if ( ! ret.second )
               ++nduplicateedges;
         }
         else
            ++nloops;
      }

      // read next line
      (void) SCIPfgets(buffer, SCIP_MAXSTRLEN, file);
      ++line;
   }

   (void) SCIPfclose(file);
   delete [] buffer;

   // create nodes
   Graph G;
   for (size_t i = 0; i < n; ++i)
   {
      Vertex v = boost::add_vertex(G);

      if ( foundweights )
      {
         boost::put(vertex_weight_t(), G, v, weights[v]);
      }
      else
         boost::put(vertex_weight_t(), G, v, 1);  // default weight is 1 for the unweighted version
   }
   assert( boost::num_vertices(G) == n );

   if ( nloops > 0 )
   {
      printf("Graph contains %zu loops (ignored)!\n", nloops);
   }

   // create edges
   size_t nedges = 0;
   std::set<std::pair<size_t, size_t> >::const_iterator eit, eend = EdgeSet.end();
   for (eit = EdgeSet.begin(); eit != eend; ++eit)
   {
      size_t s = eit->first;
      size_t t = eit->second;

      assert( ! boost::edge(s, t, G).second );

      std::pair<Edge, bool> p = boost::add_edge(s, t, G);
      ++nedges;
      if ( ! p.second )
      {
         fprintf(stderr, "Could not create edge {%lu,%lu}.\n", s, t);
         abort();
      }
   }
   assert( boost::num_edges(G) == nedges );

   if ( nduplicateedges > 0 )
      printf("Duplicate edges: %zu\n", nduplicateedges);

   if ( nedges + nduplicateedges + nloops != m )
   {
      fprintf(stderr, "%s: line %u: Found %zu edges. There should be %zu.\n", filename, line, nedges + nduplicateedges + nloops, m);
      return;
   }

   printf("Number of nodes:\t\t %zu\n", n);
   printf("Number of edges:\t\t %zu\n", nedges);

   // create tclique graph
   if ( ! tcliqueCreate(&TCG) )
      abort();

   // add vertices
   for (unsigned int i = 0; i < n; ++i)
   {
      if ( ! tcliqueAddNode(TCG, (int) i, 0) )
	 abort();
   }

   assert( tcliqueGetNEdges(TCG) == (int) (2 * (nedges - nduplicateedges) ) );
   assert( tcliqueGetNNodes(TCG) == (int) n );

   delete [] buffer;
}




/* main function */
int main(int argc, const char** argv)
{
   if (argc != 2)
   {
      std::cerr << "usage: " << argv[0] << " <file>" << std::endl;
      exit(1);
   }

   TCLIQUE_GRAPH* TCG;

   // char probname[10];
   // tcliqueLoadFile(&TCG, argv[1], 1.0, probname, 10);

   readGraph(argv[1], TCG);
   int N = tcliqueGetNNodes(TCG);
   int M = tcliqueGetNEdges(TCG);

   // initialize weights
   for (int i = 0; i < N; ++i)
      tcliqueChangeWeight(TCG, i, 1);

   // call the clique algorithm ...
   int ncliquenodes = 0;
   int* cliquenodes = new int [N]; /*lint !e737*/
   TCLIQUE_WEIGHT weight = 0;
   TCLIQUE_STATUS status;
   // 3rd to last parameter: allow to generate cliques containing 0 weight nodes (not necessary in this case)
   tcliqueMaxClique(tcliqueGetNNodes, tcliqueGetWeights, tcliqueIsEdge, tcliqueSelectAdjnodes,
      TCG, nullptr, nullptr, cliquenodes, &ncliquenodes, &weight,
      N, 0, 10000000, 0, (int) N, -1, nullptr, &status);

   std::cout << " " << std::setw(4) << N << " & " << std::setw(5) << M/2 << " & ";
   std::cout << ncliquenodes;
   if (status != TCLIQUE_OPTIMAL)
      std::cout << "*";
   std::cout << std::endl;

   delete [] cliquenodes;
   tcliqueFree(&TCG);
}
