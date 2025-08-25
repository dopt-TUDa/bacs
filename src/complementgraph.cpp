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

/**@file   complementgraph.cpp
 * @brief  compute the complement graph
 * @author Marc Pfetsch
 */

#include <set>
#include <vector>
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


/** create complemented graph */
static
void createComplementedGraph(
   const Graph&          G,                  /**< graph to be complemented */
   Graph&                CG                  /**< complemented graph to be created */
   )
{
   assert( boost::num_vertices(CG) == 0 );

   // create nodes
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(G); vit != vend; ++vit)
   {
      Vertex v = boost::add_vertex(CG);
      SCIP_Real weight = boost::get(vertex_weight_t(), G, *vit);
      boost::put(vertex_weight_t(), CG, v, weight);
   }

   std::vector<int> A(boost::num_vertices(G), -1);

   // add edges by brute force
   for (boost::tie(vit, vend) = boost::vertices(G); vit != vend; ++vit)
   {
      Vertex v = *vit;

      // mark neighbors
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, G); ait != aend; ++ait)
         A[*ait] = (int) v;

      // add non-edges
      VertexIterator vit2 = vit;
      for (++vit2; vit2 != vend; ++vit2)
      {
         if ( A[*vit2] != (int) v )
            (void) boost::add_edge(v, *vit2, CG);
      }
   }
}


/** Read graph */
int readGraph(
   const char*           filename,           /**< file name */
   Graph&                G,                  /**< graph (on exit) */
   bool&                 weighted            /**< whether the graph is weighted */
   )
{
   weighted = false;

   // open file
   SCIP_FILE* file = SCIPfopen(filename, "r");
   if ( ! file )
   {
      fprintf(stderr, "Error: could not open %s\n", filename);
      return 0;
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
         weighted = true;
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
   for (size_t i = 0; i < n; ++i)
   {
      Vertex v = boost::add_vertex(G);

      if ( weighted )
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
   std::set<std::pair<size_t, size_t> >::const_iterator sit, send = EdgeSet.end();
   for (sit = EdgeSet.begin(); sit != send; ++sit)
   {
      size_t s = sit->first;
      size_t t = sit->second;

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
      return 0;
   }
   return 1;
}


/** main function */
int main(int argc, const char** argv)
{
   if ( argc != 3 )
   {
      printf("usage: <dimacs infile> <dimacs outfile>\n");
      return 0;
   }

   // read graph and stort it in boost graph G
   Graph G;
   bool weighted;
   if ( ! readGraph(argv[1], G, weighted) )
   {
      return 0;
   }

   // create complemented graph
   printf("Complementing graph ...\n");
   Graph CG;
   createComplementedGraph(G, CG);

   printf("Size of complemented graph:\n");
   printf("Number of nodes: %zu\n", boost::num_vertices(CG));
   printf("Number of edges: %zu\n", boost::num_edges(CG));

   // write complemented graph
   SCIP_FILE* file = SCIPfopen(argv[2], "w");
   if ( ! file )
   {
      printf("Error: Could not open <%s>.\n", argv[2]);
      return 0;
   }

   printf("Writing complemented graph to <%s>.\n", argv[2]);
   (void) SCIPfprintf(file, "p edge %zu %zu\n", boost::num_vertices(CG), boost::num_edges(CG));
   // if the graph is weighted, we have to output the nodes
   VertexIterator vit, vend;
   for (boost::tie(vit,vend) = boost::vertices(CG); vit != vend; ++vit)
      (void) SCIPfprintf(file, "n %zu %g\n", *vit + 1, boost::get(vertex_weight_t(), CG, *vit));

   EdgeIterator eit, eend;
   for (boost::tie(eit,eend) = boost::edges(CG); eit != eend; ++eit)
   {
      (void) SCIPfprintf(file, "e %zu %zu\n", boost::source(*eit, CG) + 1, boost::target(*eit, CG) + 1);
   }
   (void) SCIPfclose(file);

   return 1;
}
