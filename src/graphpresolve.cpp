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

/**@file   graphpresolve.cpp
 * @brief  functions to preprocess the graph
 * @author Marc Pfetsch
 */

#include "extiterators.hpp"
#include "symmetry.h"
#include <unistd.h>
#include <sys/times.h>


//! determine connected components through BFS, ignoring nodes whose variable is fixed
SCIP_RETCODE BFSExtended(
   SCIP*                 scip,               //!< SCIP main data structure
   const Graph*          G,                  //!< graph
   const Graph*          E,                  //!< additional edges to be processed
   BACS_NODEFIXING*      fixed,              //!< array for fixings of nodes
   int*                  components,         //!< array to store the component index for each variable
   size_t&               ncomponents         //!< number of components in graph
   )
{
   assert( scip != nullptr );
   assert( fixed != nullptr );
   assert( components != nullptr );

   ncomponents = 0;

   // compute components by BFS
   size_t n = boost::num_vertices(*G);
   Vertex* Q; // queue
   SCIP_CALL( SCIPallocBufferArray(scip, &Q, (int) n) );

   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
      components[*vit] = -1;

   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      assert( v < n );

      if ( components[v] >= 0 )
         continue;

      if ( fixed[v] != BACSunfixed )
         continue;

      // mark vertex
      components[v] = (int) ncomponents;

      // init queue
      Q[0] = (Vertex) v;
      size_t startq = 0;
      size_t endq = 0;

      // loop through Q until empty
      while ( startq <= endq )
      {
         Vertex u = Q[startq++];
         assert( startq <= n );

         bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(u, G, E);
         for (; ! ait.at_end(); ++ait)
         {
            Vertex w = *ait;

            // do not connect through fixed vertices
            if ( fixed[w] == BACSunfixed )
            {
               if ( components[w] < 0 )
               {
                  components[w] = (int) ncomponents;
                  Q[++endq] = w;
                  assert( endq <= n );
               }
            }
         }
      }
      ++ncomponents;
   }

   SCIPfreeBufferArray(scip, &Q);

   return SCIP_OKAY;
}


/** compute degrees */
void BACScomputeDegrees(
   const Graph*          G,                  /**< graph to be processed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees to be computed */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo          /**< number of degree 2 nodes */
   )
{
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );

   ndegreezero = 0;
   ndegreeone = 0;
   ndegreetwo = 0;

   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      assert( v < boost::num_vertices(*G) );

      if ( fixed[v] != BACSunfixed )
         degrees[v] = SIZE_MAX;
      else
      {
         degrees[v] = 0;

         AdjacencyIterator ait, aend;
         for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            Vertex w = *ait;
            assert( w < boost::num_vertices(*G) );
            if ( fixed[w] == BACSunfixed )
               ++degrees[v];
         }

         if ( degrees[v] == 0 )
            ++ndegreezero;
         else if ( degrees[v] == 1 )
            ++ndegreeone;
         else if ( degrees[v] == 2 )
            ++ndegreetwo;
      }
   }
}


#ifndef NDEBUG
/** check whether degrees are correct w.r.t. fixings */
static
SCIP_RETCODE BACScheckDegrees(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   const Graph*          E,                  /**< possible extended graph */
   const BACS_NODEFIXING*fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees to be computed */
   size_t                ndegreezero,        /**< number of degree 0 nodes */
   size_t                ndegreeone,         /**< number of degree 1 nodes */
   size_t                ndegreetwo          /**< number of degree 2 nodes */
   )
{
   assert( scip != nullptr );
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );

   size_t n = boost::num_vertices(*G);

   // compute current degrees
   size_t* mydeg;
   SCIP_CALL( SCIPallocClearBufferArray(scip, &mydeg, n) );

   bacs::ExtEdgeIterator eit = bacs::ExtEdgeIterator(G, E);
   for (; ! eit.at_end(); ++eit)
   {
      Vertex s = eit.source();
      Vertex t = eit.target();

      assert( s < n );
      assert( t < n );

      // skip edges with one fixed variable
      if ( fixed[s] != BACSunfixed || fixed[t] != BACSunfixed )
         continue;

      ++(mydeg[s]);
      ++(mydeg[t]);
   }

   // count nodes of degree 0, 1, and 2
   size_t ndegzero = 0;
   size_t ndegone = 0;
   size_t ndegtwo = 0;
   for (size_t i = 0; i < n; ++i)
   {
      if ( fixed[i] != BACSunfixed )
      {
         assert( degrees[i] == SIZE_MAX );
         assert( mydeg[i] == 0 );
      }
      else
      {
         if ( mydeg[i] == 0 )
            ++ndegzero;
         else if ( mydeg[i] == 1 )
            ++ndegone;
         else if ( mydeg[i] == 2 )
            ++ndegtwo;

         assert( degrees[i] < n );
         assert( mydeg[i] == degrees[i] );
      }
   }
   assert( ndegzero == ndegreezero );
   assert( ndegone == ndegreeone );
   assert( ndegtwo == ndegreetwo );

   SCIPfreeBufferArray(scip, &mydeg);

   return SCIP_OKAY;
}

/** check whether fixings are consistent, i.e., do not fix neighboring nodes to 1 */
static
void BACScheckConsistentFixings(
   const Graph*          G,                  /**< graph to be processed */
   const Graph*          E,                  /**< additional edges to be processed */
   const BACS_NODEFIXING*fixed               /**< array for fixings of nodes */
   )
{
   assert( G != nullptr );
   assert( fixed != nullptr );

   bacs::ExtEdgeIterator eit = bacs::ExtEdgeIterator(G, E);
   for (; ! eit.at_end(); ++eit)
   {
      Vertex s = eit.source();
      Vertex t = eit.target();
      assert( fixed[s] != BACSfixedone || fixed[t] != BACSfixedone );
   }
}
#endif


/** fix node to 0 and update data */
inline
void BACSfixNodeZero(
   const Graph*          G,                  /**< graph to be processed */
   Vertex                v,                  /**< node to be fixed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees to be updated */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo          /**< number of degree 2 nodes */
   )
{
   assert( G != nullptr );
   assert( v < boost::num_vertices(*G) );
   assert( fixed != nullptr );
   assert( fixed[v] == BACSunfixed );
   assert( degrees != nullptr );
   assert( degrees[v] < boost::num_vertices(*G) );

   // adjust degrees of neighboring nodes
   AdjacencyIterator ait, aend;
   for (boost::tie(ait,aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
   {
      Vertex w = *ait;

      if ( fixed[w] == BACSunfixed )
      {
         assert( degrees[w] > 0 );
         assert( degrees[w] < boost::num_vertices(*G) );

         // possibly adjust counts
         if ( degrees[w] == 3 )
            ++ndegreetwo;
         else if ( degrees[w] == 2 )
         {
            assert( ndegreetwo > 0 );
            --ndegreetwo;
            ++ndegreeone;
         }
         else if ( degrees[w] == 1 )
         {
            assert( ndegreeone > 0 );
            --ndegreeone;
            ++ndegreezero;
         }

         // adjust degree of neighbor
         --degrees[w];
      }
   }

   // adjust counts w.r.t. v
   if ( degrees[v] == 2 )
   {
      assert( ndegreetwo > 0 );
      --ndegreetwo;
   }
   else if ( degrees[v] == 1 )
   {
      assert( ndegreeone > 0 );
      --ndegreeone;
   }
   else if ( degrees[v] == 0 )
   {
      assert( ndegreezero > 0 );
      --ndegreezero;
   }

   // finally fix node
   fixed[v] = BACSfixedzero;
   degrees[v] = SIZE_MAX;
}


/** fix node to 0 and update data - fast version without update of neighbors */
inline
void BACSfixNodeZeroFast(
   const Graph*          G,                  /**< graph to be processed */
   Vertex                v,                  /**< node to be fixed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees to be updated */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo          /**< number of degree 2 nodes */
   )
{
   assert( G != nullptr );
   assert( v < boost::num_vertices(*G) );
   assert( fixed != nullptr );
   assert( fixed[v] == BACSunfixed );
   assert( degrees != nullptr );

   if ( degrees[v] == 2 )
      --ndegreetwo;
   else if ( degrees[v] == 1 )
      --ndegreeone;
   else if ( degrees[v] == 0 )
      --ndegreezero;

   fixed[v] = BACSfixedzero;
   degrees[v] = SIZE_MAX;
}


/** fix node to 0 and update data - extended version */
inline
void BACSfixNodeZeroExtended(
   const Graph*          G,                  /**< graph to be processed */
   const Graph*          E,                  /**< additional edges to be processed */
   Vertex                v,                  /**< node to be fixed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees to be updated */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo          /**< number of degree 2 nodes */
   )
{
   assert( G != nullptr );
   assert( v < boost::num_vertices(*G) );
   assert( fixed != nullptr );
   assert( fixed[v] == BACSunfixed );
   assert( degrees != nullptr );

   // adjust degrees of neighboring nodes
   bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(v, G, E);
   for (; ! ait.at_end(); ++ait)
   {
      Vertex w = *ait;

      if ( fixed[w] == BACSunfixed )
      {
         assert( degrees[w] > 0 );

         // possibly adjust counts
         if ( degrees[w] == 3 )
            ++ndegreetwo;
         else if ( degrees[w] == 2 )
         {
            assert( ndegreetwo > 0 );
            --ndegreetwo;
            ++ndegreeone;
         }
         else if ( degrees[w] == 1 )
         {
            assert( ndegreeone > 0 );
            --ndegreeone;
            ++ndegreezero;
         }

         // adjust degree of neighbor
         --degrees[w];
      }
   }

   // adjust counts w.r.t. v
   if ( degrees[v] == 2 )
      --ndegreetwo;
   else if ( degrees[v] == 1 )
      --ndegreeone;
   else if ( degrees[v] == 0 )
      --ndegreezero;

   fixed[v] = BACSfixedzero;
   degrees[v] = SIZE_MAX;
}


/** fix node to 1 and update data */
inline
void BACSfixNodeOne(
   const Graph*          G,                  /**< graph to be processed */
   Vertex                v,                  /**< node to be fixed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees to be updated */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   size_t&               nfixedzero          /**< number of nodes fixed to 0 to be updated */
   )
{
   assert( G != nullptr );
   assert( v < boost::num_vertices(*G) );
   assert( fixed != nullptr );
   assert( fixed[v] == BACSunfixed );
   assert( degrees != nullptr );

   // first adjust counts w.r.t. v
   if ( degrees[v] == 2 )
      --ndegreetwo;
   else if ( degrees[v] == 1 )
      --ndegreeone;
   else if ( degrees[v] == 0 )
      --ndegreezero;

   // fix node
   fixed[v] = BACSfixedone;
   degrees[v] = SIZE_MAX;

   // fix adjacenty nodes
   AdjacencyIterator ait, aend;
   for (boost::tie(ait,aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
   {
      Vertex w = *ait;
      assert( fixed[w] != BACSfixedone );
      if ( fixed[w] == BACSunfixed )
      {
         BACSfixNodeZero(G, w, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
         ++nfixedzero;
      }
   }
}


/** fix node to 1 and update data - fast version without visiting neighbors */
inline
void BACSfixNodeOneFast(
   const Graph*          G,                  /**< graph to be processed */
   Vertex                v,                  /**< node to be fixed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees to be updated */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo          /**< number of degree 2 nodes */
   )
{
   assert( G != nullptr );
   assert( v < boost::num_vertices(*G) );
   assert( fixed != nullptr );
   assert( fixed[v] == BACSunfixed );
   assert( degrees != nullptr );

   if ( degrees[v] == 2 )
      --ndegreetwo;
   else if ( degrees[v] == 1 )
      --ndegreeone;
   else if ( degrees[v] == 0 )
      --ndegreezero;

   fixed[v] = BACSfixedone;
   degrees[v] = SIZE_MAX;
}


/** fix node to 1 and update data - extended version */
inline
void BACSfixNodeOneExtended(
   const Graph*          G,                  /**< graph to be processed */
   const Graph*          E,                  /**< additional edges to be processed */
   Vertex                v,                  /**< node to be fixed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees to be updated */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   size_t&               nfixedzero          /**< returns number of nodes fixed to 0 */
   )
{
   assert( G != nullptr );
   assert( v < boost::num_vertices(*G) );
   assert( fixed != nullptr );
   assert( fixed[v] == BACSunfixed );
   assert( degrees != nullptr );

   nfixedzero = 0;

   if ( degrees[v] == 2 )
      --ndegreetwo;
   else if ( degrees[v] == 1 )
      --ndegreeone;
   else if ( degrees[v] == 0 )
      --ndegreezero;

   fixed[v] = BACSfixedone;
   degrees[v] = SIZE_MAX;

   // fix adjacent nodes
   bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(v, G, E);
   for (; ! ait.at_end(); ++ait)
   {
      Vertex w = *ait;
      assert( fixed[w] != BACSfixedone );
      if ( fixed[w] == BACSunfixed )
      {
         BACSfixNodeZeroExtended(G, E, w, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
         ++nfixedzero;
      }
   }
}


/** check for two adjacent nodes whether the neighborhood of one node is contained in the other and possibly fix nodes */
inline
SCIP_Bool BACScheckNeighborhoodContainmentSorted(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees to be updated */
   SCIP_Real*            weights,            /**< new weights */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   Vertex                s,                  /**< first node of edge */
   Vertex                t                   /**< second node of edgoe */
   )
{
   assert( scip != nullptr );
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );
   assert( fixed[s] == BACSunfixed && fixed[t] == BACSunfixed );

   // loop over neighborhood of s node
   AdjacencyIterator sit, send;
   boost::tie(sit, send) = boost::adjacent_vertices(s, *G);

   // loop over neighborhood of t node
   AdjacencyIterator tit, tend;
   boost::tie(tit, tend) = boost::adjacent_vertices(t, *G);

   // setting up Booleans deciding whether neighborhood of s is subset of t and vice versa
   SCIP_Bool sint = TRUE;
   SCIP_Bool tins = TRUE;

   // checking objective criteria
   SCIP_Real sobj = weights[s];
   SCIP_Real tobj = weights[t];

   if ( SCIPisGT(scip, sobj, tobj) )
      tins = FALSE;
   if ( SCIPisGT(scip, tobj, sobj) )
      sint = FALSE;

   // loop over both adjacency lists to determine inclusion in either direction
   while ( sit != send || tit != tend )
   {
      // if we already have non-inclusion in both directions
      if ( ! sint && ! tins )
         break;

      // s iterator at end
      if ( sit == send )
      {
         if ( fixed[*tit] == BACSunfixed && *tit != s )
            tins = FALSE;
         ++tit;
         continue;
      }

      // t iterator at end
      if ( tit == tend )
      {
         if ( fixed[*sit] == BACSunfixed && *sit != t )
            sint = FALSE;
         ++sit;
         continue;
      }

      // advance iterators over shared neighborhood
      if ( *sit == *tit )
      {
         ++sit;
         ++tit;
         continue;
      }

      // s has neighbor which t does not have
      if ( *sit < *tit )
      {
         if ( fixed[*sit] == BACSunfixed && *sit != t )
            sint = FALSE;
         ++sit;
         continue;
      }

      // t has neighbor which s does not have
      if ( *sit > *tit )
      {
         if ( fixed[*tit] == BACSunfixed && *tit != s )
            tins = FALSE;
         ++tit;
         continue;
      }

      // previous cases are exhaustive
      abort();
   }

   // fix one variable if condition is satisfied
   if ( sint )
   {
      SCIPdebugMsg(scip, "Fixed node %lu to 0, because its neighborhood contains that of node %lu.\n", t, s);
      BACSfixNodeZero(G, t, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
      return TRUE;
   }
   else if ( tins )
   {
      SCIPdebugMsg(scip, "Fixed node %lu to 0, because its neighborhood contains that of node %lu.\n", s, t);
      BACSfixNodeZero(G, s, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
      return TRUE;
   }

   return FALSE;
}


/** check for two adjacent nodes whether the neighborhood of one node is contained in the other and possibly fix nodes
 *
 *  Slightly slower version that does not rely on the adjecency lists being sorted and takes additional edges in E into account.
 */
inline
SCIP_Bool BACScheckNeighborhoodContainmentExtended(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   const Graph*          E,                  /**< additional edges to be processed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   Vertex                s,                  /**< first node of edge */
   Vertex                t,                  /**< second node of edgoe */
   Vertex*               marked,             /**< auxiliary array for marking neighborhoods */
   size_t                nneighs             /**< number ob neighbors of s */
   )
{
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );
   assert( marked != nullptr );
   assert( fixed[s] == BACSunfixed && fixed[t] == BACSunfixed );
   assert( weights != nullptr );
   // get weights
   SCIP_Real sobj = weights[s];
   SCIP_Real tobj = weights[t];

   // loop over neighborhood of t node
   size_t ntins = 1;    // initialize with 1 to account for s
   size_t nneight = 0;
   bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(t, G, E);
   for (; ! ait.at_end(); ++ait)
   {
      if ( fixed[*ait] == BACSunfixed )
      {
         ++nneight;
         if ( marked[*ait] == s )
         {
            ++ntins;

            // if we already found all neighbors of s, we can stop
            if ( ntins == nneighs && SCIPisLE(scip, tobj, sobj) )
            {
               // make sure that we do not accidently run into the first case below
               nneight = SIZE_MAX;
               break;
            }
         }
      }
   }
   assert( ntins <= nneight );
   assert( ntins <= nneighs );

   // if all neighbors of t are in s
   if ( ntins == nneight && SCIPisLE(scip, sobj, tobj) )
   {
      SCIPdebugMsg(scip, "Fixed node %lu to 0, because its neighborhood contains that of node %lu.\n", s, t);
      BACSfixNodeZeroExtended(G, E, s, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
      return TRUE;
   }

   // if all neighbors of s are in t
   if ( ntins == nneighs && SCIPisLE(scip, tobj, sobj) )
   {
      SCIPdebugMsg(scip, "Fixed node %lu to 0, because its neighborhood contains that of node %lu.\n", t, s);
      BACSfixNodeZeroExtended(G, E, t, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
      return TRUE;
   }

   return FALSE;
}


/** neighborhood graph presolving - fast version relying on sorted adjacency lists in G */
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
   )
{
   assert( scip != nullptr );
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );

   nfixedzero = 0;

   // loop through all start vertices
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex s = *vit;

      // skip fixed nodes
      if ( fixed[s] != BACSunfixed )
         continue;

      // loop through adjacent nodes
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(s, *G); ait != aend; ++ait)
      {
         Vertex t = *ait;

         // skip fixed nodes
         if ( fixed[t] != BACSunfixed )
            continue;

         // test edge only once for t < s (we can break since the adjacent list is sorted)
         if ( t > s )
            break;

         assert( fixed[s] == BACSunfixed && fixed[t] == BACSunfixed );
         if ( BACScheckNeighborhoodContainmentSorted(scip, G, fixed, degrees, weights, ndegreezero, ndegreeone, ndegreetwo, s, t) )
         {
            assert( fixed[s] == BACSfixedzero || fixed[t] == BACSfixedzero );
            ++nfixedzero;

            // exit s-loop if s has been fixed
            if ( fixed[s] != BACSunfixed )
               break;
         }
      }
   }
}

/** neighborhood graph presolving - version that takes additional edges in E into account */
static
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
   )
{
   assert( scip != nullptr );
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );

   nfixedzero = 0;

   size_t n = boost::num_vertices(*G);
   assert( E == nullptr || n == boost::num_vertices(*E) );

   // init auxiliary array
   Vertex* marked;
   SCIP_CALL( SCIPallocBufferArray(scip, &marked, n) );
   for (size_t i = 0; i < n; ++i)
      marked[i] = Graph::null_vertex();

   // loop through all start vertices
   for (size_t i = 0; i < n; ++i)
   {
      Vertex s = (Vertex) i;

      // skip fixed nodes
      if ( fixed[s] != BACSunfixed )
         continue;

      // loop over neighborhood of s node
      size_t nneighs = 0;
      bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(s, G, E);
      for (; ! ait.at_end(); ++ait)
      {
         if ( fixed[*ait] == BACSunfixed )
         {
            marked[*ait] = s;
            ++nneighs;
         }
      }

      // loop through adjacent nodes
      ait = bacs::adjacent_vertices(s, G, E);
      for (; ! ait.at_end(); ++ait)
      {
         Vertex t = *ait;

         // skip fixed nodes
         if ( fixed[t] != BACSunfixed )
            continue;

         // test edge only once for s < t
         if ( t <= s )
            continue;

         assert( fixed[s] == BACSunfixed && fixed[t] == BACSunfixed );
         if ( BACScheckNeighborhoodContainmentExtended(scip, G, E, fixed, degrees, weights, ndegreezero, ndegreeone, ndegreetwo, s, t, marked, nneighs) )
         {
            assert( fixed[s] == BACSfixedzero || fixed[t] == BACSfixedzero );
            ++nfixedzero;

            // exit s-loop if s has been fixed
            if ( fixed[s] != BACSunfixed )
               break;
         }
      }
   }

   SCIPfreeBufferArray(scip, &marked);

   return SCIP_OKAY;
}

/** check for two non-adjacent nodes whether their neighborhoods are the same */
inline
SCIP_Bool BACScheckNeighborhoodEqualitySorted(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   Vertex                s,                  /**< first node of edge */
   Vertex                t                   /**< second node of edge */
   )
{
   assert( scip != nullptr );
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( fixed[s] == BACSunfixed && fixed[t] == BACSunfixed );

   // loop over neighborhood of s node
   AdjacencyIterator sit, send;
   boost::tie(sit, send) = boost::adjacent_vertices(s, *G);

   // loop over neighborhood of t node
   AdjacencyIterator tit, tend;
   boost::tie(tit, tend) = boost::adjacent_vertices(t, *G);

   // setting up Boolean deciding whether neighborhoods are the same
   SCIP_Bool seqt = TRUE;

   // loop over both adjacency lists
   while ( sit != send || tit != tend )
   {
      // s iterator at end
      if ( sit == send )
      {
         if ( fixed[*tit] == BACSunfixed )
         {
            seqt = FALSE;
            break;
         }
         ++tit;
         continue;
      }

      // t iterator at end
      if ( tit == tend )
      {
         if ( fixed[*sit] == BACSunfixed )
         {
            seqt = FALSE;
            break;
         }
         ++sit;
         continue;
      }

      // advance iterators over shared neighborhood
      if ( *sit == *tit )
      {
         ++sit;
         ++tit;
         continue;
      }

      // s has neighbor which t does not have
      if ( *sit < *tit  )
      {
         if ( fixed[*sit] == BACSunfixed )
         {
            seqt = FALSE;
            break;
         }
         ++sit;
         continue;
      }

      // t has neighbor which s does not have
      if ( *sit > *tit )
      {
         if ( fixed[*tit] == BACSunfixed )
         {
            seqt = FALSE;
            break;
         }
         ++tit;
         continue;
      }

      // previous cases are exhaustive
      abort();
   }
   return seqt;
}

/** check for two non-adjacent nodes whether their neighborhoods are the same and possibly merge nodes
 *
 *  Slightly slower version that does not rely on the adjecency lists being sorted and takes additional edges in E into account.
 */
inline
SCIP_Bool BACScheckNeighborhoodEqualityExtended(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   const Graph*          E,                  /**< additional edges to be processed */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   Vertex                s,                  /**< first node of edge */
   Vertex                t,                  /**< second node of edgoe */
   Vertex*               marked              /**< auxiliary array for marking neighborhoods */
   )
{
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( marked != nullptr );
   assert( fixed[s] == BACSunfixed && fixed[t] == BACSunfixed );

   // setting up Boolean deciding whether neighborhoods are the same
   SCIP_Bool seqt = TRUE;

   // loop over neighborhood of t node
   bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(t, G, E);
   for (; ! ait.at_end(); ++ait)
   {
      if ( fixed[*ait] != BACSunfixed )
         continue;

      if ( marked[*ait] != s )
      {
         seqt = FALSE;
         break;
      }

   }
   return seqt;
}

/** merge operation - fast version relying on sorted adjacency lists in G */
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
   )
{
   assert( scip != nullptr );
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );
   assert( mergemapping != nullptr );

   nmerged = 0;
   size_t n = boost::num_vertices(*G);

   // init auxiliary array
   Vertex* marked;
   SCIP_CALL( SCIPallocBufferArray(scip, &marked, n) );
   for (size_t i = 0; i < n; ++i)
      marked[i] = Graph::null_vertex();

   // loop through all vertices
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex s = *vit;

      // skip fixed nodes
      if ( fixed[s] != BACSunfixed )
         continue;

      // loop over neighborhood of s node
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(s, *G); ait != aend; ++ait)
      {
         if ( fixed[*ait] == BACSunfixed )
         {
            marked[*ait] = s;
         }
      }

      // loop through all vertices
      VertexIterator vit2 = vit;
      ++vit2;
      for (; vit2 != vend; ++vit2)
      {
         Vertex t = *vit2;

         // test nodes only once for s < t and skip already fixed nodes
         if ( fixed[t] != BACSunfixed )
            continue;

         // only consider nodes with the same degree
         if ( degrees[s] != degrees[t] )
            continue;

         // only consider non-adjacent nodes
         if ( marked[t] == s )
            continue;

         assert( fixed[s] == BACSunfixed && fixed[t] == BACSunfixed );

         // if nodes have same neighborhood contract them into one node with sum of objective values being new objective value
         if ( BACScheckNeighborhoodEqualitySorted(scip, G, fixed, s, t) )
         {
            BACSfixNodeZero(G, t, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
            weights[s] += weights[t];
            isweighted = TRUE;
            SCIPdisjointsetUnion(mergemapping, (int) s, (int) t, TRUE);
            SCIPdebugMsg(scip, "Merged node %lu into node %lu, because they have the same neighborhood. Node %lu now has objective value %f. \n", t, s, s, weights[s]);
            ++nmerged;
         }
      }
   }
   SCIPfreeBufferArray(scip, &marked);
   return SCIP_OKAY;
}

/** merge presolving - version that takes additional edges in E into account */
static
SCIP_RETCODE BACSgraphMergePresolvingExtended(
   SCIP*                 scip,               /**< SCIP pointer */
   const Graph*          G,                  /**< graph to be processed */
   const Graph*          E,                  /**< additional edges to be processed */
   SCIP_Bool&            isweighted,         /**< is graph weighted? */
   BACS_NODEFIXING*      fixed,              /**< array for fixings of nodes */
   size_t*               degrees,            /**< degrees */
   SCIP_Real*            weights,            /**< new weights */
   SCIP_DISJOINTSET*     mergemapping,       /**< disjoint set structure such that representatives are the nodes in the presolved graph */
   size_t&               ndegreezero,        /**< number of degree 0 nodes */
   size_t&               ndegreeone,         /**< number of degree 1 nodes */
   size_t&               ndegreetwo,         /**< number of degree 2 nodes */
   size_t&               nmerged             /**< return number of pairs of nodes have have been merged */
   )
{
   assert( scip != nullptr );
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );
   assert( mergemapping != nullptr );

   nmerged = 0;

   // init auxiliary array
   Vertex* marked;
   size_t n = boost::num_vertices(*G);
   assert( E == nullptr || n == boost::num_vertices(*E) );
   SCIP_CALL( SCIPallocBufferArray(scip, &marked, n) );
   for (size_t i = 0; i < n; ++i)
      marked[i] = Graph::null_vertex();

   // loop through all vertices
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex s = *vit;

      // skip fixed nodes
      if ( fixed[s] != BACSunfixed )
         continue;

      // loop over neighborhood of s node
      bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(s, G, E);
      for (; ! ait.at_end(); ++ait)
      {
         if ( fixed[*ait] == BACSunfixed )
            marked[*ait] = s;
      }

      // loop through all vertices
      VertexIterator vit2 = vit;
      ++vit2;
      for (; vit2 != vend; ++vit2)
      {
         Vertex t = *vit2;

         // test nodes only once for s < t and skip already fixed nodes
         if ( fixed[t] != BACSunfixed )
            continue;

         // only consider nodes with the same degree
         if ( degrees[s] != degrees[t] )
            continue;

         // only consider non-adjacent nodes
         if ( marked[t] == s )
            continue;

         assert( fixed[s] == BACSunfixed && fixed[t] == BACSunfixed );

         // if nodes have same neighborhood contract them into one node with sum of objective values being new objective value
         if ( BACScheckNeighborhoodEqualityExtended(scip, G, E, fixed, s, t, marked) )
         {
            BACSfixNodeZeroExtended(G, E, t, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
            SCIP_Real sobj = weights[s];
            SCIP_Real tobj = weights[t];
            weights[s] = sobj + tobj;
            isweighted = TRUE;
            SCIPdisjointsetUnion(mergemapping, (int) s, (int) t, TRUE);
            SCIPdebugMsg(scip, "Merged node %lu into node %lu, because they have the same neighborhood. Node %lu now has objective value %f. \n", t, s, s, weights[s]);
            ++nmerged;
         }
      }
   }

   SCIPfreeBufferArray(scip, &marked);

   return SCIP_OKAY;
}

/** perform simplicial presolving w.r.t. degrees 0/1 */
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
   )
{
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );

   nfixedzero = 0;
   nfixedone = 0;

   if ( ndegreezero + ndegreeone > 0 )
   {
      size_t n = boost::num_vertices(*G);

      for (size_t i = 0; i < n; ++i)
      {
         if ( degrees[i] == 0 )
         {
            assert( fixed[i] == BACSunfixed );
            BACSfixNodeOneFast(G, (Vertex) i, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
            ++nfixedone;
            assert( fixed[i] == BACSfixedone );
            assert( degrees[i] == SIZE_MAX );
         }
         else if ( degrees[i] == 1 )
         {
            assert( fixed[i] == BACSunfixed );

            Vertex s = (Vertex) i;
            SCIP_Real sobj = weights[s];

            // find neighbor
            bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(s, G, E);
            for (; ! ait.at_end(); ++ait)
            {
               if ( fixed[*ait] == BACSunfixed )
                  break;
            }
            assert( ! ait.at_end() );

            Vertex t = *ait;
            SCIP_Real tobj = weights[t];

            if ( SCIPisGE(scip, sobj, tobj) )
            {
               size_t nlocalfixedzero;
               BACSfixNodeOneExtended(G, E, s, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo, nlocalfixedzero);
               nfixedzero += nlocalfixedzero;
               ++nfixedone;
               assert( fixed[i] == BACSfixedone );
               assert( degrees[i] == SIZE_MAX );
            }
         }
      }
   }
}


/** perform cycle presolving */
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
   )
{
   assert( scip != nullptr );
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );

   nfixedzero = 0;
   nfixedone = 0;

   size_t n = boost::num_vertices(*G);
   assert( E == nullptr || n == boost::num_vertices(*E) );
   size_t m = boost::num_edges(*G);
   if ( E != nullptr )
      m += boost::num_edges(*E);

   // init arrays for DFS
   SCIP_Bool* visited;
   Vertex* stack;
   Vertex* parent;
   Vertex* markcycle;
   SCIP_CALL( SCIPallocBufferArray(scip, &visited, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stack, m) );
   SCIP_CALL( SCIPallocBufferArray(scip, &parent, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &markcycle, n) );

   // init visited and parent
   for (size_t i = 0; i < n; ++i)
   {
      visited[i] = FALSE;
      parent[i] = Graph::null_vertex();
      markcycle[i] = Graph::null_vertex();
   }

   // loop through all nodes
   for (size_t i = 0; i < n; ++i)
   {
      Vertex v = (Vertex) i;
      assert( weights[v] == 1.0 );

      // skip already visited nodes
      if ( visited[v] )
         continue;

      // skip fixed nodes
      if ( fixed[v] == BACSfixedzero || fixed[v] == BACSfixedone )
         continue;

      // mark current node as visited and put it on the stack
      stack[0] = v;
      size_t nstack = 1;

      // loop until stack is empty
      while ( nstack > 0 )
      {
         // pop node u from stack
         Vertex u = stack[--nstack];

         if ( visited[u] )
            continue;
         visited[u] = TRUE;

         // skip fixed nodes (note that some nodes might be fixed below while they are on the stack)
         if ( fixed[u] == BACSfixedzero || fixed[u] == BACSfixedone )
            continue;

         // the following assert does not need to hold, because previous fixings might have changed the degrees
         // assert( degrees[u] >= 2 );

         // check neighbors of u
         bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(u, G, E);
         for (; ! ait.at_end(); ++ait)
         {
            Vertex w = *ait;
            assert( w != u );

            // skip fixed nodes
            if ( fixed[w] == BACSfixedzero || fixed[w] == BACSfixedone )
               continue;

            // check whether new node is not visited
            if ( ! visited[w] )
            {
               // put new node on the stack and set parent
               parent[w] = u;
               stack[nstack++] = w;
               assert( nstack <= m );
            }
            // if we have found a cycle
            else if ( w != parent[u] )
            {
               // loop along the cycle to determine the number of nodes that can be fixed to 1
               size_t ncycle = 0;
               size_t nhighdeg = 0;
               size_t ncycleone = 0;  // number of nodes that can be fixed to 1
               size_t firsthighdeg = SIZE_MAX; // index of first node with degree > 2
               size_t lasthighdeg = SIZE_MAX;  // index of last node with degree > 2
               Vertex x = u;
               bool stop = false;
               while ( x != Graph::null_vertex() )
               {
                  if( fixed[x] != BACSunfixed )
                  {
                     stop = true;
                     break;
                  }

                  if ( degrees[x] > 2 )
                  {
                     // mark node to be of higher degree and in cycle
                     markcycle[x] = u;
                     if ( firsthighdeg == SIZE_MAX )
                     {
                        assert( lasthighdeg == SIZE_MAX );
                        firsthighdeg = ncycle;
                     }
                     else if ( lasthighdeg != SIZE_MAX )
                     {
                        assert( firsthighdeg < n );
                        assert( lasthighdeg < n );

                        /* The number of nodes that can be fixed to 1 between the current higher degree node and the
                         * last is the number of nodes between them plus 1 divided by 2. */
                        ncycleone += (ncycle - lasthighdeg) / 2;
                     }
                     lasthighdeg = ncycle;
                     ++nhighdeg;
                  }
                  ++ncycle;

                  // possibly stop early: if the higher degree nodes are fixed to 0, are there degree 2 nodes left?
                  if ( nhighdeg > (ncycle + 1)/ 2 )
                     break;

                  // stop if we are back to starting node
                  if ( x == w )
                     break;
                  x = parent[x];
               }
               if ( stop )
                  continue;

               assert( x != Graph::null_vertex() );

               // if we stopped early
               if ( nhighdeg > (ncycle + 1)/ 2 )
                  continue;

               /* Determine the number of nodes between the first higher degree node and the end of the cycle. */
               if ( firsthighdeg != SIZE_MAX )
               {
                  assert( firsthighdeg < n );
                  assert( lasthighdeg < n );
                  ncycleone += ((ncycle - lasthighdeg) + firsthighdeg) / 2;
               }
               else
                  ncycleone = ncycle / 2;

               // Check whether we were able to set at least ncycle / 2 rounded down many nodes to 1.
               if ( ncycleone >= ncycle / 2 )
               {
                  // determine the starting value
                  int parity = 1; // 1: fixing to 1, -1: fixing to 0
                  if ( firsthighdeg != SIZE_MAX )
                  {
                     assert( nhighdeg > 0 );
                     if ( firsthighdeg % 2 == 0 )
                        parity = -1;
                  }

                  // loop over cycle
                  x = u;
                  ncycle = 0;
                  ncycleone = 0;
                  size_t ncyclezero = 0;
                  while ( x != Graph::null_vertex() )
                  {
                     if ( markcycle[x] == u )
                     {
                        SCIPdebugMsg(scip, "fix node %lu of higher degree in cycle to 0.\n", x);
                        if ( fixed[x] == BACSunfixed )
                           BACSfixNodeZero(G, x, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
                        ++ncyclezero;
                        parity = 1; // reset parity
                     }
                     else
                     {
                        // if we want to fix to 1
                        if ( parity == 1 )
                        {
                           // w already fixed to 0 if u was fixed to 1
                           if( x == w && fixed[u] == BACSfixedone )
                           {
                              assert( fixed[x] == BACSfixedzero );
                              ++ncyclezero;
                              ++ncycle;
                              break;
                           }
                           SCIPdebugMsg(scip, "fix node %lu of degree 2 in cycle to 1.\n", x);
                           size_t k = 0;
                           BACSfixNodeOne(G, x, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo, k);
                           ++ncycleone;
                        }
                        else
                        {
                           assert( parity == -1 );
                           // otherwise fix to 0
                           SCIPdebugMsg(scip, "fix node %lu of degree 2 in cycle to 0.\n", x);
                           if ( fixed[x] == BACSunfixed )
                              BACSfixNodeZero(G, x, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
                           ++ncyclezero;
                        }
                        parity = -parity; // flip parity
                     }

                     markcycle[x] = Graph::null_vertex();
                     ++ncycle;
                     if ( x == w )
                        break;
                     x = parent[x];
                  }

                  assert( ncyclezero + ncycleone == ncycle );
                  assert( ncycleone == ncycle / 2 );
                  nfixedzero += ncyclezero;
                  nfixedone += ncycleone;

                  SCIPdebugMsg(scip, "Fixed cycle of length %lu (higher degree nodes: %lu).\n", ncycle, nhighdeg);
               }
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &markcycle);
   SCIPfreeBufferArray(scip, &visited);
   SCIPfreeBufferArray(scip, &stack);
   SCIPfreeBufferArray(scip, &parent);

   return SCIP_OKAY;
}


/** perform one complete round of presolving */
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
   )
{
   assert( scip != nullptr );
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );

   nfixedzeroneigh = 0;
   nfixedzerosimpl = 0;
   nfixedzerocycle = 0;
   nfixedonesimpl = 0;
   nfixedonecycle = 0;
   nmerged = 0;

   // first perform neighborhood presolving
   if ( neighborhoodpresol )
      BACSgraphNeighborhoodPresolving(scip, G, fixed, degrees, weights, ndegreezero, ndegreeone, ndegreetwo, nfixedzeroneigh);

   // then perform simplicial presolving
   if ( simplicialpresol )
      BACSsimplicialPresolving(scip, G, nullptr, fixed, degrees, weights, ndegreezero, ndegreeone, ndegreetwo, nfixedzerosimpl, nfixedonesimpl);

   // perform merge operation
   if ( mergepresol )
   {
      SCIP_CALL( BACSgraphMergePresolving(scip, G, isweighted, fixed, degrees, weights, mergemapping, ndegreezero, ndegreeone, ndegreetwo, nmerged) );
      if ( nmerged > 0 )
         cyclepresol = FALSE;
   }

   // run cycle presolving
   if ( cyclepresol && ndegreetwo > 0 )
   {
      SCIP_CALL( BACScyclePresolving(scip, G, nullptr, fixed, degrees, weights, ndegreezero, ndegreeone, ndegreetwo, nfixedzerocycle, nfixedonecycle) );
   }
#ifndef NDEBUG
   BACScheckConsistentFixings(G, nullptr, fixed);
   SCIP_CALL( BACScheckDegrees(scip, G, nullptr, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo) );
#endif

   return SCIP_OKAY;
}


/** perform one complete round of presolving, also using additional edges in E */
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
   )
{
   assert( scip != nullptr );
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );

   nfixedzeroneigh = 0;
   nfixedzerosimpl = 0;
   nfixedzerocycle = 0;
   nfixedonesimpl = 0;
   nfixedonecycle = 0;
   nmerged = 0;

   // first perform neighborhood presolving
   if ( neighborhoodpresol )
   {
      SCIP_CALL( BACSgraphNeighborhoodPresolvingExtended(scip, G, E, fixed, degrees, weights, ndegreezero, ndegreeone, ndegreetwo, nfixedzeroneigh) );
   }

   // then perform simplicial presolving
   if ( simplicialpresol )
      BACSsimplicialPresolving(scip, G, E, fixed, degrees, weights, ndegreezero, ndegreeone, ndegreetwo, nfixedzerosimpl, nfixedonesimpl);

   // perform merge operation
   if ( mergepresol )
   {
      SCIP_CALL( BACSgraphMergePresolvingExtended(scip, G, E, isweighted, fixed, degrees, weights, mergemapping, ndegreezero, ndegreeone, ndegreetwo, nmerged) );
      if ( nmerged > 0 )
         cyclepresol = FALSE;
   }

   // run cycle presolving
   if ( cyclepresol && ndegreetwo > 0 )
   {
      SCIP_CALL( BACScyclePresolving(scip, G, E, fixed, degrees, weights, ndegreezero, ndegreeone, ndegreetwo, nfixedzerocycle, nfixedonecycle) );
   }
#ifndef NDEBUG
   BACScheckConsistentFixings(G, E, fixed);
   SCIP_CALL( BACScheckDegrees(scip, G, E, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo) );
#endif

   return SCIP_OKAY;
}


/** perform rounds of presolving until no changes happen anymore, also using additional edges in E */
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
   )
{
   nfixedzeroneigh = 0;
   nfixedzerosimpl = 0;
   nfixedzerocycle = 0;
   nfixedonesimpl = 0;
   nfixedonecycle = 0;
   nmerged = 0;

   SCIP_Bool changed;
   do
   {
      ++nrounds;

      // perform one round of presolving
      size_t nroundfixedzeroneigh;
      size_t nroundfixedzerosimpl;
      size_t nroundfixedzerocycle;
      size_t nroundfixedonesimpl;
      size_t nroundfixedonecycle;
      size_t nroundmerged;

      SCIP_CALL( BACSpresolvingRoundExtend(scip, G, E, isweighted, neighborhoodpresol, simplicialpresol, mergepresol, cyclepresol,
            fixed, degrees, weights, mergemapping, ndegreezero, ndegreeone, ndegreetwo,
            nroundfixedzeroneigh, nroundfixedzerosimpl, nroundfixedzerocycle, nroundfixedonesimpl, nroundfixedonecycle, nroundmerged) );
      assert( neighborhoodpresol || nroundfixedzeroneigh == 0 );
      assert( simplicialpresol || nroundfixedzerosimpl + nroundfixedonesimpl == 0 );
      assert( cyclepresol || nroundfixedzerocycle + nroundfixedonecycle == 0 );

      SCIPinfoMessage(scip, nullptr, "Round %lu: neighborhood (z: %lu), merge (z: %lu), simplicial (z: %lu, o: %lu), cycle (z: %lu, o: %lu).\n",
         nrounds, nroundfixedzeroneigh, nroundmerged, nroundfixedzerosimpl, nroundfixedonesimpl, nroundfixedzerocycle, nroundfixedonecycle);

      if ( nroundmerged > 0 )
         cyclepresol = FALSE;

      nfixedzeroneigh += nroundfixedzeroneigh;
      nfixedzerosimpl += nroundfixedzerosimpl;
      nfixedzerocycle += nroundfixedzerocycle;
      nfixedonesimpl += nroundfixedonesimpl;
      nfixedonecycle += nroundfixedonecycle;
      nmerged += nroundmerged;

      if ( nroundfixedzeroneigh + nroundfixedzerosimpl + nroundfixedzerocycle + nroundfixedonesimpl + nroundfixedonecycle + nroundmerged > 0 )
         changed = TRUE;
      else
         changed = FALSE;
   }
   while ( changed );

   return SCIP_OKAY;
}


/** perform in-probing */
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
   )
{
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );
   assert( E != nullptr );

   naddededges = 0;

   size_t n = boost::num_vertices(*G);
   assert( n == boost::num_vertices(*E) );

   // local copy of fixed
   BACS_NODEFIXING* localfixed;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &localfixed, fixed, n) );
   assert( localfixed != nullptr );

   // local copy of degrees
   size_t* localdegrees;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &localdegrees, degrees, n) );
   assert( localdegrees != nullptr );
   size_t nlocaldegreezero = ndegreezero;
   size_t nlocaldegreeone = ndegreeone;
   size_t nlocaldegreetwo = ndegreetwo;

   // array to mark neighbors of current nodes
   Vertex* neighbors;
   SCIP_CALL( SCIPallocBufferArray(scip, &neighbors, n) );
   assert( neighbors != nullptr );
   for (size_t i = 0; i < n; ++i)
      neighbors[i] = Graph::null_vertex();

   // loop through all nodes
   VertexIterator vit, vend;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;

      if ( localfixed[v] != BACSunfixed )
         continue;

      // check timelimit
      if ( ! SCIPisInfinity(scip, timelimit) )
      {
         struct tms timer_end;
         (void) times(&timer_end);
         double t = (timer_end.tms_utime - timer_beg.tms_utime) / (double)sysconf(_SC_CLK_TCK);
         if ( t > timelimit )
            break;
      }

      SCIPdebugMsg(scip, "In-Probing on node %lu ...\n", v);
#ifndef NDEBUG
      for (size_t i = 0; i < n; ++i)
         assert( localfixed[i] == fixed[i] );
#endif

      // tentatively fix v to 1
      size_t nfixedzero;
      BACSfixNodeOneExtended(G, E, v, localfixed, localdegrees, nlocaldegreezero, nlocaldegreeone, nlocaldegreetwo, nfixedzero);

      // perform one round of presolving
      size_t nroundfixedzeroneigh;
      size_t nroundfixedzerosimpl;
      size_t nroundfixedzerocycle;
      size_t nroundfixedonesimpl;
      size_t nroundfixedonecycle;
      size_t nroundmerged;

      SCIP_CALL( BACSpresolvingRoundExtend(scip, G, E, isweighted, TRUE, TRUE, FALSE, FALSE, localfixed, localdegrees, weights, mergemapping,
            nlocaldegreezero, nlocaldegreeone, nlocaldegreetwo,
            nroundfixedzeroneigh, nroundfixedzerosimpl, nroundfixedzerocycle, nroundfixedonesimpl, nroundfixedonecycle, nroundmerged) );
      assert( nroundfixedzerocycle == 0 );
      assert( nroundfixedonecycle == 0 );

      // reset directly changed entries
      localfixed[v] = fixed[v];
      bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(v, G, E);
      for (; ! ait.at_end(); ++ait)
         localfixed[*ait] = fixed[*ait];

      // check for variables that have been fixed to 0
      if ( nroundfixedzeroneigh + nroundfixedzerosimpl > 0 )
      {
         // loop through nodes
         VertexIterator vit2, vend2;
         for (boost::tie(vit2, vend2) = boost::vertices(*G); vit2 != vend2; ++vit2)
         {
            Vertex w = *vit2;

            if ( localfixed[w] == BACSfixedzero && localfixed[w] != fixed[w] )
            {
               assert( fixed[w] == BACSunfixed );
               assert( localfixed[w] != BACSunfixed );

               // add non-existing edges to list
               if ( neighbors[w] != v )
               {
                  assert( ! boost::edge(v, w, *G).second );
                  assert( ! boost::edge(v, w, *E).second );

                  SCIPdebugMsg(scip, "In-Probing detected edge (%lu, %lu).\n", v, w);
#ifndef NDEBUG
                  std::pair<Edge, bool> p = boost::add_edge(v, w, *E);
                  assert( p.second );
#else
                  (void) boost::add_edge(v, w, *E);
#endif
                  ++naddededges;

                  // update degrees
                  assert( degrees[v] < SIZE_MAX );
                  assert( degrees[w] < SIZE_MAX );

                  if ( degrees[v] == 0 )
                  {
                     --ndegreezero;
                     ++ndegreeone;
                  }
                  else if ( degrees[v] == 1 )
                  {
                     --ndegreeone;
                     ++ndegreetwo;
                  }
                  else if ( degrees[v] == 2 )
                     --ndegreetwo;

                  if ( degrees[w] == 0 )
                  {
                     --ndegreezero;
                     ++ndegreeone;
                  }
                  else if ( degrees[w] == 1 )
                  {
                     --ndegreeone;
                     ++ndegreetwo;
                  }
                  else if ( degrees[w] == 2 )
                     --ndegreetwo;

                  ++degrees[v];
                  ++degrees[w];

#ifndef NDEBUG
                  SCIP_CALL( BACScheckDegrees(scip, G, E, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo) );
#endif
               }
            }
         }
      }

      // possibly reset fixed arrays
      if ( nroundfixedzeroneigh + nroundfixedzerosimpl + nroundfixedonesimpl > 0 )
      {
         BMScopyMemoryArray(localfixed, fixed, n);
      }
      // alsways copy degrees, since they have changed from neighbors to neighors of v
      BMScopyMemoryArray(localdegrees, degrees, n);
      nlocaldegreezero = ndegreezero;
      nlocaldegreeone = ndegreeone;
      nlocaldegreetwo = ndegreetwo;
#ifndef NDEBUG
      SCIP_CALL( BACScheckDegrees(scip, G, E, localfixed, localdegrees, nlocaldegreezero, nlocaldegreeone, nlocaldegreetwo) );
#endif
   }

   SCIPfreeBufferArray(scip, &neighbors);
   SCIPfreeBufferArray(scip, &localdegrees);
   SCIPfreeBufferArray(scip, &localfixed);

#ifndef NDEBUG
   SCIP_CALL( BACScheckDegrees(scip, G, E, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo) );
#endif

   return SCIP_OKAY;
}


/** SST presolving */
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
   size_t&               nmerged             /**< return number of pairs of nodes that have been merged */
   )
{
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );
   assert( weights != nullptr );
   assert( mergemapping != nullptr );
   assert( components != nullptr );

   nfixedzero = 0;
   naddededges = 0;
   nmerged = 0;

   size_t n = boost::num_vertices(*G);
   assert( n == boost::num_vertices(*E) );

   // variables for local timing
   struct tms timer_beg;
   struct tms timer_end;

   long double groupsize;
   VecVecUInt gens;
   (void) times(&timer_beg);
   SCIP_CALL( computeAutomorphismsGraph(scip, G, E, fixed, weights, FALSE, groupsize, gens) );
   (void) times(&timer_end);

   SCIPinfoMessage(scip, nullptr, "Number of generators of Aut(G):\t\t%lu\n", gens.size());
   SCIPinfoMessage(scip, nullptr, "Size of Aut(G):\t\t\t\t%.1Lg\n", groupsize);
   double t = (timer_end.tms_utime - timer_beg.tms_utime) / (double)sysconf(_SC_CLK_TCK);
   SCIPinfoMessage(scip, nullptr, "Computing symmetry time:\t\t%4.2f\n", t);

   // do nothing if there are no generators
   if ( gens.size() == 0 )
   {
      SCIPinfoMessage(scip, nullptr, "No symmetry left.\n");
      return SCIP_OKAY;
   }

   // initialize markers
   if ( ncomponents > 1 )
   {
      // sort components
      std::vector<size_t> sortedcomponents(n);
      for (size_t i = 0; i < n; ++i)
         sortedcomponents[i] = i;
      std::sort(sortedcomponents.begin(), sortedcomponents.end(), [&](int a, int b){ return components[a] < components[b]; });

      // store starting indices of each component of sortedcomponents
      std::vector<size_t> component_start_indices(ncomponents + 1);
      int oldcomp = -1;
      for (size_t i = 0; i < n; ++i)
      {
         int comp = components[sortedcomponents[i]];
         assert( comp >= oldcomp );
         if ( comp != oldcomp )
         {
            assert( comp >= 0 );
            component_start_indices[comp] = i; /*lint !e732*/
            oldcomp = comp;
         }
      }
      assert( components[sortedcomponents[n-1]] == (int) ncomponents -1 );
      component_start_indices[ncomponents] = n;

      // merge components
      size_t ncomponentsmerged = 0;
      for (size_t curgen = 0; curgen < gens.size(); ++curgen)
      {
         VertexIterator vit, vend;
         for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
         {
            Vertex v = *vit;
            if ( components[v] < 0 || fixed[v] != BACSunfixed )
               continue;

            Vertex imgv = gens[curgen][v];
            assert( components[imgv] >= 0 );
            assert( components[v] >= 0 );
            if ( components[imgv] != components[v] )
            {
               // skip if components are already merged or single nodes are fixed
               if ( fixed[imgv] != BACSunfixed || components[imgv] == components[SCIPdisjointsetFind(mergemapping, (int) v)] )
                  continue;

               // merge component components[imgv] onto component components[v]
               assert( component_start_indices[components[v]] <= component_start_indices[components[v] + 1] ); /*lint !e732*/
               for (size_t i = component_start_indices[components[v]]; i < component_start_indices[components[v] + 1]; ++i) /*lint !e732*/
               {
                  Vertex w = (Vertex) sortedcomponents[i];
                  assert( components[w] == components[v] );

                  Vertex imgw = gens[curgen][w];
                  assert( components[imgw] == components[imgv] );
                  if ( fixed[imgw] != BACSunfixed )
                     continue;

                  // get representing node
                  Vertex origin = (Vertex) (unsigned) SCIPdisjointsetFind(mergemapping, (int) w);

                  // fix image to 0 and add weights to represting node
                  BACSfixNodeZeroExtended(G, E, imgw, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
                  weights[origin] += weights[imgw];
                  isweighted = TRUE;
                  SCIPdisjointsetUnion(mergemapping, (int) origin, (int) imgw, TRUE);
                  ++nmerged;
               }
               SCIPdebugMsg(scip, "Merged component %d into component %d, because they are symmetric.\n", components[imgv], components[v]);
               ++ncomponentsmerged;
            }

            // stop if all components are merged
            if ( ncomponentsmerged == ncomponents )
               break;
         }
      }
      SCIPinfoMessage(scip, nullptr, "Number of components merged:\t\t%lu\n", ncomponentsmerged);
   }

   // compute leaders and followers
   VecVertex leaders;
   VecVecVertex followers;

   SCIP_CALL( computeLeaderFollowersSSTgensFilter(scip, G, n, gens, orbitrule, fixed, leaders, followers, degrees) );

   SCIPinfoMessage(scip, nullptr, "Number of leaders:\t\t\t%zu\n", leaders.size());
   size_t nfollowers = 0;
   for (unsigned int i = 0; i < followers.size(); ++i)
      nfollowers += followers[i].size();
   SCIPinfoMessage(scip, nullptr, "Number of followers:\t\t\t%lu\n", nfollowers);

   if ( leaders.size() > 0 )
   {
      // collect followers that are adjacent to their leader
      std::vector<Vertex> marked(n, Graph::null_vertex());
      std::vector<size_t> nleaderremoved(leaders.size(), 0);
      for (size_t lit = 0; lit < leaders.size(); ++lit)
      {
         Vertex l = leaders[lit];
         assert( fixed[l] != BACSfixedone );

         if ( fixed[l] == BACSfixedzero )
            continue;

         // mark followers
         for (size_t fit = 0; fit < followers[lit].size(); ++fit)
         {
            Vertex f = followers[lit][fit];
            assert( fixed[f] != BACSfixedone );

            if ( fixed[f] != BACSfixedzero )
               marked[f] = l;
         }

         // collect removed nodes: nodes adjacent to l, but not removed earlier
         bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(l, G, E);
         for (; ! ait.at_end(); ++ait)
         {
            Vertex u = *ait;
            if ( marked[u] == l && fixed[u] != BACSfixedzero )
            {
               assert( fixed[u] != BACSfixedone );
               BACSfixNodeZeroExtended(G, E, u, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);
               ++nfixedzero;
               ++nleaderremoved[lit];
            }
         }
      }

      // add edges
      if ( sstaddedges )
      {
         // initialize markers
         VertexIterator vit, vend;
         for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
            marked[*vit] = Graph::null_vertex();

         // if {l, f} is not an edge, add edge {u, f} for all neighbors u of l
         for (size_t lit = 0; lit < leaders.size(); ++lit)
         {
            Vertex l = leaders[lit];
            assert( fixed[l] != BACSfixedone );

            if ( fixed[l] != BACSunfixed )
               continue;

            // skip if all followers have been removed
            if ( followers[lit].size() == nleaderremoved[lit] )
               continue;

            // loop through neighbors u of l
            bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(l, G, E);
            for (; ! ait.at_end(); ++ait)
            {
               Vertex u = *ait;
               assert( fixed[u] != BACSfixedone );

               // skip already removed nodes
               if ( fixed[u] != BACSunfixed )
                  continue;

               // mark neighbors of u
               bacs::ExtAdjacencyIterator ait2 = bacs::adjacent_vertices(u, G, E);
               for (; ! ait2.at_end(); ++ait2)
                  marked[*ait2] = u;

               // loop through followers
               for (size_t fit = 0; fit < followers[lit].size(); ++fit)
               {
                  Vertex f = followers[lit][fit];
                  assert( fixed[f] != BACSfixedone );

                  // skip already removed nodes
                  if ( fixed[f] != BACSunfixed )
                     continue;

                  // add edges to followers that are not yet adjacent
                  if ( marked[f] != u )
                  {
                     assert( ! boost::edge(u, f, *G).second );
                     assert( ! boost::edge(u, f, *E).second );

                     (void) boost::add_edge(u, f, *E);
                     ++naddededges;

                     assert( degrees[u] < SIZE_MAX );
                     assert( degrees[f] < SIZE_MAX );

                     if ( degrees[u] == 0 )
                     {
                        --ndegreezero;
                        ++ndegreeone;
                     }
                     else if ( degrees[u] == 1 )
                     {
                        --ndegreeone;
                        ++ndegreetwo;
                     }
                     else if ( degrees[u] == 2 )
                        --ndegreetwo;

                     if ( degrees[f] == 0 )
                     {
                        --ndegreezero;
                        ++ndegreeone;
                     }
                     else if ( degrees[f] == 1 )
                     {
                        --ndegreeone;
                        ++ndegreetwo;
                     }
                     else if ( degrees[f] == 2 )
                        --ndegreetwo;

                     ++degrees[u];
                     ++degrees[f];
                  }
               }
            }
         }
      }
   }

#ifndef NDEBUG
   SCIP_CALL( BACScheckDegrees(scip, G, E, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo) );
#endif

   return SCIP_OKAY;
}
