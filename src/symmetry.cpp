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

/**@file   symmetry.cpp
 * @brief  Methods to deal with symmetries
 * @author Marc Pfetsch
 *
 * This code is based on the implementation in maxkcol by Christopher Hojny.
 */

#include "symmetry.h"
#include "extiterators.hpp"
#include <stack>


/* include bliss */
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4189)  // local variable is initialized but not referenced
# pragma warning(disable: 4267)  // conversion of size_t to int (at sassy/preprocessor.h:2897)
# pragma warning(disable: 4388)  // compare signed and unsigned expression
# pragma warning(disable: 4456)  // shadowed variable
# pragma warning(disable: 4430)  // missing type specifier
#endif


#ifdef BACS_WITH_NAUTY
extern "C" {
#include <nauty/nauty.h>
#include <nauty/nausparse.h>
}
#endif

#ifdef BACS_WITH_BLISS
/* include bliss */
#include <bliss/defs.hh>
#include <bliss/graph.hh>
#endif

#ifdef BACS_WITH_SBLISS
/* include bliss */
#include <bliss/defs.hh>
#include <bliss/graph.hh>
/* include sassy */
#include <sassy/preprocessor.h>
#include <sassy/tools/bliss_converter.h>
#endif

#ifdef __GNUC__
#pragma GCC diagnostic warning "-Wzero-as-null-pointer-constant"
#pragma GCC diagnostic warning "-Wunused-but-set-variable"
#pragma GCC diagnostic warning "-Wsign-compare"
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wshadow"
#endif

#ifdef _MSC_VER
# pragma warning(pop)
#endif


//! Select an orbit and its leader from a list of orbits
SCIP_RETCODE selectOrbitAndLeader(
   const Graph*          G,                  //!< graph
   const VecVecVertex&   orbits,             //!< vector of orbits (each orbit has size > 1)
   Vertex&               leader,             //!< return selected leader
   size_t&               orbitidx,           //!< return selected orbit
   ORBIT_RULE            orbitrule,          //!< rule to select orbit
   size_t*               degrees,            //!< degrees
   SCIP_Bool*            wsviolcliques       //!< cliques that violate the WS property
   )
{
   leader = Graph::null_vertex();
   orbitidx = SIZE_MAX;

   // select an orbit and its leader
   if ( orbitrule == ORBIT_RULE_MAXORBIT )
   {
      orbitidx = 0;
      leader = orbits[0][0];
      size_t orbitsize = orbits[orbitidx].size();

      for (size_t i = 1; i < orbits.size(); ++i)
      {
         // skip nodes that violate the WS property if this is the global orbitrule
         if ( wsviolcliques != nullptr )
         {
            if ( wsviolcliques[i] )
               continue;

            if ( wsviolcliques[orbitidx] )
            {
               orbitsize = orbits[i].size();
               orbitidx = i;
               leader = orbits[i][0];
            }
         }

         if ( orbits[i].size() > orbitsize )
         {
            orbitsize = orbits[i].size();
            orbitidx = i;
            leader = orbits[i][0];
         }
      }
   }
   else if ( orbitrule == ORBIT_RULE_MINORBIT )
   {
      orbitidx = 0;
      leader = orbits[0][0];
      size_t orbitsize = orbits[orbitidx].size();

      for (size_t i = 1; i < orbits.size(); ++i)
      {
         assert( orbits[i].size() >= 2 );
         if ( orbits[i].size() < orbitsize )
         {
            orbitsize = orbits[i].size();
            orbitidx = i;
            leader = orbits[i][0];
         }
      }
   }
   else if ( orbitrule == ORBIT_RULE_DEGREE )
   {
      orbitidx = 0;
      leader = orbits[0][0];
      size_t leaderdegree = degrees[orbits[0][0]];
      SCIP_Bool maxorbsizetwo = TRUE;
      SCIP_Bool stop = FALSE;

      for (size_t i = 1; i < orbits.size(); ++i)
      {
         if ( orbits[i].size() > 2 )
         {
            maxorbsizetwo = FALSE;
            break;
         }
      }

      for (size_t i = 1; i < orbits.size(); ++i)
      {
         if ( maxorbsizetwo ) // checking the existence of edges is quite expensive and should be made more efficiently if we decide to use this orbitrule
         {
            if ( orbits[i].size() == 2 && ! boost::edge(orbits[orbitidx][0], orbits[orbitidx][1], *G).second )
            {
               orbitidx = i;
               leader = orbits[i][1];
               if ( boost::edge(orbits[orbitidx][0], orbits[orbitidx][1], *G).second )
               {
                  stop = TRUE;
                  break;
               }
            }
         }
         else
         {
            for (size_t j =1; j < orbits[i].size(); ++j)
            {
               if ( degrees[orbits[i][j]] > leaderdegree )
               {
                  leaderdegree = degrees[orbits[i][j]];
                  orbitidx = i;
                  leader = orbits[i][j];
               }
            }
         }
         if ( stop )
            break;
      }
   }
   // choose the leader with the most adjacent followers (this does not take into account that some of these followers could already be adjacent followers of an earlier leader)
   else if ( orbitrule == ORBIT_RULE_MAXDEL )
   {
      orbitidx = 0;
      leader = orbits[0][0];

      int potdels = -1;
      std::vector<SCIP_Bool> orbitmarker(boost::num_vertices(*G), FALSE);

      // we only need to consider one element of each orbit since all orbit elements should have the same number of adjacent nodes in the orbit
      for (size_t i = 0; i < orbits.size(); ++i)
      {
         int potdelslocal = 0;
         Vertex v = orbits[i][0];

         // skip nodes that violate the WS property if this is the global orbitrule
         if ( wsviolcliques != nullptr )
         {
            if ( wsviolcliques[i] )
               continue;
         }

         if ( (int) orbits[i].size() - 1 <= potdels || (int) degrees[v] <= potdels )
            continue;

         // mark elements in orbit
         for (size_t j = 1; j < orbits[i].size(); ++j)
            orbitmarker[orbits[i][j]] = TRUE;

         AdjacencyIterator ait, aend;
         for (boost::tie(ait,aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
         {
            if ( ! orbitmarker[*ait] )
               continue;

            ++potdelslocal;
         }
         if ( potdelslocal > potdels )
         {
            potdels = potdelslocal;
            orbitidx = i;
            leader = v;
         }
         std::fill(orbitmarker.begin(), orbitmarker.end(), FALSE);
      }
   }
   else
   {
      assert( orbitrule == ORBIT_RULE_FIRST || orbitrule == ORBIT_RULE_STRINGENT ||orbitrule == ORBIT_RULE_STRINGENT_DEGREE || orbitrule == ORBIT_RULE_STRINGENT_MAXDEL || orbitrule == ORBIT_RULE_WSTRINGENT || orbitrule == ORBIT_RULE_WSTRINGENT_MAXDEL );

      // select the first variable of the first orbit as leader
      leader = orbits[0][0];
      orbitidx = 0;

      for (size_t i = 1; i < orbits[0].size(); ++i)
      {
         if ( orbits[0][i] < leader )
            leader = orbits[0][i];
      }
   }

   return SCIP_OKAY;
}

//! Select a stringent orbit and its leader from a list of orbits
static
SCIP_RETCODE selectOrbitAndLeaderStringent(
   const Graph*          G,                  //!< graph
   const VecVecVertex&   orbits,             //!< vector of orbits (each orbit has size > 1)
   Vertex&               leader,             //!< return selected leader
   size_t&               orbitidx,           //!< return selected orbit
   std::list<Vertex>&    preferrednodes,     //!< LIFO list of nodes for which orbits are computed preferably
   std::vector<int>      node2orbit,         //!< assigns each node its orbitidx (or -1)
   ORBIT_RULE            orbitrule,          //!< rule to select orbit
   size_t*               degrees             //!< degrees
   )
{
   leader = Graph::null_vertex();
   orbitidx = SIZE_MAX;
   size_t osize = 0;

   std::list<Vertex>::iterator pit = preferrednodes.begin();
   std::list<Vertex>::iterator pend = preferrednodes.end();
   std::list<Vertex>::iterator p = preferrednodes.begin();

   Vertex nextnode;
   while ( pit != pend )
   {
      nextnode = *pit;

      // if the node is not affected by symmetry anymore
      int idx = node2orbit[nextnode];
      if ( idx == -1 )
      {
         (void)preferrednodes.erase(pit++);
         continue;
      }

      assert( idx >= 0 );
      // choose largest orbit (also for ORBIT_RULE_DEGREE since all nodes in preferrednodes have the same degree)
      if ( orbits[(size_t)idx].size() > osize )
      {
         p = pit;
         leader = nextnode;
         orbitidx = (size_t) idx;
         osize = orbits[orbitidx].size();
      }
      ++pit;
   }

   if ( leader != Graph::null_vertex() )
   {
      assert( orbitidx != SIZE_MAX );
      assert( *p == leader );

      // remove leader and add remaining elements of orbit
#ifndef NDEBUG
      // all elements should already in preferredvars
      for (unsigned int j = 0; j < orbits[orbitidx].size(); ++j)
         assert( std::find(preferrednodes.begin(), preferrednodes.end(), orbits[orbitidx][j]) != preferrednodes.end() );
#endif
      (void)preferrednodes.erase(p);

      return SCIP_OKAY;
   }

   // otherwise, add the first orbit
   assert( preferrednodes.size() == 0 );

   if ( orbitrule == ORBIT_RULE_STRINGENT )
   {
      SCIP_CALL( selectOrbitAndLeader(G, orbits, leader, orbitidx, ORBIT_RULE_MAXORBIT, degrees, nullptr) );
   }
   else if ( orbitrule == ORBIT_RULE_STRINGENT_DEGREE )
   {
      SCIP_CALL( selectOrbitAndLeader(G, orbits, leader, orbitidx, ORBIT_RULE_DEGREE, degrees, nullptr) );
   }
   else
   {
      assert( orbitrule == ORBIT_RULE_STRINGENT_MAXDEL );
      SCIP_CALL( selectOrbitAndLeader(G, orbits, leader, orbitidx, ORBIT_RULE_MAXDEL, degrees, nullptr) );
   }

   assert( leader != Graph::null_vertex() );
   assert( orbitidx != SIZE_MAX );

   // add followers to preferrednodes
   VecVertex::const_reverse_iterator oit, oend = orbits[orbitidx].rend();
   for (oit = orbits[orbitidx].rbegin(); oit != oend; ++oit)
   {
      if ( *oit == leader )
         continue;
      preferrednodes.push_back(*oit);
   }

   return SCIP_OKAY;
}

// Select a weakly stringent orbit and its leader from a list of orbits
static
SCIP_RETCODE selectOrbitAndLeaderWStringent(
   const Graph*          G,                  //!< graph
   const VecVecVertex&   orbits,             //!< vector of orbits (each orbit has size > 1)
   Vertex&               leader,             //!< return selected leader
   size_t&               orbitidx,           //!< return selected orbit
   std::vector<int>&     orbitnodes,         //!< list of nodes in orbits of leaders and leadercounter of the leader that caused their addition to the set
   const std::vector<int>& node2orbit,       //!< assigns each node its orbitidx (or -1)
   ORBIT_RULE            orbitrule,          //!< rule to select orbit
   size_t*               degrees,            //!< degrees
   SCIP_Bool*            wsviolcliques,      //!< cliques that violate the WS property
   size_t                leadercounter       //!< number of leaders already selected
   )
{
   leader = Graph::null_vertex();
   orbitidx = SIZE_MAX;

   // need to mark the nodes that violate the WS property
   SCIP_Bool foundwsleader = FALSE;
   SCIP_Bool foundviolation;
   while ( ! foundwsleader )
   {
      if ( orbitrule == ORBIT_RULE_WSTRINGENT_MAXDEL )
         SCIP_CALL( selectOrbitAndLeader(G, orbits, leader, orbitidx, ORBIT_RULE_MAXDEL, degrees, wsviolcliques) );
      else if ( orbitrule == ORBIT_RULE_WSTRINGENT )
         SCIP_CALL( selectOrbitAndLeader(G, orbits, leader, orbitidx, ORBIT_RULE_MAXORBIT, degrees, wsviolcliques) );
      assert( wsviolcliques[orbitidx] == FALSE );

      // if it is the first leader it is ws
      if ( leadercounter == 0 )
         break;

      // check whether the chosen leader is valid for the weakly stringent orbitrule
      VecVertex::const_iterator it, itend = orbits[orbitidx].end();
      for (it = orbits[orbitidx].begin(); it != itend; ++it)
      {
         Vertex u = *it;
         foundviolation = FALSE;
         AdjacencyIterator ait, aend;
         for (boost::tie(ait,aend) = boost::adjacent_vertices(u, *G); ait != aend; ++ait)
         {
            Vertex v = *ait;
            if ( orbitnodes[v] == -1 )
               continue;

            // check whether v is in orbit of size 1; then we can view it as fixed and can ignore it
            if ( node2orbit[v] == -1 )
               continue;

            // check whether u and v were in the same orbit at some point; this means they were added to orbitnodes in the same run
            if ( orbitnodes[u] != orbitnodes[v] )
            {
               foundviolation = TRUE;
               break;
            }
         }

         if ( foundviolation )
         {
            foundwsleader = FALSE;
            wsviolcliques[orbitidx] = TRUE;
            break;
         }
      }
      if ( foundviolation == FALSE )
         foundwsleader = TRUE;
   }

   assert( leader != Graph::null_vertex() );
   assert( orbitidx != SIZE_MAX );

   // add followers to orbitnodes
   VecVertex::const_iterator oit, oitend = orbits[orbitidx].end();
   for (oit = orbits[orbitidx].begin(); oit != oitend; ++oit)
   {
      if ( orbitnodes[*oit] == -1 )
         orbitnodes[*oit] = (int) leadercounter;
   }

   return SCIP_OKAY;
}

//! Compute leaders and followers for SST cuts heuristically by filtering generators
/*lint --e{715}*/
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
   )
{
   assert( G != nullptr );
   assert( degrees != nullptr );

   // initialize data structures to store which generators are still active (based on filtering)
   size_t ngens = gens.size();
   size_t nactivegens = ngens;
   std::vector<SCIP_Bool> genIsActive(ngens, TRUE);

   // to keep track which elements have already been added to an orbit in an iteration for a specific stabilizer
   std::vector<SCIP_Bool> elemadded(nelems, FALSE);

   // for stringent SST cuts, we need to keep track of preferred variables for which we compute the orbits
   std::list<Vertex> preferrednodes;

   // for weakly stringent SST cuts, we need to keep track of variables in orbits of leaders and when they have been added to this set
   std::vector<int> orbitnodes;
   if ( orbitrule == ORBIT_RULE_WSTRINGENT_MAXDEL || orbitrule == ORBIT_RULE_WSTRINGENT )
   {
      orbitnodes.resize(boost::num_vertices(*G));
      std::fill(orbitnodes.begin(), orbitnodes.end(), -1);
   }

   // iteratively compute orbits as long as filtered stabilizer is non-trivial
   SCIP_Bool printorbsize = TRUE;
   size_t leadercounter = 0;
   while ( nactivegens > 0 )
   {
      // collect the orbits of the current stabilizer
      VecVecVertex orbits;
      std::vector<int> node2orbit(nelems, -1);

      // compute orbits for each element of the ground set
      for (size_t v = 0; v < nelems; ++v)
      {
         VecVertex orbit;

         // don't compute orbit of already treated elements
         if ( elemadded[v] || fixed[v] != BACSunfixed )
            continue;

         orbit.push_back(v);
         elemadded[v] = TRUE;

         // iterate over all elements in the partially built orbit and compute their images
         size_t curpos = 0;
         while ( curpos < orbit.size() )
         {
            Vertex curelem;

            curelem = orbit[curpos];

            // compute orbits using generators
            for (size_t i = 0; i < ngens && nactivegens > 0; ++i)
            {
               if ( ! genIsActive[i] )
                  continue;

               Vertex image = (Vertex) gens[i][curelem];

               // if we have found a new element of the orbit
               if ( ! elemadded[image] && fixed[image] == BACSunfixed )
               {
                  orbit.push_back(image);
                  elemadded[image] = TRUE;
               }
            }
            ++curpos;
         }

         // add the orbit if it is non-trivial
         if ( orbit.size() > 1 )
         {
            orbits.push_back(orbit);
            for (VecVertex::iterator it = orbit.begin(); it != orbit.end(); ++it)
               node2orbit[*it] = (int) orbits.size() - 1;
         }
      }

      if ( orbits.size() == 0 )
         return SCIP_OKAY;

      if ( printorbsize )
      {
         // compute average orbit size and percentage of variables on which symmetry group is active
         size_t orbitsum = 0;
         for (size_t i = 0; i < orbits.size(); ++i)
            orbitsum += orbits[i].size();
         SCIP_Real avorbsize = (orbits.size() > 0) ? orbitsum / (SCIP_Real) orbits.size() : 0.0;
         SCIP_Real symvarfrac = (boost::num_vertices(*G) > 0) ? orbitsum / (SCIP_Real) boost::num_vertices(*G) : 0.0;
         SCIPinfoMessage(scip, nullptr, "Average orbit size:\t\t\t%.1f\n", avorbsize);
         SCIPinfoMessage(scip, nullptr, "Fraction of variables with symmetry:\t%.1f\n", symvarfrac);
         printorbsize = FALSE;
      }

      // select the orbit and its leader
      Vertex curleader;
      size_t curorbit;

      if ( orbitrule == ORBIT_RULE_STRINGENT || orbitrule == ORBIT_RULE_STRINGENT_DEGREE || orbitrule == ORBIT_RULE_STRINGENT_MAXDEL )
      {
         SCIP_CALL( selectOrbitAndLeaderStringent(G, orbits, curleader, curorbit, preferrednodes, node2orbit, orbitrule, degrees) );
      }
      else if ( orbitrule == ORBIT_RULE_WSTRINGENT_MAXDEL || orbitrule == ORBIT_RULE_WSTRINGENT )
      {
         SCIP_Bool* wsviolcliques;
         SCIP_CALL( SCIPallocClearBufferArray(scip, &wsviolcliques, orbits.size()) );
         SCIP_CALL( selectOrbitAndLeaderWStringent(G, orbits, curleader, curorbit, orbitnodes, node2orbit, orbitrule, degrees, wsviolcliques, leadercounter) );
         SCIPfreeBufferArray(scip, &wsviolcliques);
      }
      else
      {
         SCIP_CALL( selectOrbitAndLeader(G, orbits, curleader, curorbit, orbitrule, degrees, nullptr) );
      }

      VecVertex curfollowers;
      for (size_t i = 0; i < orbits[curorbit].size(); ++i)
      {
         if ( orbits[curorbit][i] != curleader )
            curfollowers.push_back(orbits[curorbit][i]);
      }

      // store orbit and leader
      followers.push_back(curfollowers);
      leaders.push_back(curleader);
      ++leadercounter;

      // clean auxiliary data structure
      for (size_t i = 0; i < nelems; ++i)
         elemadded[i] = FALSE;

      // filter generators
      for (size_t i = 0; i < ngens && nactivegens > 0; ++i)
      {
         if ( ! genIsActive[i] )
            continue;

         if ( gens[i][curleader] != curleader )
         {
            genIsActive[i] = FALSE;
            --nactivegens;
         }
      }
   }

   return SCIP_OKAY;
}


/* To define the callback function for nauty/traces we need the following global variables to store the information. */
static std::list<std::vector<unsigned int> >   generator_list_;        //!< generators of the graph automorphism group

#ifdef BACS_WITH_NAUTY
/** auxiliary function that nauty calls when it finds a generator
 *
 *  If a permutation is found, we put it into the global list of generators.
 */
static
void writeGenerator(
   int,
   int*                  perm,               //!< Generator that nauty found
   int*,
   int,
   int,
   int                   nNodes              //!< Number of nodes in the graph
   )
{
#ifndef NDEBUG
   for (int i = 0; i < nNodes; ++i)
   {
      assert( perm[i] >= 0 );
   }
#endif
   // use iterator constructor: copy permutation given at perm (up to position n-1)
   (void)generator_list_.emplace_back(std::vector<unsigned int>(perm, perm + nNodes) );
}

//! Compute generators of automorphism group of a graph with nauty
SCIP_RETCODE computeAutomorphismsGraphNauty(
   SCIP*                 scip,               //!< SCIP instance
   const Graph*          G,                  //!< graph
   const Graph*          E,                  //!< additional graph
   const BACS_NODEFIXING* fixed,             //!< array to mark fixed nodes
   const SCIP_Real*      weights,            //!< current weights
   long double&          groupsize,          //!< group size
   VecVecUInt&           generators          //!< returns generators
   )
{
   assert( G != nullptr );
   assert( E != nullptr );

   size_t n = boost::num_vertices(*G);

   assert( n > 0 );

   // init options
   const char* mallocString = "malloc";
   static DEFAULTOPTIONS_SPARSEGRAPH(options);
   statsblk stats;

   // init callback functions for nauty (accumulate the group generators found by nauty)
   options.writeautoms = FALSE;
   options.userautomproc = writeGenerator;
   // options.userlevelproc = grouplevelproc;
   options.defaultptn = FALSE; /* use color classes */

   // memory allocation for nauty
   DYNALLSTAT(int, lab, lab_sz);
   DYNALLSTAT(int, ptn, ptn_sz);
   DYNALLSTAT(int, orbits, orbits_sz);

   sparsegraph sg;
   SG_INIT(sg);

   // compute degrees of nodes
   size_t* degrees;
   SCIP_CALL( SCIPallocClearBufferArray(scip, &degrees, n) );
   SCIP_Real* symweights;
   SCIP_CALL( SCIPallocBufferArray(scip, &symweights, n) );


   VertexIterator vit, vend;
   size_t nedges = 0;
   size_t nfixed = 0;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      if ( fixed[v] == BACSunfixed )
      {
         symweights[v] = weights[v];
         assert( symweights[v] > 0.0 );

         bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(v, G, E);
         for (; ! ait.at_end(); ++ait)
         {
            Vertex w = *ait;
            if ( fixed[w] == BACSunfixed )
            {
               ++degrees[v];
               ++nedges;
            }
         }
      }
      else
         symweights[v] = - (SCIP_Real) (1 + nfixed++);  // mark fixed nodes with distinct negative numbers
   }
   assert( nedges <= 2 * (boost::num_edges(*G) + boost::num_edges(*E)) );

   // init space
   SG_ALLOC(sg, n, nedges, mallocString); /*lint !e774*/
   nauty_check(WORDSIZE, (int) nedges, (int) n, NAUTYVERSIONID);

   sg.nv = (int) n;
   sg.nde = nedges;

   DYNALLOC1(int, lab, lab_sz, n, mallocString);   /*lint !e727*/
   DYNALLOC1(int, ptn, ptn_sz, n, mallocString);   /*lint !e727*/

   // fill in array with colors for variables
   for (size_t i = 0; i < n; ++i)
      lab[i] = (int) i;

   // sort nodes according to weights
   SCIPsortRealInt(symweights, lab, (int) n);

   // fill in degrees and color classes
   std::vector<int> pos(n);
   size_t cnt = 0;
   for (size_t i = 0; i < n; ++i)
   {
      sg.d[i] = (int) degrees[i];   // degree of node i
      sg.v[i] = cnt;                // position of edges for node i
      pos[i] = (int) cnt;           // also store position
      cnt += degrees[i];

      if ( i < n - 1 && SCIPisEQ(scip, symweights[i], symweights[i+1]) )
         ptn[i] = 1;  // color class does not end
      else
         ptn[i] = 0;  // color class ends
   }

   SCIPfreeBufferArray(scip, &degrees);
   SCIPfreeBufferArray(scip, &symweights);

   // fill in graph
   for (boost::tie(vit,vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;

      if ( fixed[v] == BACSunfixed )
      {
         // create edges in nauty graph
         int p = pos[v];

         bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(v, G, E);
         for (; ! ait.at_end(); ++ait)
         {
            Vertex w = *ait;
            if ( fixed[w] == BACSunfixed )
               sg.e[p++] = (int) w;
         }
         assert( v == n - 1 || p == (int) sg.v[v+1] );
      }
   }

   DYNALLOC1(int, orbits, orbits_sz, n, mallocString); /*lint !e727*/

   // call nauty. Calls auxiliary method writeGenerator(Traces).
   sparsenauty(&sg, lab, ptn, orbits, &options, &stats, nullptr);

   // store number of generators
   assert( stats.numgenerators >= 0 );

   // put generator_list into generators and transposition generators into transpositions
   std::list<std::vector<unsigned int> >::const_iterator lit, lend = generator_list_.end();
   for (lit = generator_list_.begin(); lit != lend; ++lit)
   {
#ifndef NDEBUG
      // whether the output generator contains all numbers 0, ..., n-1
      std::vector<bool> containsZeroToN(n, false);
      for (size_t i = 0; i < lit->size(); ++i)
         containsZeroToN[ (*lit)[i] ] = true;
      for (unsigned int i = 0; i < n; ++i)
         assert( containsZeroToN[i] );
#endif

      // store permutation
      generators.push_back(*lit);
   }

   SG_FREE(sg); /*lint !e774*/

   // too many graph symmetry generators ...
   if ( generators.size() > n * n )
   {
      SCIPwarningMessage(scip, "The number of generators of Aut(G) is large: %lu.\n\n", generators.size());
   }

   // get the size of the group (stats.grapsize1 is a double)
   groupsize = stats.grpsize1 * exp10l(stats.grpsize2);

   return SCIP_OKAY;
}
#endif

#ifdef BACS_WITH_BLISS
/** callback function for bliss */
static
void blisshook(
   void*                 user_param,         /**< parameter supplied at call to bliss */
   unsigned int          dim,                /**< size of aut vector */
   const unsigned int*   aut                 /**< automorphism */
   )
{
   assert( aut != nullptr );
   (void)generator_list_.emplace_back(std::vector<unsigned int>(aut, aut + dim) );
}

//! Compute generators of automorphism group of a graph with bliss
static
SCIP_RETCODE computeAutomorphismsGraphBliss(
   SCIP*                 scip,               //!< SCIP instance
   const Graph*          G,                  //!< graph
   const Graph*          E,                  //!< additional graph
   const BACS_NODEFIXING* fixed,             //!< array to mark fixed nodes
   const SCIP_Real*      weights,            //!< current weights
   long double&          groupsize,          //!< group size
   VecVecUInt&           generators          //!< returns generators
   )
{
   assert( G != nullptr );
   assert( fixed != nullptr );

   size_t n = boost::num_vertices(*G);

   // compute degrees of nodes
   int* idx;
   SCIP_CALL( SCIPallocBufferArray(scip, &idx, n) );
   SCIP_Real* symweights;
   SCIP_CALL( SCIPallocBufferArray(scip, &symweights, n) );

   // collect weights
   VertexIterator vit, vend;
   size_t nfixed = 0;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      idx[v] = (int) v;  // init for sorting below (needs to be int)

      if ( fixed[v] == BACSunfixed )
      {
         symweights[v] = weights[v];
         assert( symweights[v] > 0.0 );
      }
      else
      {
         symweights[v] = -1.0;
         ++nfixed;
      }
   }

   // sort weights
   SCIPsortRealInt(symweights, idx, (int) n);

   // determine colors
   size_t* colors;
   SCIP_CALL( SCIPallocBufferArray(scip, &colors, n) );
   size_t color = 0;
   SCIP_Real lastweight = symweights[0];
   size_t cnt = 0;
   // colors of fixed nodes come first with distinct values, then all other nodes with values depending on distinct weights
   for (size_t i = 0; i < n; ++i)
   {
      assert( idx[i] >= 0 );
      size_t k = (size_t) idx[i];
      assert( k < n );
      if ( fixed[k] == BACSunfixed )
      {
         assert( symweights[i] >= 0.0 );
         if ( ! SCIPisEQ(scip, symweights[i], lastweight) )
            ++color;
         lastweight = symweights[i];
         colors[k] = color + nfixed; // shift colors
      }
      else
      {
         assert( symweights[i] < 0.0 );
         colors[k] = cnt++;
      }
   }
   assert( cnt == nfixed );

   // create bliss graph
   bliss::Graph blissgraph(0);

   // create nodes
   for (size_t i = 0; i < n; ++i)
   {
      (void) blissgraph.add_vertex(colors[i]);
   }

   // fill in edges
   bacs::ExtEdgeIterator eit = bacs::ExtEdgeIterator(G, E);
   for (; ! eit.at_end(); ++eit)
   {
      Vertex s = eit.source();
      Vertex t = eit.target();

      if ( fixed[s] == BACSunfixed && fixed[t] == BACSunfixed )
      {
         blissgraph.add_edge(s, t);
      }
   }

   SCIPfreeBufferArray(scip, &colors);
   SCIPfreeBufferArray(scip, &idx);
   SCIPfreeBufferArray(scip, &symweights);


   /* compute automorphisms */
   bliss::Stats stats;

   /* Prefer splitting partition cells corresponding to variables over those corresponding
    * to inequalities. This is because we are only interested in the action
    * of the automorphism group on the variables, and we need a base for this action */
   blissgraph.set_splitting_heuristic(bliss::Graph::shs_f);
   /* disable component recursion as advised by Tommi Junttila from bliss */
   blissgraph.set_component_recursion(false);

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   /* lambda function to have access to data and pass it to the blisshook above */
   auto reportglue = [&](unsigned int dim, const unsigned int* aut) {
      blisshook(nullptr, dim, aut);
   };

   /* start search */
   blissgraph.find_automorphisms(stats, reportglue, nullptr);
#else

   /* start search */
   blissgraph.find_automorphisms(stats, blisshook, nullptr);
#endif

#ifdef SCIP_OUTPUT
   (void) stats.print(stdout);
#endif

   // determine group size
   groupsize = stats.get_group_size_approx();

   // put generator_list into generators
   std::list<std::vector<unsigned int> >::const_iterator lit, lend = generator_list_.end();
   for (lit = generator_list_.begin(); lit != lend; ++lit)
   {
#ifndef NDEBUG
      // whether the output generator contains all numbers 0, ..., n-1
      std::vector<bool> containsZeroToN(n, false);
      for (size_t i = 0; i < lit->size(); ++i)
         containsZeroToN[ (*lit)[i] ] = true;
      for (unsigned int i = 0; i < n; ++i)
         assert( containsZeroToN[i] );
#endif

      // store permutation
      generators.push_back(*lit);
   }

   return SCIP_OKAY;
}
#endif


#ifdef BACS_WITH_SBLISS
/** callback function for sassy */
static
void sassyhook(
   int                   dim,                /**< size of aut vector */
   const int*            aut                 /**< automorphism */
   )
{
   assert( aut != nullptr );
   generator_list_.emplace_back(std::vector<unsigned int>(aut, aut + dim) );
}

//! Compute generators of automorphism group of a graph with bliss and preprocessing by sassy
static
SCIP_RETCODE computeAutomorphismsGraphSBliss(
   SCIP*                 scip,               //!< SCIP instance
   const Graph*          G,                  //!< graph
   const Graph*          E,                  //!< additional graph
   const BACS_NODEFIXING* fixed,             //!< array to mark fixed nodes
   const SCIP_Real*      weights,            //!< current weights
   long double&          groupsize,          //!< group size
   VecVecUInt&           generators          //!< returns generators
   )
{
   assert( G != nullptr );
   assert( fixed != nullptr );

   size_t n = boost::num_vertices(*G);

   // compute degrees of nodes
   int* idx;
   SCIP_CALL( SCIPallocBufferArray(scip, &idx, n) );
   SCIP_Real* symweights;
   SCIP_CALL( SCIPallocBufferArray(scip, &symweights, n) );

   // collect weights
   VertexIterator vit, vend;
   size_t nfixed = 0;
   for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
   {
      Vertex v = *vit;
      idx[v] = (int) v;  // init for sorting below (needs to be int)

      if ( fixed[v] == BACSunfixed )
      {
         symweights[v] = weights[v];
         assert( symweights[v] > 0.0 );
      }
      else
      {
         symweights[v] = -1.0;
         ++nfixed;
      }
   }

   // sort weights
   SCIPsortRealInt(symweights, idx, (int) n);

   // determine colors
   size_t* colors;
   size_t* degrees;
   SCIP_CALL( SCIPallocBufferArray(scip, &colors, n) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &degrees, n) );
   size_t color = 0;
   SCIP_Real lastweight = symweights[0];
   size_t cnt = 0;
   size_t nedges = 0;
   // colors of fixed nodes come first with distinct values, then all other nodes with values depending on distinct weights
   for (size_t i = 0; i < n; ++i)
   {
      assert( idx[i] >= 0 );
      Vertex v = (Vertex) idx[i]; /*lint !e571*/
      assert( v < n );
      if ( fixed[v] == BACSunfixed )
      {
         assert( symweights[i] > 0.0 );
         if ( ! SCIPisEQ(scip, symweights[i], lastweight) )
            ++color;
         lastweight = symweights[i];
         colors[v] = color + nfixed; // shift colors

         bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(v, G, E);
         for (; ! ait.at_end(); ++ait)
         {
            Vertex w = *ait;
            if ( fixed[w] == BACSunfixed )
            {
               ++degrees[v];
               if ( v < w )
                  ++nedges;
            }
         }
      }
      else
      {
         assert( symweights[i] <= 0.0 );
         colors[v] = cnt++;
      }
   }
   assert( cnt == nfixed );

   // create dejavu graph
   sassy::static_graph sassygraph;

   sassygraph.initialize_graph((unsigned int) n, (unsigned int) nedges);

   // create nodes
   for (size_t i = 0; i < n; ++i)
   {
      (void) sassygraph.add_vertex((int) colors[i], (int) degrees[i]);
   }

   // fill in edges
   bacs::ExtEdgeIterator eit = bacs::ExtEdgeIterator(G, E);
   for (; ! eit.at_end(); ++eit)
   {
      Vertex s = eit.source();
      Vertex t = eit.target();

      if ( fixed[s] == BACSunfixed && fixed[t] == BACSunfixed )
      {
         if ( s < t )
            sassygraph.add_edge((unsigned int) s, (unsigned int) t);
         else
            sassygraph.add_edge((unsigned int) t, (unsigned int) s);
      }
   }

   SCIPfreeBufferArray(scip, &degrees);
   SCIPfreeBufferArray(scip, &colors);
   SCIPfreeBufferArray(scip, &idx);
   SCIPfreeBufferArray(scip, &symweights);

   /* set up sassy preprocessor */
   sassy::preprocessor sassy;

   /* turn off some preprocessing that generates redudant permutations */
   sassy::configstruct sconfig;
   sconfig.CONFIG_PREP_DEACT_PROBE = true;
   sassy.configure(&sconfig);

   /* lambda function to have access to data and pass it to sassyhook above */
   sassy::sassy_hook sassyglue = [&](int dim, const int* p, int nsupp, const int* suppa) /*lint --e{715}*/ {
      sassyhook(dim, p);
   };

   /* call sassy to reduce graph */
   sassy.reduce(&sassygraph, &sassyglue);

   /* create bliss graph */
   bliss::Graph blissgraph(0);

   /* convert sassy to bliss graph */
   convert_sassy_to_bliss(&sassygraph, &blissgraph);

   /* compute automorphisms */
   bliss::Stats stats;

   /* Prefer splitting partition cells corresponding to variables over those corresponding
    * to inequalities. This is because we are only interested in the action
    * of the automorphism group on the variables, and we need a base for this action */
   blissgraph.set_splitting_heuristic(bliss::Graph::shs_f);
   /* disable component recursion as advised by Tommi Junttila from bliss */
   blissgraph.set_component_recursion(false);

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   auto hook = [&](unsigned int dim, const unsigned int* aut) {
      sassy.bliss_hook(dim, aut);
   };

   /* start search */
   blissgraph.find_automorphisms(stats, hook, nullptr);
#else

   /* start search */
   blissgraph.find_automorphisms(stats, sassy::preprocessor::bliss_hook, (void*) &sassy);
#endif

   // determine group size
   groupsize = sassy.base * powl(10.0, sassy.exp) * stats.get_group_size_approx();

   // put generator_list into generators
   std::list<std::vector<unsigned int> >::const_iterator lit, lend = generator_list_.end();
   for (lit = generator_list_.begin(); lit != lend; ++lit)
   {
#ifndef NDEBUG
      // whether the output generator contains all numbers 0, ..., n-1
      std::vector<bool> containsZeroToN(n, false);
      for (size_t i = 0; i < lit->size(); ++i)
         containsZeroToN[ (*lit)[i] ] = true;
      for (unsigned int i = 0; i < n; ++i)
         assert( containsZeroToN[i] );
#endif

      // store permutation
      generators.push_back(*lit);
   }

   return SCIP_OKAY;
}
#endif


#ifndef NDEBUG
//! check whether generators are actually automorphisms
void checkAutomorphismsGraph(
   SCIP*                 scip,               //!< SCIP instance
   const Graph*          G,                  //!< graph
   const Graph*          E,                  //!< additional graph
   const BACS_NODEFIXING* fixed,             //!< array to mark fixed nodes
   const VecVecUInt&     generators,         //!< returns generators
   SCIP_Bool             silent              //!< whether no output should be produced
   )
{
   assert( G != nullptr );
   assert( fixed != nullptr );

   size_t n = boost::num_vertices(*G);
   for (size_t i = 0; i < generators.size(); ++i)
   {
      const std::vector<unsigned int>& gen = generators[i];
      assert( gen.size() == n );

      // loop over generator to check that fixed nodes are invariant
      for (size_t j = 0; j < n; ++j)
      {
         unsigned int img = gen[j];
         if ( fixed[j] != BACSunfixed )
            assert( j == img );
      }

      // loop over edges
      bacs::ExtEdgeIterator eit = bacs::ExtEdgeIterator(G, E);
      for (; ! eit.at_end(); ++eit)
      {
         Vertex s = eit.source();
         Vertex t = eit.target();

         Vertex simg = (Vertex) gen[s];
         Vertex timg = (Vertex) gen[t];
         assert( simg < n && timg < n );

         // check that permuted edge exists
         if ( fixed[s] != BACSunfixed && fixed[t] != BACSunfixed )
         {
            // check is only necessary if s or t are changed
            if ( s != simg || t != timg )
               assert( boost::edge(simg, timg, *G).second );
         }
      }
   }
   if ( ! silent )
      SCIPinfoMessage(scip, nullptr, "Symmetry check passed.\n");
}
#endif


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
   )
{
   generators.clear();
   generator_list_.clear();
   groupsize = 0.0;

#ifdef BACS_WITH_NAUTY
   SCIPinfoMessage(scip, nullptr, "Using nauty to compute symmetries ...\n");
   SCIP_CALL( computeAutomorphismsGraphNauty(scip, G, E, fixed, weights, groupsize, generators) );
#endif

#ifdef BACS_WITH_BLISS
   SCIPinfoMessage(scip, nullptr, "Using bliss to compute symmetries ...\n");
   SCIP_CALL( computeAutomorphismsGraphBliss(scip, G, E, fixed, weights, groupsize, generators) );
#endif

#ifdef BACS_WITH_SBLISS
   SCIPinfoMessage(scip, nullptr, "Using sassy + bliss to compute symmetries ...\n");
   SCIP_CALL( computeAutomorphismsGraphSBliss(scip, G, E, fixed, weights, groupsize, generators) );
#endif

#ifndef NDEBUG
   checkAutomorphismsGraph(scip, G, E, fixed, generators, silent);
#endif

   return SCIP_OKAY;
}
