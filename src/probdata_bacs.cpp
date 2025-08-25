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

/**@file   probdata_bacs.cpp
 * @brief  handling of problem data
 * @author Marc Pfetsch
 */

#include "probdata_bacs.h"
#include "cons_clique.h"
#include "graphpresolve.h"
#include "extiterators.hpp"
#include "type_symmetry.h"

#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#include <boost/graph/connected_components.hpp>
#pragma GCC diagnostic warning "-Wzero-as-null-pointer-constant"
#pragma GCC diagnostic warning "-Wshadow"

#include <unistd.h>
#include <sys/times.h>

 // default values for parameters
#define DEFAULT_PRINTPROBSTATS     FALSE     //!< Print problem statistics?
#define DEFAULT_GRAPHPRESOLVING    FALSE     //!< Perform graph presolving?
#define DEFAULT_GRAPHINPROBING     FALSE     //!< Perform graph in-probing?
#define DEFAULT_MERGEPRESOL        TRUE      //!< Perform merge presolving?
#define DEFAULT_NEIGHBORHOODPRESOL TRUE      //!< Perform graph neighborhood presolving?
#define DEFAULT_SIMPLICIALPRESOL   TRUE      //!< Perform graph simplicial presolving?
#define DEFAULT_CYCLEPRESOL        TRUE      //!< Perform graph cycle presolving?
#define DEFAULT_DEGTWOCONTRACT     TRUE      //!< Perform contraction of degree 2 paths?
#define DEFAULT_SSTPRESOLVING      FALSE     //!< Perfrom SST presolving?
#define DEFAULT_SSTADDEDGES        FALSE     //!< Add edges during SST presolving?
#define DEFAULT_SSTREPEAT          TRUE      //!< Repeat SST presolving if successful
#define DEFAULT_SSTORBITRULE  ORBIT_RULE_STRINGENT   //!< orbits selection (0: maximum length, 1: minimum length, 2: stringent, 3: first variable, 4: maximum degree, 5: stringent degree, 6: max adjacent followers, 7: weakly stringent max variable fix, 8: stringent max variable fix, 9: weakly stringent)
#define DEFAULT_ONLYPRESOL         FALSE     //!< Only perform graph presolving?


//--------------------------------------------------------------------------------------------
//--------------------------------- event handler --------------------------------------------
//--------------------------------------------------------------------------------------------

#define EVENTHDLR_NAME         "localdegrees"
#define EVENTHDLR_DESC         "report bound changes for updating the local degrees"

//! event execution
SCIP_DECL_EVENTEXEC(BACSeventExecLocalDegrees)
{
   assert( scip != nullptr );
   assert( eventhdlr != nullptr );
   assert( event != nullptr );
   assert( eventdata != nullptr );

   SCIP_EVENTTYPE eventtype = SCIPeventGetType(event);
   assert( eventtype == SCIP_EVENTTYPE_LBTIGHTENED || eventtype == SCIP_EVENTTYPE_LBRELAXED || eventtype == SCIP_EVENTTYPE_UBTIGHTENED || eventtype == SCIP_EVENTTYPE_UBRELAXED );

   SCIP_VAR* var = SCIPeventGetVar(event);
   assert( var != nullptr );
   assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );
   assert( ! SCIPvarIsNegated(var) );

   SCIPdebugMsg(scip, "Event for changed bound of variable <%s> with new bound %f, type = %" SCIP_EVENTTYPE_FORMAT ".\n", SCIPvarGetName(var), SCIPeventGetNewbound(event), eventtype);

   // determine the bound that would correspond to a newly fixed variable
   SCIP_Real fixedbnd;
   if ( eventtype == SCIP_EVENTTYPE_LBTIGHTENED || eventtype == SCIP_EVENTTYPE_LBRELAXED )
      fixedbnd = 1.0;
   else
      fixedbnd = 0.0;

   // take care of negated variables
   if ( SCIPvarIsNegated(var) )
   {
      var = SCIPvarGetNegatedVar(var);
      assert( var != nullptr );
      // flip bound
      fixedbnd = 1.0 - fixedbnd;
   }
   assert( ! SCIPvarIsNegated(var) );

   SCIP_PROBDATA* probdata = (SCIP_PROBDATA*) eventdata;
   assert( probdata != nullptr );
   assert( probdata->eventhdlr != nullptr );
   assert( probdata->localdegrees != nullptr );
   const Graph* G = probdata->G;
   Vertex v = (Vertex) SCIPvarGetData(var);
   assert( SCIPvarGetTransVar(probdata->vars[v]) == var );

   if ( SCIPisEQ(scip, SCIPeventGetNewbound(event), fixedbnd) )
   {
      assert( SCIPvarGetUbLocal(var) < 0.5 || SCIPvarGetLbLocal(var) > 0.5 );

      // if variable is fixed to 0 or 1, reduce local degree of unfixed neighboring nodes
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;
         if ( probdata->localdegrees[w] < 0 )
            continue;

         --(probdata->localdegrees[w]);
         assert( probdata->localdegrees[w] >= 0 );
         assert( probdata->nlocaledges > 0 );
         --(probdata->nlocaledges);
      }
      probdata->localdegrees[v] = -1;

      // update number of local fixings
      if ( fixedbnd == 1 )
         ++(probdata->nlocalones);
      else
         ++(probdata->nlocalzeros);

      assert( probdata->nlocalones + probdata->nlocalzeros <= probdata->n );
   }
   else
   {
      // if variable is unfixed, increase local degree of unfixed neighboring nodes
      assert( probdata->localdegrees[v] == -1 );
      probdata->localdegrees[v] = 0;
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(v, *G); ait != aend; ++ait)
      {
         Vertex w = *ait;
         if ( probdata->localdegrees[w] < 0 )
            continue;

         ++(probdata->localdegrees[w]);
         assert( probdata->localdegrees[w] > 0 );

         ++(probdata->nlocaledges);
         assert( probdata->nlocaledges <= probdata->m );

         // also increase degree of unfixed variable
         ++probdata->localdegrees[v];
      }

      // update number of local one-fixings
      if ( fixedbnd == 1 )
      {
         assert( probdata->nlocalones > 0 );
         --(probdata->nlocalones);
      }
      else
      {
         assert( probdata->nlocalzeros > 0 );
         --(probdata->nlocalzeros);
      }
   }
   assert( probdata->n >= probdata->nlocalones + probdata->nlocalzeros );
   size_t nfreenodes = probdata->n - probdata->nlocalones - probdata->nlocalzeros;

   if ( nfreenodes <= 1 )
   {
      // if there is at most one vertex left, there should be no edges
      assert( probdata->nlocaledges == 0 );
      probdata->localdensity = 1.0;
   }
   else
      probdata->localdensity = (double) probdata->nlocaledges / (nfreenodes * (nfreenodes - 1) / 2);

   assert( SCIPisGE(scip, probdata->localdensity, 0.0) );
   assert( SCIPisLE(scip, probdata->localdensity, 1.0) );

   return SCIP_OKAY;
}



//! init problem data (allocate memory, init variables)
SCIP_RETCODE BACSinitProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< problem data structure */
   )
{
   SCIP_CALL( SCIPallocBlockMemory(scip, probdata) );
   (*probdata)->Gorig = nullptr;
   (*probdata)->G = nullptr;
   (*probdata)->changedpresol = FALSE;
   (*probdata)->origfixings = nullptr;
   (*probdata)->origtopresol = nullptr;
   (*probdata)->origtopresolparity = nullptr;
   (*probdata)->mergemapping = nullptr;
   (*probdata)->objoffset = 0.0;
   (*probdata)->norig = 0;
   (*probdata)->morig = 0;
   (*probdata)->n = 0;
   (*probdata)->m = 0;
   (*probdata)->unweighted = TRUE;
   (*probdata)->isolated = nullptr;
   (*probdata)->CG = nullptr;
   (*probdata)->vars = nullptr;

   SCIP_CALL( SCIPaddBoolParam(scip, "bacs/printprobstats", "Print problem statistics?",
         &(*probdata)->printprobstats, TRUE, DEFAULT_PRINTPROBSTATS, nullptr, nullptr) );
   SCIP_CALL( SCIPaddBoolParam(scip, "bacs/graphpresolving", "Perform graph presolving?",
         &(*probdata)->graphpresolving, TRUE, DEFAULT_GRAPHPRESOLVING, nullptr, nullptr) );
   SCIP_CALL( SCIPaddBoolParam(scip, "bacs/graphinprobing", "Perform graph in-probing?",
         &(*probdata)->graphinprobing, TRUE, DEFAULT_GRAPHINPROBING, nullptr, nullptr) );
   SCIP_CALL( SCIPaddBoolParam(scip, "bacs/mergepresol", "Perform merge presolving?",
      &(*probdata)->mergepresol, TRUE, DEFAULT_MERGEPRESOL, nullptr, nullptr) );
   SCIP_CALL( SCIPaddBoolParam(scip, "bacs/neighborhoodpresol", "Perform graph neighborhood presolving?",
         &(*probdata)->neighborhoodpresol, TRUE, DEFAULT_NEIGHBORHOODPRESOL, nullptr, nullptr) );
   SCIP_CALL( SCIPaddBoolParam(scip, "bacs/simplicialpresol", "Perform graph simplicial presolving?",
         &(*probdata)->simplicialpresol, TRUE, DEFAULT_SIMPLICIALPRESOL, nullptr, nullptr) );
   SCIP_CALL( SCIPaddBoolParam(scip, "bacs/cyclepresol", "Perform graph cycle presolving?",
         &(*probdata)->cyclepresol, TRUE, DEFAULT_CYCLEPRESOL, nullptr, nullptr) );
   SCIP_CALL( SCIPaddBoolParam(scip, "bacs/degtwocontract", "Perform contraction of degree 2 paths?",
         &(*probdata)->degtwocontract, TRUE, DEFAULT_DEGTWOCONTRACT, nullptr, nullptr) );
   SCIP_CALL( SCIPaddBoolParam(scip, "bacs/sstpresolving", "Perform SST presolving?",
         &(*probdata)->sstpresolving, TRUE, DEFAULT_SSTPRESOLVING, nullptr, nullptr) );
   SCIP_CALL( SCIPaddBoolParam(scip, "bacs/sstaddedges", "Add edges during SST presolving?",
         &(*probdata)->sstaddedges, TRUE, DEFAULT_SSTADDEDGES, nullptr, nullptr) );
   SCIP_CALL( SCIPaddBoolParam(scip, "bacs/sstrepeat", "Repeat SST presolving if successful?",
         &(*probdata)->sstrepeat, TRUE, DEFAULT_SSTREPEAT, nullptr, nullptr) );
   SCIP_CALL( SCIPaddIntParam(scip, "bacs/sstorbitrule", "SST orbit selection (0: maximum length, 1: minimum length, 2: stringent, 3: first variable, 4: maximum degree, 5: stringent degree, 6: max adjacent followers, 7: weakly stringent max variable fix, 8: stringent max variable fix, 9: weakly stringent)",
         &(*probdata)->sstorbitrule, TRUE, DEFAULT_SSTORBITRULE, 0, 9, nullptr, nullptr) );
   SCIP_CALL( SCIPaddBoolParam(scip, "bacs/onlypresol", "Only perform graph presolving?",
      &(*probdata)->onlypresol, TRUE, DEFAULT_ONLYPRESOL, nullptr, nullptr) );

   (*probdata)->degrees = nullptr;
   (*probdata)->localdegrees = nullptr;
   (*probdata)->nlocalones = 0;
   (*probdata)->nlocalzeros = 0;
   (*probdata)->nlocaledges = 0;
   (*probdata)->localdensity = 1.0;

   // include event handler
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, BACSeventExecLocalDegrees, (SCIP_EVENTHDLRDATA*) *probdata) );

   // get event handler
   (*probdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if ( (*probdata)->eventhdlr == nullptr )
   {
      SCIPerrorMessage("event handler for local degrees not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   (*probdata)->cliqueindex = nullptr;
   (*probdata)->cliquesizes = nullptr;
   (*probdata)->ncliques = 0;
   (*probdata)->maxcliqueindex = 0;
   (*probdata)->cliquenode = nullptr;

   return SCIP_OKAY;
}


//! free problem data
SCIP_RETCODE BACSfreeProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< problem data structure */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );

   if ( (*probdata)->Gorig != nullptr )
      delete (*probdata)->Gorig;
   if ( (*probdata)->changedpresol && (*probdata)->G != nullptr )
      delete (*probdata)->G;
   if ( (*probdata)->isolated != nullptr )
      delete (*probdata)->isolated;
   if ( (*probdata)->CG != nullptr )
      delete (*probdata)->CG;

   for (unsigned int i = 0; i < (*probdata)->n; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->vars, (*probdata)->n);

   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->degrees, (*probdata)->n);
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->localdegrees, (*probdata)->n);
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->origfixings, (*probdata)->norig);
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->origtopresol, (*probdata)->norig);
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->origtopresolparity, (*probdata)->norig);
   if ( (*probdata)->mergemapping != nullptr )
      SCIPfreeDisjointset(scip, &(*probdata)->mergemapping);

   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->cliqueindex, (*probdata)->n);
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->cliquesizes, (*probdata)->n);

   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}


//! compute degrees for global graph and store degrees in probdata, setup isolated
static
SCIP_RETCODE computeDegrees(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( probdata->degrees == nullptr );

   const Graph* G = probdata->G;
   size_t n = probdata->n;
   assert( n == boost::num_vertices(*G) );

   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(probdata->degrees), n) );

   EdgeIterator eit, eend;
   for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
   {
      Vertex s = boost::source(*eit, *G);
      Vertex t = boost::target(*eit, *G);

      assert( s < n );
      assert( t < n );

      ++(probdata->degrees[s]);
      ++(probdata->degrees[t]);
   }

   assert( probdata->isolated == nullptr );
   probdata->isolated = new std::vector<bool>(n, false);
   size_t nisolated = 0;
   for (size_t i = 0; i < n; ++i)
   {
      if ( probdata->degrees[i] == 0 )
      {
         (*probdata->isolated)[i] = true;
         ++nisolated;
      }
   }
   if ( nisolated == 0 )
   {
      delete probdata->isolated;
      probdata->isolated = nullptr;
   }

   return SCIP_OKAY;
}


/** check whether local degrees are correct */
SCIP_RETCODE BACScheckLocalDegrees(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( probdata->localdegrees != nullptr );

   const Graph* G = probdata->G;
   size_t n = probdata->n;
   assert( n == boost::num_vertices(*G) );

   // compute current degrees
   size_t nedges = 0;
   int* mydeg;
   SCIP_CALL( SCIPallocClearBufferArray(scip, &mydeg, n) );

   EdgeIterator eit, eend;
   for (boost::tie(eit, eend) = boost::edges(*G); eit != eend; ++eit)
   {
      Vertex s = boost::source(*eit, *G);
      Vertex t = boost::target(*eit, *G);

      assert( s < n );
      assert( t < n );

      // skip edges with one fixed variable
      if ( SCIPvarGetUbLocal(probdata->vars[s]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[s]) > 0.5 )
         continue;

      if ( SCIPvarGetUbLocal(probdata->vars[t]) < 0.5 || SCIPvarGetLbLocal(probdata->vars[t]) > 0.5 )
         continue;

      ++(mydeg[s]);
      ++(mydeg[t]);
      ++nedges;
   }
   assert( nedges == probdata->nlocaledges );

   size_t nf0 = 0;
   size_t nf1 = 0;
   for (size_t i = 0; i < n; ++i)
   {
      if ( SCIPvarGetUbLocal(probdata->vars[i]) < 0.5 )
      {
         ++nf0;
         assert( probdata->localdegrees[i] == -1 );
      }
      else if ( SCIPvarGetLbLocal(probdata->vars[i]) > 0.5 )
      {
         ++nf1;
         assert( probdata->localdegrees[i] == -1 );
      }
      else
         assert( mydeg[i] == probdata->localdegrees[i] );
   }
   assert( nf0 == probdata->nlocalzeros );
   assert( nf1 == probdata->nlocalones );
   assert( n >= nf0 + nf1 );

   // check density
   if ( n - nf0 - nf1 <= 1 )
   {
      // at most one free node
      assert( nedges == 0 );
      assert( SCIPisEQ(scip, probdata->localdensity, 1.0) );
   }
   else
      assert( SCIPisEQ(scip, probdata->localdensity, (double) nedges / ((n - nf0 - nf1) * (n - nf0 - nf1 - 1) / 2)) );

   SCIPfreeBufferArray(scip, &mydeg);

   return SCIP_OKAY;
}

//! get next number from string s
static
SCIP_Longint getNextNumber(
   char*&                s                   //!< string pointer (modified)
   )
{
   // skip whitespace
   while ( isspace((unsigned char) *s) )
      ++s;

   // get number
   SCIP_Longint tmp = std::atol(s);

   // skip number
   while ( (*s != 0) && (! isspace((unsigned char) *s)) )
      ++s;

   return tmp;
}


//! check whether the graph is triangle free
static
bool graphIsTriangleFree(
   const Graph&          G                   /**< graph */
   )
{
   // prepare clique marker
   size_t n = boost::num_vertices(G);
   std::vector<int> cand(n, -1);

   // loop through all edges
   size_t edgenum = 0;
   EdgeIterator eit, eend;
   for (boost::tie(eit, eend) = boost::edges(G); eit != eend; ++eit)
   {
      Edge e = *eit;
      Vertex source = boost::source(e, G);
      Vertex target = boost::target(e, G);

      // mark neighbors of source
      AdjacencyIterator ait, aend;
      for (boost::tie(ait, aend) = boost::adjacent_vertices(source, G); ait != aend; ++ait)
      {
         Vertex w = *ait;
         cand[w] = (int) edgenum;
      }

      // check whether we can extend the clique
      for (boost::tie(ait, aend) = boost::adjacent_vertices(target, G); ait != aend; ++ait)
      {
         Vertex w = *ait;
         assert( w != target );

         // if we found a common neighbor
         if ( cand[w] == (int) edgenum )
            return false;
      }
      ++edgenum;
   }
   return true;
}


//! read graph
SCIP_RETCODE BACSreadGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   const char*           filename            /**< name of graph file to read */
   )
{
   // should be run on fresh problem data
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( filename != nullptr );
   assert( probdata->Gorig == nullptr );
   assert( probdata->vars == nullptr );

   // open file
   SCIP_FILE* file = SCIPfopen(filename, "r");
   if ( ! file )
   {
      SCIPerrorMessage("Error: could not open %s\n", filename);
      return SCIP_NOFILE;
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Reading file %s ...\n", filename);

   // get buffer for reading input lines
   char* buffer = nullptr;
   SCIP_CALL( SCIPallocBufferArray(scip, &buffer, SCIP_MAXSTRLEN + 1) );

   // get number of nodes and edges: find line starting with 'p'
   unsigned int line = 1;
   char* str = buffer;
   do
   {
      (void) SCIPfgets(buffer, SCIP_MAXSTRLEN, file);
      ++line;

      // skip whitespace
      str = buffer;
      while ( isspace((unsigned char) *str) )
         ++str;
   }
   while ( ! SCIPfeof(file) && *str != 'p' );

   if ( SCIPfeof(file) )
   {
      SCIPerrorMessage("%s: Could not find line starting with 'p'.\n", filename);
      SCIPfreeBufferArray(scip, &buffer);
      (void) SCIPfclose(file);
      return SCIP_INVALIDDATA;
   }

   // skip 'p'
   ++str;

   // skip whitespace
   while ( isspace((unsigned char) *str) )
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
         SCIPerrorMessage("%s: line %u: Line starting with 'p' must continue with 'edge' or 'col'!\n%s\n", filename, line, buffer);
         SCIPfreeBufferArray(scip, &buffer);
         (void) SCIPfclose(file);
         return SCIP_INVALIDDATA;
      }
   }

   // read sizes
   SCIP_Longint nn = getNextNumber(str);
   SCIP_Longint mm = getNextNumber(str);

   if ( nn <= 0 )
   {
      SCIPerrorMessage("%s: line %u: Number of nodes must be positive!\n%s\n", filename, line, buffer);
      SCIPfreeBufferArray(scip, &buffer);
      (void) SCIPfclose(file);
      return SCIP_INVALIDDATA;
   }

   if ( mm < 0 )
   {
      SCIPerrorMessage("%s: line %u: Number of edges must be nonnegative!\n%s\n", filename, line, buffer);
      SCIPfreeBufferArray(scip, &buffer);
      (void) SCIPfclose(file);
      return SCIP_INVALIDDATA;
   }
   size_t n = (size_t) nn;
   size_t m = (size_t) mm;

   // set to store edges to avoid duplicates
   std::set<std::pair<size_t, size_t> > EdgeSet;

   // start to read nodes and edges into set
   (void) SCIPfgets(buffer, SCIP_MAXSTRLEN, file);
   ++line;
   assert( ! SCIPfeof(file) );

   // detect loops and duplicate edges
   std::vector<SCIP_Real> weights(n, 0.0); // default value for the weighted version is 0
   bool foundweights = false;
   bool foundnegativeweights = false;
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
         // check node
         if ( node < 0 || node >= nn )
         {
            SCIPerrorMessage("%s: line %u: Node number invalid!\n%s\n", filename, line, buffer);
            SCIPfreeBufferArray(scip, &buffer);
            return SCIP_INVALIDDATA;
         }
         assert( node >= 0 );

         // skip spaces
         while ( isspace((unsigned char) *str) )
            ++str;

         // check wether weight is present
         if ( *str == '\0' )
         {
            SCIPerrorMessage("%s: line %u: weight not specified!\n%s\n", filename, line, buffer);
            SCIPfreeBufferArray(scip, &buffer);
            return SCIP_INVALIDDATA;
         }

         SCIP_Real weight = atof(str);
         weights[node] = weight;  /*lint !e732*/
         foundweights = true;
         if ( SCIPisLT(scip, weight, 0.0) )
            foundnegativeweights = true;

         // check for infinite weights
         if ( SCIPisInfinity(scip, weight) )
         {
            SCIPerrorMessage("%s: line %u: weight for node is considered infinite; cannot currently deal with this.\n%s\n", filename, line, buffer);
            SCIPfreeBufferArray(scip, &buffer);
            return SCIP_INVALIDDATA;
         }
      }
      else if ( str[0] == 'e' )
      {
         ++str;
         SCIP_Longint source = getNextNumber(str) - 1;
         SCIP_Longint target = getNextNumber(str) - 1;

         // check source
         if ( source < 0 || source >= nn || target < 0 || target >= nn )
         {
            SCIPerrorMessage("%s: line %u: Node number invalid!\n%s\n", filename, line, buffer);
            SCIPfreeBufferArray(scip, &buffer);
            return SCIP_INVALIDDATA;
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

   // create nodes
   probdata->Gorig = new Graph();
   Graph& G = *probdata->Gorig;
   probdata->unweighted = TRUE;
   for (size_t i = 0; i < n; ++i)
   {
      Vertex v = boost::add_vertex(G);

      if ( foundweights )
      {
         if ( ! SCIPisEQ(scip, weights[v], 1.0) )
            probdata->unweighted = FALSE;
         boost::put(vertex_weight_t(), G, v, weights[v]);
      }
      else
         boost::put(vertex_weight_t(), G, v, 1);  // default weight is 1 for the unweighted version
   }
   assert( boost::num_vertices(G) == n );

   if ( nloops > 0 )
   {
      SCIPwarningMessage(scip, "Graph contains %zu loops (ignored)!\n", nloops);
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
         SCIPerrorMessage("Could not create edge {%lu,%lu}.\n", s, t);
         SCIPfreeBufferArray(scip, &buffer);
         return SCIP_INVALIDDATA;
      }
   }
   assert( boost::num_edges(G) == nedges );

   if ( nduplicateedges > 0 )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, nullptr, "Duplicate edges: %zu\n", nduplicateedges);

   if ( nedges + nduplicateedges + nloops != m )
   {
      SCIPerrorMessage("%s: line %u: Found %zu edges. There should be %zu.\n", filename, line, nedges + nduplicateedges + nloops, m);
      SCIPfreeBufferArray(scip, &buffer);
      return SCIP_INVALIDDATA;
   }

   SCIPfreeBufferArray(scip, &buffer);

   probdata->norig = n;
   probdata->morig = nedges;

   SCIPinfoMessage(scip, nullptr, "Number of nodes:\t\t %zu\n", probdata->norig);
   SCIPinfoMessage(scip, nullptr, "Number of edges:\t\t %zu\n", probdata->morig);

   if ( probdata->unweighted )
      SCIPinfoMessage(scip, nullptr, "Is weighted:\t\t\t False\n");
   else
      SCIPinfoMessage(scip, nullptr, "Is weighted:\t\t\t True\n");

   // print statistics
   if ( probdata->printprobstats )
   {
      // components stores the index of the component for each vertex
      std::vector<int> components(boost::num_vertices(G));
      // compute connected components in graph
      int ncomp = connected_components(G, &components[0]);

      SCIP_Real avdeg = (SCIP_Real) nedges / n * 2;
      size_t maxdeg = 0;
      size_t deg = 0;
      size_t nisolatednodes = 0;
      for (size_t i = 0; i < n; ++i)
      {
         deg = boost::degree(i, G);
         if ( deg > maxdeg )
            maxdeg = deg;
         if ( deg == 0 )
            ++nisolatednodes;
      }
      assert( SCIPisLE(scip, avdeg, maxdeg) );

      SCIP_Real nsquared = (SCIP_Real) n * (n - 1);
      SCIP_Real density = nedges / (nsquared / 2);

      assert( SCIPisLE(scip, density, 1.0) );
      assert( SCIPisGE(scip, density, 0.0) );

      if ( foundnegativeweights )
         SCIPinfoMessage(scip, nullptr, "Has negative weights:\t\t True\n");
      else
         SCIPinfoMessage(scip, nullptr, "Has negative weights:\t\t False\n");
      SCIPinfoMessage(scip, nullptr, "Number of loops:\t\t %zu\n", nloops);
      SCIPinfoMessage(scip, nullptr, "Number of isolated nodes:\t %zu\n", nisolatednodes);
      SCIPinfoMessage(scip, nullptr, "Density of graph:\t\t %.4f\n", density);
      SCIPinfoMessage(scip, nullptr, "Number of connected components:\t %d\n", ncomp);
      SCIPinfoMessage(scip, nullptr, "Maximum degree:\t\t\t %zu\n", maxdeg);
      SCIPinfoMessage(scip, nullptr, "Average degree:\t\t\t %.1f\n", avdeg);
      if ( graphIsTriangleFree(G) )
         SCIPinfoMessage(scip, nullptr, "Graph is triangle free.\n");
      else
         SCIPinfoMessage(scip, nullptr, "Graph is not triangle free.\n");
   }

   return SCIP_OKAY;
}



//! find neighbors of degree 2 node
static
void BACSfindNeighborsDegreeTwo(
   const Graph*          G,                  /**< graph */
   const Graph*          E,                  /**< extended graph to be processed */
   Vertex                v,                  /**< degree 2 node */
   BACS_NODEFIXING*      fixed,              /**< array to mark fixed nodes */
   size_t*               degrees,            /**< degrees */
   Vertex&               s,                  /**< returns one neighbor */
   Vertex&               t                   /**< returns other neighbor */
   )
{
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );
   assert( degrees[v] == 2 );

   // init
   s = Graph::null_vertex();
   t = Graph::null_vertex();

   // loop through neighbors of given node
   bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(v, G, E);
   for (; ! ait.at_end(); ++ait)
   {
      Vertex w = *ait;
      if ( fixed[w] != BACSunfixed )
         continue;

      if ( s == Graph::null_vertex() )
         s = w;
      else
      {
         assert( t == Graph::null_vertex() );
         t = w;
      }
   }
   assert( s != Graph::null_vertex() );
   assert( t != Graph::null_vertex() );
}


/** collect nodes on path of degree 2 nodes, starting at @a v and going along @a w
 *
 *  The path can close itself, in which case we have found a cycle. The node @a v is part of the path if and only if we
 *  have found a cycle.
 */
static
void BACScollectPathDegreeTwo(
   const Graph*          G,                  /**< graph */
   const Graph*          E,                  /**< extended graph to be processed */
   Vertex                v,                  /**< start degree 2 node */
   Vertex                w,                  /**< next node */
   BACS_NODEFIXING*      fixed,              /**< array to mark fixed nodes */
   size_t*               degrees,            /**< degrees */
   std::vector<Vertex>&  path,               /**< returns the path */
   SCIP_Bool&            foundcycle          /**< returns wether we have a cycle */
   )
{
   assert( G != nullptr );
   assert( fixed != nullptr );
   assert( degrees != nullptr );
   assert( degrees[v] == 2 );

   foundcycle = FALSE;
   path.clear();

   // initialize path, start with w
   path.push_back(w);

   // follow path by looping as long as we have degree 2 nodes
   Vertex neighbor1 = v;
   Vertex curnode = w;
   while ( degrees[curnode] == 2 )
   {
      assert( boost::get(vertex_weight_t(), *G, curnode) == boost::get(vertex_weight_t(), *G, v) );

      // check whether we found a cycle
      if ( curnode == v )
      {
         // found cycle
         foundcycle = TRUE;
         break;
      }

      // find other neighbor of s
      Vertex neighbor2 = Graph::null_vertex();
      bacs::ExtAdjacencyIterator ait = bacs::adjacent_vertices(curnode, G, E);
      for (; ! ait.at_end(); ++ait)
      {
         Vertex u = *ait;
         if ( fixed[u] != BACSunfixed )
            continue;

         if ( u == neighbor1 )
            continue;

         assert( neighbor2 == Graph::null_vertex() );
         neighbor2 = u;
      }
      assert( neighbor2 != Graph::null_vertex() );
      assert( neighbor2 != curnode );

      neighbor1 = curnode;
      curnode = neighbor2;

      path.push_back(curnode);
   }
   assert( curnode != Graph::null_vertex() );
   assert( degrees[curnode] != 2 || foundcycle );
}


//! perform graph presolving
SCIP_RETCODE BACSgraphPresolving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   SCIP_Real             timelimit           /**< time limit for presolving */
   )
{
   // possibly skip graph presolving
   if ( ! probdata->graphpresolving || probdata->norig == 0 )
   {
      probdata->G = probdata->Gorig;
      probdata->changedpresol = FALSE;
      probdata->n = probdata->norig;
      probdata->m = probdata->morig;
      return SCIP_OKAY;
   }

   // start time measurement
   struct tms timer_beg;
   (void) times(&timer_beg);

   const Graph* G = probdata->Gorig;
   size_t n = probdata->norig;

#ifndef NDEBUG
   {
      // assert AdjacencyIterator gives nodes in ascending order as this is crucial for following neighborhood search
      VertexIterator vit, vend;
      for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
      {
         AdjacencyIterator nit, nend;
         boost::tie(nit, nend) = boost::adjacent_vertices(*vit, *G);
         Vertex nprev;
         while ( nit != nend )
         {
            nprev = *nit;
            ++nit;
            assert( nit == nend || nprev < *nit );
         }
      }
   }
#endif

   SCIP_CALL( SCIPcreateDisjointset(scip, &probdata->mergemapping, (int) n) );
   SCIP_DISJOINTSET* mergemapping = probdata->mergemapping;

   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->origfixings, n) );
   BACS_NODEFIXING* fixed = probdata->origfixings;

   SCIP_Real* curweights;
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &curweights, n) );

   for (size_t i = 0; i < n; ++i)
      curweights[i] = boost::get(vertex_weight_t(), *G, (Vertex) i);

   SCIPinfoMessage(scip, nullptr, "\n");
   SCIPinfoMessage(scip, nullptr, "------------------------------------------------------------------------------------\n");
   SCIPinfoMessage(scip, nullptr, "Running graph presolving ...\n");

   Graph E;   // graph of new nodes
   size_t nfixedzero = 0;
   size_t nfixedone = 0;
   size_t naddededges = 0;
   size_t nfixedzeroneigh = 0;
   size_t nfixedzerosimpl = 0;
   size_t nfixedzerocycle = 0;
   size_t nfixedonesimpl = 0;
   size_t nfixedonecycle = 0;
   size_t nfixedzerosst = 0;
   size_t naddededgessst = 0;
   size_t naddededgesinprob = 0;

   SCIP_Bool isweighted;
   if ( probdata->unweighted )
      isweighted = FALSE;
   else
      isweighted = TRUE;

   // first make sure that all unfixed nodes have positive weight
   if ( isweighted )
   {
      VertexIterator vit, vend;
      for (boost::tie(vit, vend) = boost::vertices(*G); vit != vend; ++vit)
      {
         Vertex v = *vit;
         // nodes with nonpositive weights can be fixed to be 0
         if ( SCIPisLE(scip, curweights[v], 0.0) )
         {
            fixed[v] = BACSfixedzero;
            ++nfixedzero;
         }
      }
      SCIPinfoMessage(scip, nullptr, "Fixed %lu nodes to zero with nonpositive weight.\n", nfixedzero);
   }

   // determine which presolving methods to run
   SCIP_Bool neighborhoodpresol = probdata->neighborhoodpresol;
   SCIP_Bool simplicialpresol = probdata->simplicialpresol;
   SCIP_Bool cyclepresol = probdata->cyclepresol;
   SCIP_Bool mergepresol = probdata->mergepresol;
   if ( isweighted )
      cyclepresol = FALSE;

   size_t* degrees = nullptr;
   size_t ndegreezero;
   size_t ndegreeone;
   size_t ndegreetwo;
   size_t nrounds = 0;
   size_t nmerged = 0;

   // check time
   if ( ! SCIPisInfinity(scip, timelimit) )
   {
      struct tms timer_end;
      (void) times(&timer_end);
      double t = (timer_end.tms_utime - timer_beg.tms_utime) / (double)sysconf(_SC_CLK_TCK);
      if ( t > timelimit )
         goto TERMINATE;
   }

   // get degrees
   SCIP_CALL( SCIPallocBufferArray(scip, &degrees, n) );
   BACScomputeDegrees(G, fixed, degrees, ndegreezero, ndegreeone, ndegreetwo);

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

      SCIP_CALL( BACSpresolvingRound(scip, G, isweighted, neighborhoodpresol, simplicialpresol, mergepresol, cyclepresol,
            fixed, degrees, curweights, mergemapping, ndegreezero, ndegreeone, ndegreetwo, nroundfixedzeroneigh, nroundfixedzerosimpl, nroundfixedzerocycle,
            nroundfixedonesimpl, nroundfixedonecycle, nroundmerged) );

      SCIPinfoMessage(scip, nullptr, "Round %lu: neighborhood (z: %lu), merge (z: %lu), simplicial (z: %lu, o: %lu), cycle (z: %lu, o: %lu).\n",
         nrounds, nroundfixedzeroneigh, nroundmerged, nroundfixedzerosimpl, nroundfixedonesimpl, nroundfixedzerocycle, nroundfixedonecycle);

      if ( nroundmerged > 0 )
         cyclepresol = FALSE;

      nfixedzero += nroundfixedzeroneigh + nroundfixedzerosimpl + nroundfixedzerocycle + nroundmerged;
      nfixedone += nroundfixedonesimpl + nroundfixedonecycle;
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

      assert( nfixedzero + nfixedone <= n );
      if ( nfixedzero + nfixedone == n )
         break;
   }
   while ( changed );

   // if there are unfixed nodes and we want to perform inprobing or SST presolving
   if ( (nfixedzero + nfixedone < n) && (probdata->graphinprobing || probdata->sstpresolving) )
   {
      for (size_t i = 0; i < n; ++i)
      {
         (void) boost::add_vertex(E);
         boost::put(vertex_weight_t(), E, (Vertex) i, curweights[i]);
      }

      // run in-probing
      if ( probdata->graphinprobing )
      {
         size_t nlocaladdededges;
         SCIP_CALL( BACSgraphInProbing(scip, timer_beg, timelimit, G, isweighted, fixed, degrees, curweights, mergemapping,
               ndegreezero, ndegreeone, ndegreetwo, &E, nlocaladdededges) );
         SCIPinfoMessage(scip, nullptr, "In-probing found %lu edges.\n", nlocaladdededges);
         assert( nlocaladdededges == boost::num_edges(E) );
         naddededges += nlocaladdededges;
         naddededgesinprob += nlocaladdededges;

         // run extended presolving if we found edges
         if ( nlocaladdededges > 0 )
         {
            size_t nlocalfixedzeroneigh;
            size_t nlocalfixedzerosimpl;
            size_t nlocalfixedzerocycle;
            size_t nlocalfixedonesimpl;
            size_t nlocalfixedonecycle;
            size_t nlocalmerged;

            SCIP_CALL( BACSpresolvingRoundsExtend(scip, G, &E, isweighted, neighborhoodpresol, simplicialpresol, mergepresol, cyclepresol, fixed, degrees, curweights,
                  mergemapping, ndegreezero, ndegreeone, ndegreetwo,
                  nrounds, nlocalfixedzeroneigh, nlocalfixedzerosimpl, nlocalfixedzerocycle, nlocalfixedonesimpl, nlocalfixedonecycle, nlocalmerged) );

            if ( nlocalmerged > 0 )
               cyclepresol = FALSE;

            nfixedzero += nlocalfixedzeroneigh + nlocalfixedzerosimpl + nlocalfixedzerocycle + nlocalmerged;
            nfixedone += nlocalfixedonesimpl + nlocalfixedonecycle;
            nfixedzeroneigh += nlocalfixedzeroneigh;
            nfixedzerosimpl += nlocalfixedzerosimpl;
            nfixedzerocycle += nlocalfixedzerocycle;
            nfixedonesimpl += nlocalfixedonesimpl;
            nfixedonecycle += nlocalfixedonecycle;
            nmerged += nlocalmerged;
         }
      }
      assert( nfixedzero + nfixedone <= n );

      // perform SST presolving
      if ( probdata->sstpresolving && (nfixedzero + nfixedone < n) )
      {
         // allocate memory for sorted components
         int* components;
         SCIP_CALL( SCIPallocBufferArray(scip, &components, n) );
         size_t ncomponents;

         do
         {
            // determine connected components by BFS, ignoring fixed nodes
            SCIP_CALL( BFSExtended(scip, G, &E, fixed, components, ncomponents) );

            size_t nlocalfixedzero;
            size_t nlocaladdededges;
            size_t nlocalmerged;
            SCIPinfoMessage(scip, nullptr, "\nSST presolving ...\n");
            SCIP_CALL( BACSpresolvingSST(scip, probdata->sstaddedges, G, isweighted, (ORBIT_RULE) probdata->sstorbitrule, fixed, degrees, curweights,
                  mergemapping, components, ncomponents, ndegreezero, ndegreeone, ndegreetwo,
                  &E, nlocalfixedzero, nlocaladdededges, nlocalmerged) );
            SCIPinfoMessage(scip, nullptr, "SST presolving removed %lu nodes and added %lu edges.\n\n",
               nlocalfixedzero, nlocaladdededges);
            nfixedzero += nlocalfixedzero;
            nfixedzerosst += nlocalfixedzero;
            naddededges += nlocaladdededges;
            naddededgessst += nlocaladdededges;
            nmerged += nlocalmerged;

            // run extended presolving if SST presolving was successful
            changed = FALSE;
            if ( nlocalfixedzero + nlocaladdededges > 0 )
            {
               size_t nlocalfixedzeroneigh;
               size_t nlocalfixedzerosimpl;
               size_t nlocalfixedzerocycle;
               size_t nlocalfixedonesimpl;
               size_t nlocalfixedonecycle;

               SCIP_CALL( BACSpresolvingRoundsExtend(scip, G, &E, isweighted, neighborhoodpresol, simplicialpresol, mergepresol, cyclepresol, fixed, degrees, curweights,
                     mergemapping, ndegreezero, ndegreeone, ndegreetwo, nrounds,
                     nlocalfixedzeroneigh, nlocalfixedzerosimpl, nlocalfixedzerocycle, nlocalfixedonesimpl, nlocalfixedonecycle, nlocalmerged) );

               if ( nlocalmerged > 0 )
                  cyclepresol = FALSE;

               nfixedzero += nlocalfixedzeroneigh + nlocalfixedzerosimpl + nlocalfixedzerocycle + nlocalmerged;
               nfixedone += nlocalfixedonesimpl + nlocalfixedonecycle;
               nfixedzeroneigh += nlocalfixedzeroneigh;
               nfixedzerosimpl += nlocalfixedzerosimpl;
               nfixedzerocycle += nlocalfixedzerocycle;
               nfixedonesimpl += nlocalfixedonesimpl;
               nfixedonecycle += nlocalfixedonecycle;
               nmerged += nlocalmerged;

               if ( probdata->sstrepeat && nlocalfixedzero + nlocaladdededges + nlocalfixedzeroneigh + nlocalfixedzerosimpl + nlocalfixedzerocycle + nlocalfixedonesimpl + nlocalfixedonecycle + nlocalmerged > 0 )
                  changed = TRUE;
            }
            assert( nfixedzero + nfixedone <= n );
            if ( nfixedzero + nfixedone == n )
               break;
         }
         while ( changed );

         SCIPfreeBufferArray(scip, &components);

      }
      assert( nfixedzero + nfixedone <= n );
   }

 TERMINATE:
   // if reductions have been found, create presolved graph
   size_t ncontracted = 0;
   if ( nfixedzero + nfixedone + ndegreetwo + boost::num_edges(E) + nmerged > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->origtopresol, n) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->origtopresolparity, n) );
      int* origtopresol = probdata->origtopresol;
      int* origtopresolparity = probdata->origtopresolparity;
      for (size_t i = 0; i < n; ++i)
      {
         origtopresol[i] = -1;
         origtopresolparity[i] = 0;
      }
      probdata->G = new Graph();
      Graph* H = probdata->G;

      // first treat degree 2 nodes along paths (cycles have been treated above)
      if ( ndegreetwo > 0 && ! isweighted && probdata->degtwocontract )
      {
         // make sure there are no degree 0 or 1 nodes (might happen if presolving above is incomplete)
         while ( ndegreezero + ndegreeone > 0 )
         {
            size_t nlocalfixedzero;
            size_t nlocalfixedone;

            if ( boost::num_vertices(E) > 0 )
               BACSsimplicialPresolving(scip, G, &E, fixed, degrees, curweights, ndegreezero, ndegreeone, ndegreetwo, nlocalfixedzero, nlocalfixedone);
            else
               BACSsimplicialPresolving(scip, G, nullptr, fixed, degrees, curweights, ndegreezero, ndegreeone, ndegreetwo, nlocalfixedzero, nlocalfixedone);
            nfixedzero += nlocalfixedzero;
            nfixedone += nlocalfixedone;
         }
         assert( ndegreezero + ndegreeone == 0 );

         for (size_t i = 0; i < n; ++i)
         {
            Vertex v = (Vertex) i;
            if ( degrees[v] != 2 || origtopresol[v] >= 0 || fixed[v] != BACSunfixed )
               continue;
            assert( origtopresol[v] == -1 );
            assert( curweights[v] == 1.0 );

            // find the two neighbors of current node (only w.r.t. G because E was computed later)
            Vertex s;
            Vertex t;
            BACSfindNeighborsDegreeTwo(G, nullptr, v, fixed, degrees, s, t);

            // collect nodes on the path to the beginning (only w.r.t. G because E was computed later)
            std::vector<Vertex> path1;
            SCIP_Bool foundcycle;
            BACScollectPathDegreeTwo(G, nullptr, v, s, fixed, degrees, path1, foundcycle);
            if ( foundcycle )  // cycles could have been treated above
               continue;

            // find path to the end (only w.r.t. G because E was computed later)
            std::vector<Vertex> path2;
            BACScollectPathDegreeTwo(G, nullptr, v, t, fixed, degrees, path2, foundcycle);
            assert( ! foundcycle );

            // ignore cycles with one non-degree-two node (this should be solved in the cycle presolving)
            if ( path1.back() == path2.back() )
               continue; // to improve performance one could skip all degree two nodes in path1 and path2 now

            assert( path1.size() > 0 );
            assert( path2.size() > 0 );
            assert( degrees[path1.back()] > 2 );
            assert( degrees[path2.back()] > 2 );
            size_t pathsize = path1.size() + path2.size() + 1;
            assert( pathsize >= 3 );

            // ignore odd cycles with adjacent end nodes (this should be solved in the cycle presolving)
            if ( pathsize % 2 == 1 && boost::edge(path1.back(), path2.back(), *G).second )
               continue; // to improve performance one could skip all degree two nodes in path1 and path2 now

            // treat one end node
            Vertex d;
            if ( origtopresol[path1.back()] < 0 )
            {
               // create representative node in new graph
               d = boost::add_vertex(*H);
               boost::put(vertex_weight_t(), *H, d, 1.0);
               origtopresol[path1.back()] = (int) d;
               ++ncontracted;
            }
            else
            {
               assert( origtopresol[path1.back()] >= 0 );
               d = (Vertex) origtopresol[path1.back()]; /*lint !e571*/
            }
            assert( origtopresol[path1.back()] >= 0 );

            // treat other end node
            if ( origtopresol[path2.back()] < 0 )
            {
               Vertex e;
               if ( pathsize % 2 == 0 )
               {
                  // create representative node in new graph
                  e = boost::add_vertex(*H);
                  boost::put(vertex_weight_t(), *H, e, 1.0);
               }
               else
                  e = d;
               origtopresol[path2.back()] = (int) e;
               ++ncontracted;
            }
            assert( origtopresol[path2.back()] >= 0 );

            // mark current node
            assert( degrees[v] == 2 );
            origtopresol[v] = (int) d;
            ++ncontracted;
            int startparity = 1;
            if ( path1.size() % 2 == 1 )
               startparity = -1;
            origtopresolparity[v] = startparity;

            // mark path1 w.r.t. d
            std::vector<Vertex>::const_iterator pit;
            std::vector<Vertex>::const_iterator pend = path1.end();
            int parity = -startparity;
            for (pit = path1.begin(); pit != pend; ++pit)
            {
               if ( degrees[*pit] > 2 )
                  break;
               assert( degrees[*pit] == 2 );
               assert( *pit < n );
               origtopresol[*pit] = (int) d;
               ++ncontracted;
               origtopresolparity[*pit] = parity;
               parity = -parity;
            }

            // mark path2 w.r.t. d
            parity = -startparity;
            pend = path2.end();
            for (pit = path2.begin(); pit != pend; ++pit)
            {
               if ( degrees[*pit] > 2 )
                  break;
               assert( degrees[*pit] == 2 );
               assert( *pit < n );
               origtopresol[*pit] = (int) d;
               ++ncontracted;
               origtopresolparity[*pit] = parity;
               parity = -parity;
            }

            // update objective
            if ( pathsize % 2 == 0 )
               probdata->objoffset += (int) pathsize / 2 - 1;
            else
               probdata->objoffset += (int) ((pathsize - 1) / 2);

         }
      }

      // transfer nodes to new graph
      // check whether presolved graph is unweighted
      for (size_t i = 0; i < n; ++i)
      {
         Vertex v = (Vertex) i;

         // skip already treated nodes
         if ( origtopresol[v] >= 0 )
            continue;

         if ( fixed[v] == BACSunfixed )
         {
            origtopresol[v] = (int) boost::add_vertex(*H);
            boost::put(vertex_weight_t(), *H, origtopresol[v], curweights[v]);

            if ( ! SCIPisEQ(scip, curweights[v], 1.0) )
               isweighted = TRUE;
         }
         else if ( fixed[v] == BACSfixedone )
            probdata->objoffset += curweights[v]; // update objective
         // do nothing for nodes fixe to zero
      }

      // transfer edges - to ensure sorting first collect edges in a set
      std::set<std::pair<size_t, size_t> > EdgeSet;
      bacs::ExtEdgeIterator eit = bacs::ExtEdgeIterator(G, &E);
      for (; ! eit.at_end(); ++eit)
      {
         Vertex s = eit.source();
         Vertex t = eit.target();

         if ( origtopresol[s] >= 0 && origtopresol[t] >= 0 )
         {
            assert( fixed[s] == BACSunfixed );
            assert( fixed[t] == BACSunfixed );

            size_t sh = (size_t) origtopresol[s];
            size_t th = (size_t) origtopresol[t];

            // skip loops (can appear due to contraction)
            if ( sh != th )
            {
               if ( sh < th )
                  (void) EdgeSet.insert(std::make_pair(sh, th));  // we ignore duplicate edges here
               else
                  (void) EdgeSet.insert(std::make_pair(th, sh));  // we ignore duplicate edges here
            }
         }
      }

      // now create edges
      size_t nedges = 0;
      std::set<std::pair<size_t, size_t> >::const_iterator sit, send = EdgeSet.end();
      for (sit = EdgeSet.begin(); sit != send; ++sit)
      {
         size_t s = sit->first;
         size_t t = sit->second;

         assert( ! boost::edge(s, t, *H).second );

#ifndef NDEBUG
         std::pair<Edge, bool> p = boost::add_edge(s, t, *H);
         assert( p.second );
#else
         (void) boost::add_edge(s, t, *H);
#endif
         ++nedges;
      }
      assert( boost::num_edges(*H) == nedges );

      // update number of nodes and edges
      probdata->n = boost::num_vertices(*H);
      probdata->m = boost::num_edges(*H);
      probdata->changedpresol = TRUE;
      if ( isweighted )
         probdata->unweighted = FALSE;
      else
         probdata->unweighted = TRUE;
   }
   else
   {
      // otherwise just copy graph
      probdata->G = probdata->Gorig;
      probdata->changedpresol = FALSE;
      probdata->n = probdata->norig;
      probdata->m = probdata->morig;
   }

   SCIPfreeBufferArrayNull(scip, &degrees);
   SCIPfreeBlockMemoryArray(scip, &curweights, n);

   SCIPinfoMessage(scip, nullptr, "\nTotal fixings: neighborhood (z: %lu), merge (z: %lu), simplicial (z: %lu, o: %lu), cycle (z: %lu, o: %lu), SST (z: %lu).\n",
      nfixedzeroneigh, nmerged, nfixedzerosimpl, nfixedonesimpl, nfixedzerocycle, nfixedonecycle, nfixedzerosst);
   SCIPinfoMessage(scip, nullptr, "Total added edges: inprobing: %lu, SST: %lu.\n", naddededgesinprob, naddededgessst);
   SCIPinfoMessage(scip, nullptr, "Graph presolving fixed %lu nodes to 0 and %lu nodes to 1, contracted %lu nodes, merged %lu nodes, and added %lu edges.\n\n",
      nfixedzero, nfixedone, ncontracted, nmerged, naddededges);
   assert( probdata->n == boost::num_vertices(*probdata->G) );
   assert( probdata->m == boost::num_edges(*probdata->G) );
   SCIPinfoMessage(scip, nullptr, "Number of nodes of presolved graph:\t %zu\n", probdata->n);
   SCIPinfoMessage(scip, nullptr, "Number of edges of presolved graph:\t %zu\n", probdata->m);
   if ( probdata->n == 0 && probdata->m == 0 )
      SCIPinfoMessage(scip, nullptr, "Done after presolving.\n");
   if ( probdata->changedpresol )
   {
      if ( ! isweighted )
         SCIPinfoMessage(scip, nullptr, "Is weighted:\t\t\t\t False\n");
      else
         SCIPinfoMessage(scip, nullptr, "Is weighted:\t\t\t\t True\n");
   }

   if ( probdata->printprobstats && probdata->n > 0 )
   {
      // components stores the index of the component for each vertex
      std::vector<int> components(probdata->n);
      // compute connected components in graph
      int ncomp = connected_components(*probdata->G, &components[0]);

      SCIP_Real avdeg = (SCIP_Real) probdata->m / probdata->n * 2;
      size_t maxdeg = 0;
      size_t deg = 0;
      size_t nisolatednodes = 0;
      for (size_t i = 0; i < probdata->n; ++i)
      {
         deg = boost::degree(i, *probdata->G);
         if ( deg > maxdeg )
            maxdeg = deg;
         if ( deg == 0 )
            ++nisolatednodes;
      }
      assert( SCIPisLE(scip, avdeg, maxdeg) );

      SCIP_Real binom = (SCIP_Real) probdata->n * (probdata->n - 1) / 2;
      SCIP_Real density = (binom > 0) ? probdata->m / binom : 0.0;

      assert( SCIPisLE(scip, density, 1.0) );
      assert( SCIPisGE(scip, density, 0.0) );

      SCIPinfoMessage(scip, nullptr, "\n");
      SCIPinfoMessage(scip, nullptr, "Number of isolated nodes:\t %zu\n", nisolatednodes);
      SCIPinfoMessage(scip, nullptr, "Density of graph:\t\t %.4f\n", density);
      SCIPinfoMessage(scip, nullptr, "Number of connected components:\t %d\n", ncomp);
      SCIPinfoMessage(scip, nullptr, "Maximum degree:\t\t\t %zu\n", maxdeg);
      SCIPinfoMessage(scip, nullptr, "Average degree:\t\t\t %.1f\n", avdeg);
      if ( graphIsTriangleFree(*probdata->G) )
         SCIPinfoMessage(scip, nullptr, "Graph is triangle free.\n");
      else
         SCIPinfoMessage(scip, nullptr, "Graph is not triangle free.\n");
   }

   SCIPinfoMessage(scip, nullptr, "------------------------------------------------------------------------------------\n");
   SCIPinfoMessage(scip, nullptr, "\n");

#ifdef BACS_WRITE_GRAPH
   // open file
   SCIP_FILE* file = SCIPfopen("presolved.dimacs.gz", "w");
   if ( ! file )
   {
      SCIPerrorMessage("Error: could not open %s\n", "presolved.dimacs");
      return SCIP_NOFILE;
   }

   SCIPfprintf(file, "p edge %zu %zu\n", probdata->n, probdata->m);
   EdgeIterator myeit, myend;
   for (boost::tie(myeit,myend) = boost::edges(*probdata->G); myeit != myend; ++myeit)
   {
      Vertex s = boost::source(*myeit, *probdata->G);
      Vertex t = boost::target(*myeit, *probdata->G);
      SCIPfprintf(file, "e %zu %zu\n", s + 1, t + 1);
   }
   SCIPfclose(file);
#endif

   return SCIP_OKAY;
}


//! init data structures
SCIP_RETCODE BACSinitDatastructures(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( probdata->G != nullptr );
   assert( probdata->vars == nullptr );
   assert( probdata->cliqueindex == nullptr );
   assert( probdata->cliquesizes == nullptr );

   // prepare local degrees
   SCIP_CALL( computeDegrees(scip, probdata) );
   assert( probdata->localdegrees == nullptr );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &probdata->localdegrees, probdata->degrees, probdata->n) );

   probdata->nlocaledges = probdata->m;
   if ( probdata->n > 1 )
      probdata->localdensity = (double) probdata->m / ( probdata->n * (probdata->n - 1) / 2 );
   else
      probdata->localdensity = 1.0;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->cliqueindex, probdata->n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->cliquesizes, probdata->n) );
   assert( probdata->cliqueindex != nullptr );
   assert( probdata->cliquesizes != nullptr );

   return SCIP_OKAY;
}


//! create complemented graph
SCIP_RETCODE BACScreateComplementedGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   )
{
   assert( scip != nullptr );
   assert( probdata != nullptr );
   assert( probdata->G != nullptr );

   const Graph& G = *probdata->G;
   assert( probdata->n == boost::num_vertices(G) );

   // create complemented graph with probdata->n nodes
   probdata->CG = new Graph(probdata->n);
   Graph& CG = *probdata->CG;
   assert( boost::num_vertices(CG) == boost::num_vertices(G) );

   std::vector<int> A(boost::num_vertices(G), -1);

   // add edges by brute force
   VertexIterator vit, vend;
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

   return SCIP_OKAY;
}


//! delete data (will automatically be called when SCIP problem is destructed)
static
SCIP_DECL_PROBDELORIG(BACSdelorig)
{
   assert( scip != nullptr );
   assert( probdata != nullptr );

   SCIP_CALL( BACSfreeProblem(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed)
 */
SCIP_DECL_PROBTRANS(BACSprobTrans)
{
   assert( sourcedata != nullptr );
   assert( sourcedata->eventhdlr != nullptr );

   // catch events
   size_t n = sourcedata->n;
   for (size_t i = 0; i < n; ++i)
   {
      SCIP_VAR* var = sourcedata->vars[i];
      SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->vars[i], &var) );
      SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBCHANGED,
            sourcedata->eventhdlr, (SCIP_EVENTDATA*) sourcedata, nullptr) );
   }
   *targetdata = sourcedata;

   return SCIP_OKAY;
}

//! setup problem (create variables and constraints)
SCIP_RETCODE BACSsetupProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   )
{
   assert( scip != nullptr );
   assert( probname != nullptr );
   assert( probdata != nullptr );
   assert( probdata->G != nullptr );
   assert( probdata->vars == nullptr );

   // create problem
   SCIP_CALL( SCIPcreateProb(scip, probname, BACSdelorig, BACSprobTrans, nullptr, nullptr, nullptr, nullptr, probdata) );
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
   SCIP_CALL( SCIPaddOrigObjoffset(scip, probdata->objoffset) );

   const Graph& G = *probdata->G;
   size_t n = probdata->n;
   assert( n == boost::num_vertices(G) );

   // create variables
   char name[SCIP_MAXSTRLEN];
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->vars, (int) n) );
   for (size_t i = 0; i < n; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x%u", i);

      SCIP_Real weight = boost::get(vertex_weight_t(), G, (Vertex) i);

      // create variables - the variable data is set to the index (= node) of the variable
      SCIP_CALL( SCIPcreateVar(scip, &probdata->vars[i], name, 0.0, 1.0, weight, SCIP_VARTYPE_BINARY, TRUE, FALSE, nullptr, nullptr, nullptr, nullptr, (SCIP_VARDATA*) i) );
      SCIP_CALL( SCIPaddVar(scip, probdata->vars[i]) );
   }

   // create clique constraint
   assert( probdata->isolated == nullptr || probdata->isolated->size() == n );
   SCIP_CONS* cons;
   SCIP_CALL( BACScreateConsClique(scip, &cons, "clique", probdata->G, probdata->isolated, (unsigned int) n, probdata->vars,
      TRUE, TRUE, TRUE, TRUE, TRUE) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

//! check whether best solution of SCIP is a stable set
SCIP_RETCODE BACScheckBestSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   )
{
   SCIP_SOL* bestsol = SCIPgetBestSol(scip);
   if ( bestsol != nullptr )
   {
      const Graph& G = *probdata->G;
      assert( probdata->n == boost::num_vertices(G) );
      assert( probdata->m == boost::num_edges(G) );

      SCIP_Bool valid = TRUE;
      EdgeIterator eit, eend;
      for (boost::tie(eit, eend) = boost::edges(G); eit != eend; ++eit)
      {
         Vertex s = boost::source(*eit, G);
         Vertex t = boost::target(*eit, G);
         SCIP_Real vals = SCIPgetSolVal(scip, bestsol, probdata->vars[s]);
         SCIP_Real valt = SCIPgetSolVal(scip, bestsol, probdata->vars[t]);

         if ( SCIPisFeasEQ(scip, vals, 1.0) && SCIPisFeasEQ(scip, valt, 1.0) )
         {
            SCIPinfoMessage(scip, nullptr, "Adjacent nodes %lu and %lu are both 1 in best solution.\n", s, t);
            valid = FALSE;
         }
      }

      if ( valid )
         SCIPinfoMessage(scip, nullptr, "Best solution passed check for being a stable set.\n");
      else
      {
         SCIPerrorMessage("Best solution did not pass check for being a stable set.\n");
         SCIPABORT();
         return SCIP_ERROR;
      }
   }

   return SCIP_OKAY;
}

//! check whether best solution of SCIP corresponds to a stable set in the original graph
SCIP_RETCODE BACScheckBestSolOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data structure */
   )
{
   // skip if no graph presolving was performed
   if ( ! probdata->changedpresol )
      return SCIP_OKAY;

   SCIP_SOL* bestsol = SCIPgetBestSol(scip);
   if ( bestsol != nullptr )
   {
      const Graph& Gorig = *probdata->Gorig;
      size_t norig = probdata->norig;
      assert( probdata->norig == boost::num_vertices(Gorig) );
      assert( probdata->morig == boost::num_edges(Gorig) );

      // init original solution
      SCIP_Bool* origsol;
      SCIP_CALL( SCIPallocBufferArray(scip, &origsol, probdata->norig) );

      // transfer the solution to the original graph
      SCIP_Real obj = 0.0;
      for (size_t i = 0; i < norig; ++i)
      {
         if ( probdata->origfixings[i] == BACSfixedzero )
         {
            if ( SCIPdisjointsetFind(probdata->mergemapping, (int) i) != (int) i )
               continue;

            assert( probdata->origtopresol[i] < 0 );
            origsol[i] = FALSE;
         }
         else if ( probdata->origfixings[i] == BACSfixedone )
         {
            assert( probdata->origtopresol[i] < 0 );
            origsol[i] = TRUE;
            SCIP_Real weight = boost::get(vertex_weight_t(), Gorig, (Vertex) i);
            obj += weight;
         }
         else
         {
            if ( SCIPdisjointsetFind(probdata->mergemapping, (int) i) != (int) i )
            {
               assert( SCIPdisjointsetGetComponentCount(probdata->mergemapping) != (int) probdata->norig );
               continue;
            }

            assert( probdata->origfixings[i] == BACSunfixed );
            assert( probdata->origtopresol[i] >= 0 );

            size_t idx = (size_t) probdata->origtopresol[i];
            assert( idx < probdata->n );

            SCIP_Real val = SCIPgetSolVal(scip, bestsol, probdata->vars[idx]);
            assert( SCIPisFeasEQ(scip, val, 0.0) || SCIPisFeasEQ(scip, val, 1.0) );

            // if node was contracted
            if ( probdata->origtopresolparity[i] != 0 )
            {
               if ( probdata->origtopresolparity[i] == -1 )
                  val = 1.0 - val;
            }

            if ( SCIPisFeasZero(scip, val) )
               origsol[i] = FALSE;
            else
            {
               origsol[i] = TRUE;
               SCIP_Real weight = boost::get(vertex_weight_t(), Gorig, (Vertex) i);
               obj += weight;
            }
         }
      }

      if ( SCIPdisjointsetGetComponentCount(probdata->mergemapping) != (int) probdata->norig )
      {
         // process merged nodes
         size_t image;
         for (size_t i = 0; i < norig; ++i)
         {
            assert( SCIPdisjointsetFind(probdata->mergemapping, (int) i) >= 0 );
            image = (size_t) SCIPdisjointsetFind(probdata->mergemapping, (int) i);

            if ( image == i )
               continue;

            origsol[i] = origsol[image];
            if ( origsol[i] )
            {
               SCIP_Real weight = boost::get(vertex_weight_t(), Gorig, (Vertex) i);
               obj += weight;
            }
         }
      }

      // check value
      SCIP_Bool valid = TRUE;
      if ( ! SCIPisSumEQ(scip, obj, SCIPgetSolOrigObj(scip, bestsol)) )
      {
         SCIPinfoMessage(scip, nullptr, "\nOriginal solution has different objective %f as computed solution %f.\n", obj, SCIPgetSolOrigObj(scip, bestsol));
         valid = FALSE;
      }

      // check validity
      EdgeIterator eit, eend;
      for (boost::tie(eit, eend) = boost::edges(Gorig); eit != eend; ++eit)
      {
         Vertex s = boost::source(*eit, Gorig);
         Vertex t = boost::target(*eit, Gorig);

         if ( origsol[s] && origsol[t] )
         {
            SCIPinfoMessage(scip, nullptr, "Adjacent nodes %lu and %lu are both 1 in original solution.\n", s, t);
            valid = FALSE;
         }
      }

      if ( valid )
         SCIPinfoMessage(scip, nullptr, "Original solution passed check for being a stable set.\n");
      else
      {
         SCIPerrorMessage("Original solution did not pass check for being a stable set.\n");
         SCIPABORT();
         return SCIP_ERROR;
      }

      SCIPfreeBufferArray(scip, &origsol);
   }

   return SCIP_OKAY;
}
