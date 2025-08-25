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

/**@file   extiterators.hpp
 * @brief  extended iterators over two graphs
 * @author Marc Pfetsch
 */

#ifndef EXTITERATORS_HPP
#define EXTITERATORS_HPP

#include "graph.h"

namespace bacs
{
   /** extended adjacency iterator */
   struct ExtAdjacencyIterator
   {
      inline ExtAdjacencyIterator() : _G1(nullptr), _G2(nullptr), _first(true)
      {};

      inline ExtAdjacencyIterator(const Graph* G1, const Graph* G2, Vertex v)
      {
         _G1 = G1;
         _G2 = G2;
         boost::tie(_A1, _A1end) = boost::adjacent_vertices(v, *G1);
         if ( G2 != nullptr )
            boost::tie(_A2, _A2end) = boost::adjacent_vertices(v, *G2);
         _first = true;
      }

      bool at_end() const
      {
         if ( _first )
         {
            if ( _A1 == _A1end )
               return true;
         }
         else
         {
            if ( _A2 == _A2end )
               return true;
         }
         return false;
      }

      Vertex
      operator*() const
      {
         if ( _first )
            return *_A1;
         return *_A2;
      }

      ExtAdjacencyIterator&
      operator++() &
      {
         if ( _first )
         {
            ++_A1;
            if ( _A1 == _A1end && _G2 != nullptr )
               _first = false;
         }
         else
            ++_A2;
         return *this;
      }

      ExtAdjacencyIterator
      operator++(int) &&
      {
         ExtAdjacencyIterator _tmp = *this;
         if ( _first )
         {
            ++_A1;
            if ( _A1 == _A1end && _G2 != nullptr )
               _first = false;
         }
         else
            ++_A2;

         return _tmp;
      }

      friend bool
      operator==(const ExtAdjacencyIterator& _X, const ExtAdjacencyIterator& _Y)
      {
         if ( _X._first != _Y._first )
            return false;
         if ( _X._first )
            return ( _X._A1 == _Y._A1 );
         else
            return ( _X._A2 == _Y._A2 );
      }

      friend bool
      operator!=(const ExtAdjacencyIterator& _X, const ExtAdjacencyIterator& _Y)
      {
         if ( _X._first != _Y._first )
            return true;
         if ( _X._first )
            return ( _X._A1 != _Y._A1 );
         else
            return ( _X._A2 != _Y._A2 );
      }

      // member variables
      const Graph*       _G1;                /**< first graph */
      const Graph*       _G2;                /**< second graph */
      AdjacencyIterator  _A1;                /**< iterator for first graph */
      AdjacencyIterator  _A2;                /**< iterator for second graph */
      AdjacencyIterator  _A1end;             /**< end iterator for first graph */
      AdjacencyIterator  _A2end;             /**< end iterator for second graph */
      bool               _first;             /**< wether we are in the first graph */
   };

   inline ExtAdjacencyIterator
   adjacent_vertices(Vertex v, const Graph* G1, const Graph* G2)
   {
      assert( G1 != nullptr );
      return ExtAdjacencyIterator(G1, G2, v);
   }


   /** extended edge iterator */
   struct ExtEdgeIterator
   {
      inline ExtEdgeIterator() : _G1(nullptr), _G2(nullptr), _first(true)
      {};

      inline ExtEdgeIterator(const Graph* G1, const Graph* G2)
      {
         _G1 = G1;
         _G2 = G2;
         boost::tie(_E1, _E1end) = boost::edges(*G1);
         if ( G2 != nullptr )
            boost::tie(_E2, _E2end) = boost::edges(*G2);
         _first = true;
      }

      bool at_end() const
      {
         if ( _first )
         {
            if ( _E1 == _E1end )
               return true;
         }
         else
         {
            if ( _E2 == _E2end )
               return true;
         }
         return false;
      }

      Edge
      operator*() const
      {
         if ( _first )
            return *_E1;
         return *_E2;
      }

      ExtEdgeIterator&
      operator++() &
      {
         if ( _first )
         {
            ++_E1;
            if ( _E1 == _E1end && _G2 != nullptr )
               _first = false;
         }
         else
            ++_E2;
         return *this;
      }

      ExtEdgeIterator
      operator++(int) &&
      {
         ExtEdgeIterator _tmp = *this;
         if ( _first )
         {
            ++_E1;
            if ( _E1 == _E1end && _G2 != nullptr )
               _first = false;
         }
         else
            ++_E2;

         return _tmp;
      }

      Vertex source() const
      {
         if ( _first )
            return boost::source(*_E1, *_G1);
         assert( _G2 != nullptr );
         return boost::source(*_E2, *_G2);
      }

      Vertex target() const
      {
         if ( _first )
            return boost::target(*_E1, *_G1);
         assert( _G2 != nullptr );
         return boost::target(*_E2, *_G2);
      }

      friend bool
      operator==(const ExtEdgeIterator& _X, const ExtEdgeIterator& _Y)
      {
         if ( _X._first != _Y._first )
            return false;
         if ( _X._first )
            return ( _X._E1 == _Y._E1 );
         else
            return ( _X._E2 == _Y._E2 );
      }

      friend bool
      operator!=(const ExtEdgeIterator& _X, const ExtEdgeIterator& _Y)
      {
         if ( _X._first != _Y._first )
            return true;
         if ( _X._first )
            return ( _X._E1 != _Y._E1 );
         else
            return ( _X._E2 != _Y._E2 );
      }

      // member variables
      const Graph*       _G1;                /**< first graph */
      const Graph*       _G2;                /**< second graph */
      EdgeIterator       _E1;                /**< iterator for first graph */
      EdgeIterator       _E2;                /**< iterator for second graph */
      EdgeIterator       _E1end;             /**< end iterator for first graph */
      EdgeIterator       _E2end;             /**< end iterator for second graph */
      bool               _first;             /**< wether we are in the first graph */
   };

   inline ExtEdgeIterator
   edges(const Graph* G1, const Graph* G2)
   {
      assert( G1 != nullptr );
      return ExtEdgeIterator(G1, G2);
   }
}

#endif
