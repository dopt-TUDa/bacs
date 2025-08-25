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

/**@file   vectorpool.hpp
 * @brief  pool that stores vectors of unsigned ints
 * @author Marc Pfetsch
 *
 * The downside of this data structure is that it is outside of the memory functions of SCIP, so its storage
 * requirements are not tracked.
 */

#ifndef VECTORPOOL_HPP
#define VECTORPOOL_HPP

#include <vector>
#include <unordered_set>

//! hash function for vectors
struct vectorPtrHashFunction
{
   //!< compute hash function for set stored in @a V
   std::size_t operator()(
      const std::vector<unsigned int>* V     //!< vector
      ) const
   {
      std::size_t sum = 0;
      std::vector<unsigned int>::const_iterator vit, vend = V->end();
      for (vit = V->begin(); vit != vend; ++vit)
         sum += *vit;
      return sum;
   }
};


//! function for comparing two vectors
struct vectorPtrHashEqual
{
   //!< check whether sets stored in @a b1 and @a b2 are equal
   bool operator()(
      const std::vector<unsigned int>* b1,   //!< first set
      const std::vector<unsigned int>* b2    //!< second set
      ) const
   {
      return *b1 == *b2;
   }
};


// Definition of the hash set that stores vectors
typedef std::unordered_set<const std::vector<unsigned int> *, vectorPtrHashFunction, vectorPtrHashEqual> vectorPtrHash;




/**@class vectorpool
 * @brief Provides a storage for storing unsigned int vectors with hash containment testing
 */
class vectorpool
{
public:
   typedef std::vector<unsigned int>* value_type;
   typedef const std::vector<unsigned int>* const_reference;
   typedef std::vector<unsigned int>* reference;
   typedef std::vector<value_type>::const_iterator const_iterator;
   typedef std::vector<value_type>::iterator iterator;
   typedef size_t size_type;

   //! initialize empty pool
   vectorpool(): H_(10000)
   {
   }

   //! initialize empty pool and reserve space for @a n vectors
   vectorpool(
      size_type          n,                  //!< size of space to reserve
      size_type          s                   //!< size of hash
      ): H_(s)
   {
      V_.reserve(n);
   }

   //! copy constructor
   vectorpool(
      const vectorpool&  V                   //!< pool to be copied
      ) : H_(V.hash_size())
   {
      const_iterator it, iend = V.end();
      for (it = V.begin(); it != iend; ++it)
      {
         std::vector<unsigned int>* v = new std::vector<unsigned int>(**it);
         V_.push_back(v);
         (void) H_.insert(v);
      }  /*lint !e429*/
   }

   //! destructor
   ~vectorpool()
   {
      const_iterator it, iend = V_.end();
      for (it = V_.begin(); it != iend; ++it)
      {
         delete *it;
      }
   }

   //! assignment operator
   vectorpool& operator=(
      const vectorpool&  V                   //!< pool to be assigned
      ) & /*lint --e{1539}*/
   {
      if ( this == std::addressof(V) )
         return *this;

      V_.clear();
      H_.clear();
      const_iterator it, iend = V.end();
      for (it = V.begin(); it != iend; ++it)
      {
         std::vector<unsigned int>* v = new std::vector<unsigned int>(**it);
         V_.push_back(v);
         (void) H_.insert(v);
      }  /*lint !e429*/
      return *this;
   }

   //! number of vectors stored in pool
   size_type size() const
   {
      return V_.size();
   }

   //! size of hash table
   size_type hash_size() const
   {
      return H_.size();
   }

   //! check whether v is present in pool
   bool contains(
      const std::vector<unsigned int>& v     //!< set to be checked
      ) const
   {
      return H_.find(&v) != H_.end();
   }

   //! insert @a v - will be freed by this class
   bool insert(
      const std::vector<unsigned int>& v     //!< set to be inserted
      )
   {
      // insert v if not already present
      if ( H_.find(&v) == H_.end() )
      {
         std::vector<unsigned int>* vv = new std::vector<unsigned int>(v);
         V_.push_back(vv);
         (void) H_.insert(vv);
         return true;
      }  /*lint !e429*/
      return false;
   }

   //! return begin iterator
   iterator begin()
   {
      return V_.begin();
   }

   //! return end iterator
   iterator end()
   {
      return V_.end();
   }

   //! return const begin iterator
   const_iterator begin() const
   {
      return V_.begin();
   }

   //! return const end iterator
   const_iterator end() const
   {
      return V_.end();
   }

   //! return vector at position @a i
   reference
   operator[](
      size_type          i                   //!< position in pool
      )
   {
      return V_[i];
   }

   //! return vector at position @a i
   const_reference
   operator[](
      size_type          i                   //!< position in pool
      ) const
   {
      return V_[i];
   }

private:
   std::vector<std::vector<unsigned int>* > V_;   //!< array to store sets stored in vectors
   vectorPtrHash H_;                              //!< hash map
};

#endif
