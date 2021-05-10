/*----------------------------------------------------------------------------*///
// Created by ledouxf on 12/20/18.
//

/*----------------------------------------------------------------------------*/
#ifndef GMDS_BIT_VECTOR_H
#define GMDS_BIT_VECTOR_H
/*----------------------------------------------------------------------------*/
// STL File headers
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cassert>
/*----------------------------------------------------------------------------*/
#include "CommonTypes.h"
/*----------------------------------------------------------------------------*/
namespace gmds{

    class BitVector
    {
    public:

        typedef int size_type;


        class iterator {
        public:

            friend class BitVector;

            using self_type = iterator;
            using iterator_category = std::forward_iterator_tag;
            using value_type = TInt;
            using difference_type = int;
            using pointer = TInt*;
            using reference = TInt&;

            iterator(BitVector* AContainer):container_(AContainer),current_(0){
                // we value the location of the first real item
                while(current_!=container_->m_top && !container_->m_bits[current_])
                    current_++;
            }

            iterator(BitVector* AContainer, size_type ACurrent):container_(AContainer),current_(ACurrent){}
            iterator(const iterator& AIt):container_(AIt.container_),current_(AIt.current_){}


            self_type operator++() {
                // move to the first real successor of index_
                do{ current_++;
                }while(current_!=container_->m_top && !container_->m_bits[current_]);
            return *this;
            }

            reference operator*()  { return current_; }
            pointer operator->() { return &current_; }
            bool operator==(const self_type& rhs) { return (container_== rhs.container_ && current_==rhs.current_);}
            bool operator!=(const self_type& rhs) { return (container_!= rhs.container_ || current_!=rhs.current_);}
        private:
            BitVector* container_;
            size_type current_;
        };

        class const_iterator
        {
        public:
            using self_type= const_iterator;
            using iterator_category = std::forward_iterator_tag;
            using value_type = TInt;
            using difference_type = int;
            using pointer = TInt*;
            using reference = TInt&;


            const_iterator(const BitVector* AContainer):container_(AContainer),current_(0){
                // we value the location of the first real item
                while(current_!=container_->m_top && !container_->m_bits[current_])
                    current_++;
            }

            const_iterator(const BitVector* AContainer, size_type ACurrent)
                    :container_(AContainer),current_(ACurrent){}
            const_iterator(const const_iterator& AIt)
                    :container_(AIt.container_),current_(AIt.current_){}


            self_type operator++() {
                // move to the first real successor of index_
                do{ current_++;
                }while(current_!=container_->m_top && !container_->m_bits[current_]);
                return *this;
            }
            reference operator*() { return current_; }
            pointer operator->() { return &(current_); }
            bool operator==(const self_type& rhs) { return (container_== rhs.container_ && current_==rhs.current_);}
            bool operator!=(const self_type& rhs) { return (container_!= rhs.container_ || current_!=rhs.current_);}
        private:
            const BitVector* container_;
            size_type current_;
        };



/*------------------------------------------------------------------------*/
/** \brief  Default Constructor
 */
        BitVector(const int capacity = GChunkSize);

/*------------------------------------------------------------------------*/
/** \brief  Copy constructor. Note the algorithm is linear in its size.
 */
        BitVector(const BitVector &AVec);

/*------------------------------------------------------------------------*/
/** \brief  Destructor. Note the algorithm is linear in the number of holes
 * 			in the container.
 */
        ~BitVector();

        /*------------------------------------------------------------------------*/
/** \brief  Gives the number of elements stored in the container.
 */
        size_type size() const;

/*------------------------------------------------------------------------*/
/** \brief  Gives the next right index after the last used item.
 */
        size_type top() const;

/*------------------------------------------------------------------------*/
/** \brief  Gives the size of the container, i.e, the effective number of
 * 			items that are available in the container. It corresponds to
 * 			the number of stored elements + the number of free spaces.
 */
        size_type capacity() const;

/*------------------------------------------------------------------------*/
/** \brief  Indicates if the container is empty.
 */
        bool empty() const;

/*------------------------------------------------------------------------*/
/** \brief  Clear the container and resizes it to 2.
 */
        void clear();
/*------------------------------------------------------------------------*/
/** \brief  resizes the container to ASize. size and top are  both
 * 			equal to ASize. If ASize>current(size) it is just extended,
 * 			new items are marked free. Otherwise, all the items between
 * 			ASize and current(size) are lost.
 */
        void resize(const size_type& ASize);

        /*------------------------------------------------------------------------*/
/** \brief  make the container full of elements putting all the marks to 1
 * 		    and top_ equals to size_;
 */
        void fillAll();

/*------------------------------------------------------------------------*/
/** \brief  This method is necessary when you want to regularize the
 * 			container after several assignements. Indeed assignement method
 * 			has the advantage to check nothing before insertion. The problem
 * 			is then that the container is no more coherent. This method fix
 * 			the container.
 */
        void update();

        /*------------------------------------------------------------------------*/
/** \brief  Give the value stored in the AIndex item but do not allow to
 * 			modify it.
 *
 * 			Warning, in release mode, no test is performed to ensure the
 * 			choice of AIndex.
 */
        bool  operator[](const size_type AIndex) const;


/*------------------------------------------------------------------------*/
/** \brief  Provides the next index where an item can be added. It also
            reserve the necessary space for it and update some internal
            data. It must be used with the assign method.
 */
        size_type  selectNewBit();

        /*------------------------------------------------------------------------*/
        /** \brief  Unselect the bit stored in AIndex.
         *
         * 			Warning - In release mode no test is performed on AIndex. In
         * 			debug mode, AIndex must be in [0, size[
         */
        void unselect(const TInt& AIndex);

        /*------------------------------------------------------------------------*/
        /** \brief  Indicates if the AIndex item is available.
         */
        bool isAvailable(const TInt& AIndex) const;
        /*------------------------------------------------------------------------*/
        /** \brief  Indicates if the AIndex item is out of the container.
         */
        bool isOutOfContainer(const TInt& AIndex) const;
        /*------------------------------------------------------------------------*/
        /** \brief  Provides an iterator on the first element of the container.
         */

        /*------------------------------------------------------------------------*/
        /** \brief Serialize the variable into stream str. Warning this method does
         * 		   not support typename T where pointers would be present.
         */
        void serialize(std::ostream& stream);

        /*------------------------------------------------------------------------*/
        /** \brief Unserialize the variable from stream str. Warning this method
         * 		   does not support typename T where pointers would be present.
         */
        void unserialize(std::istream& stream);

        iterator begin()
        {
            return iterator(this);
        }

        iterator end()
        {
            return iterator(this, m_top);
        }

        const_iterator begin() const
        {
            return const_iterator(this);
        }

        const_iterator end() const
        {
            return const_iterator(this,m_top);
        }


        /*------------------------------------------------------------------------*/
        /** \brief  Provides the next index where an item can be added. It also
        reserve the necessary space for it and update some internal
        data. It must be used with the assign method.
        */
        inline TInt getFreeIndex(){
            TInt free=0;
            TInt cap=m_capacity;

            if(m_free_stack.empty()){
                free= m_top;
                m_top++;
                if(m_top>=cap){
                    //				m_bits.reserve((cap+1)*2);
                    m_bits.resize((cap+1)*2,false);
                    m_capacity = (cap+1)*2;
                }
            }
            else{
                free = m_free_stack.back();
                m_free_stack.pop_back();
            }

            return free;
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Assign the bit AIndex to true. If AIndex is out of the vector
         * 			bounds nothing is specified to avoid the insertion.
         */
        inline void assign(const TInt& AIndex){
            if(m_bits[AIndex]==0) {
                m_size++;
                m_bits[AIndex]=1;
            }
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Return the vector of the boolean mark.
         *
         * \return the vector of the boolean mark
         */
        std::vector<bool> getBits() {
            return m_bits;
        }


///*------------------------------------------------------------------------*/
///** \brief  Compact the container by removing all free items.
// */
//void compact();
///*------------------------------------------------------------------------*/
///** \brief  Compact the container by removing all free items and provide
// * 			the movement of elements in AMove. Elt in location i is moved
// * 			to location AMove[i].
// *
// */
//void compactWithMemory(std::vector<int>& AMove);
//
///*------------------------------------------------------------------------*/
///** \brief Serialize the variable into stream str. Warning this method does
// * 		   not support typename T where pointers would be present.
// */
//inline void serialize(std::ostream& stream);
//
///*------------------------------------------------------------------------*/
///** \brief Unserialize the variable from stream str. Warning this method
// * 		   does not support typename T where pointers would be present.
// */
//inline void unserialize(std::istream& stream);
//
///*------------------------------------------------------------------------*/
///** \brief  Return the vector of the boolean mark.
// *
// * \return the vector of the boolean mark
// */
//std::vector<bool> getMarks() {
//    return mark_;
//}

    protected:
/* top_ is the index of the the right next item after the last used)*/
        size_type m_top;
        size_type m_size;
        size_type m_capacity;
        std::vector<size_type> m_free_stack;
        std::vector<bool> m_bits;
    };
}
#endif //GMDS_SMART_VECTOR_H