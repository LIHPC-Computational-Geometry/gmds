/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/3/19.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_INDEXEDVECTOR_H
#define GMDS_INDEXEDVECTOR_H
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
/* \class IndexedVector
 *
 * \brief Provide a container storing a collection of T objects.
 */
/*----------------------------------------------------------------------------*/
    template<typename T>
    class IndexedVector {

    public:
        /*------------------------------------------------------------------------*/
        /** \brief  Default Constructor
         */
        IndexedVector(const int capacity = GChunkSize) {
            m_vec.resize(capacity);
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Copy constructor. Note the algorithm is linear in its size.
         */
        IndexedVector(const IndexedVector<T> &AVec)
                : m_vec(AVec.m_vec) {}

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor. Note the algorithm is linear in the number of holes
         * 			in the container.
         */
        ~IndexedVector() { m_vec.clear(); }
        /*------------------------------------------------------------------------*/
        /** \brief  Gives the size of the container, i.e, the effective number of
         * 			items that are available in the container. It corresponds to
         * 			the number of stored elements + the number of free spaces.
         */
        TInt capacity() const { return m_vec.size(); }
        /*------------------------------------------------------------------------*/
        /** \brief  Clear the container and resizes it to 2.
         */
        void clear() {
            m_vec.clear();
            m_vec.resize(GChunkSize);
        }

        /*------------------------------------------------------------------------*/
        /** \brief  resizes the container to ASize. size and top are  both
         * 			equal to ASize. If ASize>current(size) it is just extended,
         * 			new items are marked free. Otherwise, all the items between
         * 			ASize and current(size) are lost.
         */
        void resize(const TInt &ASize) {
            m_vec.resize(ASize);
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Give the value stored in the AIndex item but do not allow to
         * 			modify it.
         *
         * 			Warning, in release mode, no test is performed to ensure the
         * 			choice of AIndex.
         */
        T const &operator[](const TInt &AIndex) const {
#ifdef __DEBUG__
            if(AIndex>=capacity())
            throw GMDSException("Bad index in IndexedVector<T>::operator[]");
#endif //__DEBUG__

            return m_vec[AIndex];
        }
        /*------------------------------------------------------------------------*/
        /** \brief  Give the value stored in the AIndex item and allow to
         * 			modify it.
         *
         * 			Warning, in release mode, no test is performed to ensure the
         * 			choice of AIndex. Moreover, in order to improve performances,
         * 			the container structure is partially corrupted. Thus an update
         * 			is necessary at the end (especially for the iterators)
         */
        T &operator[](const TInt &AIndex) {
#ifdef __DEBUG__
            if(AIndex>=capacity())
            throw GMDSException("Bad index in IndexedVector<T>::operator[]");
#endif //__DEBUG__
            return m_vec[AIndex];
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Assign AElt to index AIndex. If this index is still occupied,
         * 			it contents is replaced. If AIndex is out of the vector
         * 			limits, nothing is specified to avoid the insertion.
         */
        void assign(T &AElt, const TInt &AIndex) {
            if ((unsigned int)AIndex >= m_vec.size()) {
                TInt prev_size = m_vec.size();
                m_vec.resize(prev_size * 2);
            }
            m_vec[AIndex] = AElt;
        }


        /*------------------------------------------------------------------------*/
        /** \brief Serialize
         */
        void serialize(std::ostream &stream) {
            int container_size = m_vec.size();
            stream.write((char *) &container_size, sizeof(int));

            for (int i = 0; i < container_size; i++) {
                T elt = m_vec[i];
                stream.write((char *) &elt, sizeof(T));
            }
        }

        /*------------------------------------------------------------------------*/
        /** \brief Unserialize
         */
        void unserialize(std::istream &stream) {
            int nb_items = 0;

            stream.read((char *) &nb_items, sizeof(int));
            clear();
            resize(nb_items);
            for (int i = 0; i < nb_items; i++) {
                T elt;
                stream.read((char *) &elt, sizeof(T));
                m_vec[i] = elt;
            }
        }

    protected:
        std::vector<T> m_vec;
    };
}
/*----------------------------------------------------------------------------*/

#endif //GMDS_INDEXEDVECTOR_H
