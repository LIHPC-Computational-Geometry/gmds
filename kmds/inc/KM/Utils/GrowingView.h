/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    GrowingView.h
 *  \author  F. LEDOUX
 *  \date    03/13/2017
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_GROWING_VIEW_H_
#define KMDS_GROWING_VIEW_H_
/*----------------------------------------------------------------------------*/
// Kokkos File Headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// KMDS File Headers
#include "KTypes.h"
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
/** \class GrowingView
 *  \brief Provides a view with atomic top counter to make it grows. The param
 *         T must fit Kokkos::View requirements
 */
/*----------------------------------------------------------------------------*/
template <typename T>
class GrowingView
{
 public:
        /*------------------------------------------------------------------------*/
        /** \brief Default constructor.
         *
         *  \param[in] AName view name for Kokkos debug purpose
         *  \param[in] ASize max capacity of this view
         */
        GrowingView(const std::string AName, const int ASize)
         : m_view(AName, ASize)
         , m_capacity(ASize)
         , m_top(0)
        {
                ;
        }

        /*------------------------------------------------------------------------*/
        /** \brief Getter on the element stored at index \p AIndex
         *
         *  \param[in] AIndex element index we want to access to
         */
        T
        get(const TInt32 AIndex) const
        {
                return m_view(AIndex);
        }

        /*------------------------------------------------------------------------*/
        /** \brief Setter on the element stored at index \p AIndex
         *
         *  \param[in] AIndex element index we want to change the value of
         *  \param[in] AValue the value to assign
         */
        void
        set(const TInt32 AIndex, const T AValue)
        {
                m_view(AIndex) = AValue;
        }

        /*------------------------------------------------------------------------*/
        /** \brief Add an element at the top of the container. We assume that the
         *         container has enough space to store \p AElem
         *
         *  \param[in] AElem  the element to add in the container
         */
        TInt32
        push_back(const T AElem)
        {
                TInt32 i = Kokkos::atomic_fetch_add(&m_top, 1);
                m_view(i) = AElem;
                return i;
        }

        /*------------------------------------------------------------------------*/
        /** \brief Add several elements at the top of the container. We assume that the
         *         container has enough space to store \p ANbElem more
         *
         *  \param[in] ANbElem  the number of elements to add in the container
         *
         *  \return the index of the first element of this batch
         */
        TInt32
        addElems(const TInt32 ANbElem)
        {
                TInt32 i = Kokkos::atomic_fetch_add(&m_top, ANbElem);
                return i;
        }

        /*------------------------------------------------------------------------*/
        /** \brief Returs the number of element that are stored
         */
        TInt32
        getNbElems() const
        {
                return m_top;
        }
        /*------------------------------------------------------------------------*/
        /** \brief Returns the capacity of the container
         */
        TInt32
        capacity() const
        {
                return m_capacity;
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Resize the container
         */
        void
        resize(const TInt32 ASize)
        {
                Kokkos::resize(m_view, ASize);
                m_capacity = ASize;
                m_top = ASize;
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Reserve space for the container
         */
        void
        reserve(const TInt32 ASize)
        {
                Kokkos::resize(m_view, ASize);
                m_capacity = ASize;
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Clear the container
         */
        void
        clear()
        {
            m_top = 0;
        }

        /*------------------------------------------------------------------------*/
        /** \brief  Set the top of the container
         */
        void
        setTop(const TInt32 ASize)
        {
            m_top = ASize;
        }


private:
        Kokkos::View<T*> m_view;
        TInt32 m_capacity;
        TInt32 m_top;
};

/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_EXCEPTION_H_ */
/*----------------------------------------------------------------------------*/
