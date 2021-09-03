/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
 * FContainer.cpp
 *
 *  Created on: 03/08/2017
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <map>
/*----------------------------------------------------------------------------*/
#include <KM/DS/FContainer.h>
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
using namespace kmds;
/*----------------------------------------------------------------------------*/
FContainer::FContainer(Mesh* AMesh, const TInt32 ACapacity)
 : m_mesh(AMesh)
 , m_first("F_FIRST",  ACapacity)
 , m_nb_nodes("F_NB_NODES", ACapacity)
 , m_F2N("F_F2N", 5 * (TSize) ACapacity)
 , m_top(0)
 , m_top_F2N(0)
{
        Kokkos::parallel_for(ACapacity, KOKKOS_LAMBDA(const TSize i) { m_first(i) = NullTSize; });
}
/*----------------------------------------------------------------------------*/
FContainer::~FContainer()
{
}
/*----------------------------------------------------------------------------*/
bool
FContainer::has(const TCellID AID) const
{
        return m_first(AID) != NullTSize;
}
/*----------------------------------------------------------------------------*/
TCellID
FContainer::addTriangle()
{
        TSize i = Kokkos::atomic_fetch_add(&m_top, 1);
        TSize first_node = Kokkos::atomic_fetch_add(&m_top_F2N, 3);
        m_first(i) = first_node;
        m_nb_nodes(i) = 3;
        return i;
} /*----------------------------------------------------------------------------*/
TCellID
FContainer::addQuad()
{
        TSize i = Kokkos::atomic_fetch_add(&m_top, 1);
        TSize first_node = Kokkos::atomic_fetch_add(&m_top_F2N, 4);
        m_first(i) = first_node;
        m_nb_nodes(i) = 4;
        return i;
} /*----------------------------------------------------------------------------*/
TCellID
FContainer::addPentagon()
{
        TSize i = Kokkos::atomic_fetch_add(&m_top, 1);
        TSize first_node = Kokkos::atomic_fetch_add(&m_top_F2N, 5);
        m_first(i) = first_node;
        m_nb_nodes(i) = 5;
        return i;
}
/*----------------------------------------------------------------------------*/
TCellID
FContainer::addTriangles(const TInt32 ANb)
{
        TSize x = Kokkos::atomic_fetch_add(&m_top, ANb);
        TSize y = Kokkos::atomic_fetch_add(&m_top_F2N, 3 * (TSize) ANb);
        for (auto i = x; i < x + ANb; i++) {
                m_first(i) = y;
                m_nb_nodes(i) = 3;
                y = y + 3;
        }
        return x;
}
/*----------------------------------------------------------------------------*/
TCellID
FContainer::addQuads(const TInt32 ANb)
{
        TSize x = Kokkos::atomic_fetch_add(&m_top, ANb);
        TSize y = Kokkos::atomic_fetch_add(&m_top_F2N, 4 * (TSize) ANb);
        for (auto i = x; i < x + ANb; i++) {
                m_first(i) = y;
                m_nb_nodes(i) = 4;
                y = y + 4;
        }
        return x;
}
/*----------------------------------------------------------------------------*/
TCellID
FContainer::addPentagons(const TInt32 ANb)
{
        TSize x = Kokkos::atomic_fetch_add(&m_top, ANb);
        TSize y = Kokkos::atomic_fetch_add(&m_top_F2N, 5 * (TSize) ANb);
        for (auto i = x; i < x + ANb; i++) {
                m_first(i) = y;
                m_nb_nodes(i) = 5;
                y = y + 5;
        }
        return x;
}
/*----------------------------------------------------------------------------*/
TCellID
FContainer::newTriangle(const TCellID AN1, const TCellID AN2, const TCellID AN3)
{
        TSize i = addTriangle();
        TSize first = m_first(i);
        m_F2N(first) = AN1;
        m_F2N(first + 1) = AN2;
        m_F2N(first + 2) = AN3;
        return i;
}
/*----------------------------------------------------------------------------*/
TCellID
FContainer::newQuad(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4)
{
        TSize i = addQuad();
        TSize first = m_first(i);
        m_F2N(first) = AN1;
        m_F2N(first + 1) = AN2;
        m_F2N(first + 2) = AN3;
        m_F2N(first + 3) = AN4;
        return i;
}
/*----------------------------------------------------------------------------*/
TCellID
FContainer::newPentagon(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4, const TCellID AN5)
{
        TSize i = addPentagon();
        TSize first = m_first(i);
        m_F2N(first) = AN1;
        m_F2N(first + 1) = AN2;
        m_F2N(first + 2) = AN3;
        m_F2N(first + 3) = AN4;
        m_F2N(first + 4) = AN5;
        return i;
}
/*----------------------------------------------------------------------------*/
void
FContainer::remove(const TCellID AId)
{
        m_first(AId) = NullTSize;
}
/*----------------------------------------------------------------------------*/
void
FContainer::get(const TCellID AI, Kokkos::View<TCellID*>& AN) const
{
        TSize first = m_first(AI);
        int nb = m_nb_nodes(AI);
        AN = Kokkos::subview(m_F2N, std::make_pair(first, first + nb));
}
/*----------------------------------------------------------------------------*/
void
FContainer::get(const TCellID AI, TCellID AIds[5], int* ASize) const
{
        TSize first = m_first(AI);
        int nb = m_nb_nodes(AI);

        *ASize = nb;
        for(int i=0; i<nb; i++) {
                AIds[i] = m_F2N(first+i);
        }
}
/*----------------------------------------------------------------------------*/
TInt32
FContainer::getNbNodes(const TCellID AI) const
{
        return m_nb_nodes(AI);
}
/*----------------------------------------------------------------------------*/
void
FContainer::set(const TCellID AI, const Kokkos::View<TCellID*>& AN)
{
        TSize first = m_first(AI);
        int nb = m_nb_nodes(AI);
        if (AN.size() <= nb) {
                setUnsafe(AI, AN);
        } else {
                // We don't have enough space to put the new node ids
                // face id doesn't change but first and nb must point on new
                // memory locations in F2N
                // it means so that we let some new holes in the m_F2N FContainer
                nb = AN.size();
                // we get a new first item in F2N
                first = Kokkos::atomic_fetch_add(&m_top_F2N, nb);
                // we update data of face AI now
                m_first(AI) = first;
                m_nb_nodes(AI) = nb;
                for (auto i = 0; i < nb; i++) {
                        m_F2N(first + i) = AN(i);
                }
        }
}
/*----------------------------------------------------------------------------*/
void
FContainer::setUnsafe(const TCellID AI, const Kokkos::View<TCellID*>& AN)
{
        TSize first = m_first(AI);
        int nb = AN.size();
        // this allocation is useful when we replace a quad by a triangle
        // for instance
        m_nb_nodes(AI) = nb;
        for (auto i = 0; i < nb; i++) {
                m_F2N(first + i) = AN(i);
        }
}
/*----------------------------------------------------------------------------*/
void
FContainer::setUnsafe(const TCellID AI, TCellID* AIds, int ASize)
{
        TSize first = m_first(AI);
        int nb = ASize;
        // this allocation is useful when we replace a quad by a triangle
        // for instance
        m_nb_nodes(AI) = nb;
        for (auto i = 0; i < nb; i++) {
                m_F2N(first + i) = AIds[i];
        }
}
/*----------------------------------------------------------------------------*/
TSize
FContainer::capacity() const
{
        return m_first.extent(0);
}
/*----------------------------------------------------------------------------*/
TSize
FContainer::top() const
{
        return m_top;
}
/*----------------------------------------------------------------------------*/
TSize
FContainer::nbCells() const
{
        TSize n = m_top;
        TSize s = 0;
        Kokkos::parallel_reduce(n,
                                KOKKOS_LAMBDA(const TSize i, TSize& ls) {
                                        if (this->has(i)) {
                                                ls = ls + 1;
                                        }
                                },
                                s);
        return s;
}
/*----------------------------------------------------------------------------*/
TSize
FContainer::nbTriangles() const
{
        TSize n = m_top;
        TSize s = 0;
        Kokkos::parallel_reduce(n,
                                KOKKOS_LAMBDA(const TSize i, TSize& ls) {
                                        if (this->has(i) && m_nb_nodes(i) == 3) {
                                                ls = ls + 1;
                                        }
                                },
                                s);
        return s;
}
/*----------------------------------------------------------------------------*/
TSize
FContainer::nbQuads() const
{
        TSize n = m_top;
        TSize s = 0;
        Kokkos::parallel_reduce(n,
                                KOKKOS_LAMBDA(const TSize i, TSize& ls) {
                                        if (this->has(i) && m_nb_nodes(i) == 4) {
                                                ls = ls + 1;
                                        }
                                },
                                s);
        return s;
}
/*----------------------------------------------------------------------------*/
TSize
FContainer::nbPentagons() const
{
        TSize n = m_top;
        TSize s = 0;
        Kokkos::parallel_reduce(n,
                                KOKKOS_LAMBDA(const TSize i, TSize& ls) {
                                        if (this->has(i) && m_nb_nodes(i) == 5) {
                                                ls = ls + 1;
                                        }
                                },
                                s);
        return s;
}
/*----------------------------------------------------------------------------*/
void
FContainer::resize(const TInt32 ASize)
{
        TSize prev_capacity = capacity();

        if(prev_capacity < ASize) {

                Kokkos::resize(m_first, ASize);
                Kokkos::resize(m_nb_nodes, ASize);
                Kokkos::resize(m_F2N, (TSize) ASize * 5);

                // We must assign NULLTSize to new empty items
                TSize from = prev_capacity;
                TSize nb = ASize - prev_capacity;
                Kokkos::View<TSize *> slice = Kokkos::subview(m_first, std::make_pair(from, prev_capacity + nb));
                Kokkos::parallel_for(nb, KOKKOS_LAMBDA(const TSize i) {
                    slice(i) = NullTSize;
                });
        }
}
/*----------------------------------------------------------------------------*/
void
FContainer::removeAll()
{

    TSize n = m_top;
    m_top = 0;
    Kokkos::parallel_for(n,
                         KOKKOS_LAMBDA(const TSize i) {
                             if (this->has(i)) {
                                 this->remove(i);
                             }
                         });
}
/*----------------------------------------------------------------------------*/
