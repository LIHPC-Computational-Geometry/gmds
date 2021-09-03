/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
 * EContainer.cpp
 *
 *  Created on: 14/02/2018
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include "KM/DS/Mesh.h"
#include "KM/DS/EContainer.h"
/*----------------------------------------------------------------------------*/
using namespace kmds;
/*----------------------------------------------------------------------------*/
EContainer::EContainer(Mesh* AMesh, const TInt32 ACapacity)
        : m_mesh(AMesh)
        , m_cells("N_ID", ACapacity)
        , m_E2N("N_E2N", ACapacity)
        , m_top(0)
{
    Kokkos::parallel_for(ACapacity, KOKKOS_LAMBDA(const TSize i) { m_cells(i) = false; });
}
/*----------------------------------------------------------------------------*/
EContainer::~EContainer()
{
}
/*----------------------------------------------------------------------------*/
bool
EContainer::has(const TCellID AID) const
{
    return m_cells(AID);
}
/*----------------------------------------------------------------------------*/
TCellID
EContainer::add()
{
    TSize i = Kokkos::atomic_fetch_add(&m_top, 1);
    m_cells(i) = true;
    return i;
}
/*----------------------------------------------------------------------------*/
TCellID
EContainer::add_unsafe()
{
    TSize i = m_top;
    m_top++;
    m_cells(i) = true;
    return i;
}
/*----------------------------------------------------------------------------*/
TCellID
EContainer::add(const TInt32 ANb)
{
    TSize x = Kokkos::atomic_fetch_add(&m_top, ANb);
    for (auto i = x; i < x + ANb; i++) {
        m_cells(i) = true;
    }
    return x;
}
/*----------------------------------------------------------------------------*/
void
EContainer::remove(const TCellID AId)
{
    m_cells(AId) = false;
}
/*----------------------------------------------------------------------------*/
TSize
EContainer::capacity() const
{
    return m_cells.extent(0);
}
/*----------------------------------------------------------------------------*/
TSize
EContainer::top() const
{
    return m_top;
}
/*----------------------------------------------------------------------------*/
TSize
EContainer::nbCells() const
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
void
EContainer::set(const TCellID AId, const TCellID AN1, const TCellID AN2)
{
    m_E2N(AId, 0) = AN1;
    m_E2N(AId, 1) = AN2;
}
/*----------------------------------------------------------------------------*/
TCellID
EContainer::newEdge(const TCellID AN1, const TCellID AN2)
{
    TCellID i = add();
    m_E2N(i, 0) = AN1;
    m_E2N(i, 1) = AN2;
    return i;
}
/*----------------------------------------------------------------------------*/
TCellID
EContainer::newEdge_unsafe(const TCellID AN1, const TCellID AN2)
{
    TCellID i = add_unsafe();
    m_E2N(i, 0) = AN1;
    m_E2N(i, 1) = AN2;
    return i;
}
/*----------------------------------------------------------------------------*/
void
EContainer::get(const TCellID AI, Kokkos::View<TCellID*>& AN) const
{
    AN = Kokkos::subview(m_E2N, AI, Kokkos::ALL());
}
/*----------------------------------------------------------------------------*/
void
EContainer::get(const TCellID AId, TCellID &AN1, TCellID &AN2) const
{
    AN1 = m_E2N(AId, 0);
    AN2 = m_E2N(AId, 1);
}
/*----------------------------------------------------------------------------*/
void
EContainer::resize(const TInt32 ASize)
{
    TSize prev_capacity = capacity();

    if(prev_capacity < ASize) {

        Kokkos::resize(m_cells, ASize);
        Kokkos::resize(m_E2N, ASize);

        // We must assign NULLId to new empty items
        TSize from = prev_capacity;
        TSize nb = ASize - prev_capacity;
        Kokkos::View<bool *> slice = Kokkos::subview(m_cells, std::make_pair(from, prev_capacity + nb));
        Kokkos::parallel_for(nb, KOKKOS_LAMBDA(const TSize i) {
            slice(i) = false;
        });

    }
}
/*----------------------------------------------------------------------------*/
void
EContainer::removeAll()
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
