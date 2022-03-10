/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
 * NContainer.cpp
 *
 *  Created on: 03/08/2017
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <KM/DS/Mesh.h>
#include <KM/DS/NContainer.h>
/*----------------------------------------------------------------------------*/
using namespace kmds;
/*----------------------------------------------------------------------------*/
NContainer::NContainer(Mesh* AMesh, const TInt32 ACapacity)
 : m_mesh(AMesh)
 , m_cells("N_ID", ACapacity)
 , m_xyz("N_XYZ", ACapacity)
 , m_top(0)
{
        Kokkos::parallel_for(ACapacity, KOKKOS_LAMBDA(const TSize i) { m_cells(i) = false; });
}
/*----------------------------------------------------------------------------*/
NContainer::~NContainer()
{
}
/*----------------------------------------------------------------------------*/
bool
NContainer::has(const TCellID AID) const
{
        return m_cells(AID);
}
/*----------------------------------------------------------------------------*/
TCellID
NContainer::add()
{
        TSize i = Kokkos::atomic_fetch_add(&m_top, 1);
        m_cells(i) = true;
        return i;
}
/*----------------------------------------------------------------------------*/
TCellID
NContainer::add(const TInt32 ANb)
{
        TSize x = Kokkos::atomic_fetch_add(&m_top, ANb);
        for (auto i = x; i < x + ANb; i++) {
                m_cells(i) = true;
        }
        return x;
}
/*----------------------------------------------------------------------------*/
void
NContainer::remove(const TCellID AId)
{
        m_cells(AId) = false;
}
/*----------------------------------------------------------------------------*/
TSize
NContainer::capacity() const
{
        return m_cells.extent(0);
}
/*----------------------------------------------------------------------------*/
TSize
NContainer::top() const
{
        return m_top;
}
/*----------------------------------------------------------------------------*/
TSize
NContainer::nbCells() const
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
NContainer::set(const TCellID AId, const TCoord AX, const TCoord AY, const TCoord AZ)
{
        m_xyz(AId, 0) = AX;
        m_xyz(AId, 1) = AY;
        m_xyz(AId, 2) = AZ;
}
/*----------------------------------------------------------------------------*/
TCellID
NContainer::newNode(const TCoord AX, const TCoord AY, const TCoord AZ)
{
        TCellID i = add();
        m_xyz(i, 0) = AX;
        m_xyz(i, 1) = AY;
        m_xyz(i, 2) = AZ;
        return i;
}
/*----------------------------------------------------------------------------*/
void
NContainer::get(const TCellID AId, TCoord& AX, TCoord& AY, TCoord& AZ) const
{
        AX = m_xyz(AId, 0);
        AY = m_xyz(AId, 1);
        AZ = m_xyz(AId, 2);
}
/*----------------------------------------------------------------------------*/
void
NContainer::resize(const TInt32 ASize)
{
        TSize prev_capacity = capacity();

        if(prev_capacity < ASize) {

                Kokkos::resize(m_cells, ASize);
                Kokkos::resize(m_xyz, ASize);

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
NContainer::removeAll()
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
