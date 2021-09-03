/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
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
