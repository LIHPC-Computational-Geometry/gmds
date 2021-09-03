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
 * RContainer.cpp
 *
 *  Created on: 06/22/2017
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include <map>
/*----------------------------------------------------------------------------*/
#include <KM/DS/RContainer.h>
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
using namespace kmds;
/*----------------------------------------------------------------------------*/
RContainer::RContainer(Mesh* AMesh, const TInt32 ACapacity)
 : m_mesh(AMesh)
 , m_first("R_FIRST", ACapacity)
 , m_nb_nodes("R_NB_NODES", ACapacity)
 , m_R2N("R_R2N", 8 * (TSize) ACapacity)
 , m_top(0)
 , m_top_R2N(0)
{
        Kokkos::parallel_for(ACapacity, KOKKOS_LAMBDA(const TSize i) { m_first(i) = NullTSize; });
}
/*----------------------------------------------------------------------------*/
RContainer::~RContainer()
{
}
/*----------------------------------------------------------------------------*/
bool
RContainer::has(const TCellID AID) const
{
        return m_first(AID) != NullTSize;
}
/*----------------------------------------------------------------------------*/
TCellID
RContainer::addTetrahedron()
{
        TSize i = Kokkos::atomic_fetch_add(&m_top, 1);
        TSize first_node = Kokkos::atomic_fetch_add(&m_top_R2N, 4);
        m_first(i) = first_node;
        m_nb_nodes(i) = 4;
        return i;
}
/*----------------------------------------------------------------------------*/
TCellID
RContainer::addHexahedron()
{
        TSize i = Kokkos::atomic_fetch_add(&m_top, 1);
        TSize first_node = Kokkos::atomic_fetch_add(&m_top_R2N, 8);
        m_first(i) = first_node;
        m_nb_nodes(i) = 8;
        return i;
}
/*----------------------------------------------------------------------------*/
TCellID
RContainer::addPyramid()
{
        TSize i = Kokkos::atomic_fetch_add(&m_top, 1);
        TSize first_node = Kokkos::atomic_fetch_add(&m_top_R2N, 5);
        m_first(i) = first_node;
        m_nb_nodes(i) = 5;
        return i;
}
/*----------------------------------------------------------------------------*/
TCellID
RContainer::addPrism3()
{
        TSize i = Kokkos::atomic_fetch_add(&m_top, 1);
        TSize first_node = Kokkos::atomic_fetch_add(&m_top_R2N, 6);
        m_first(i) = first_node;
        m_nb_nodes(i) = 6;
        return i;
}
/*----------------------------------------------------------------------------*/
TCellID
RContainer::addTetrahedra(const TInt32 ANb)
{
        TSize x = Kokkos::atomic_fetch_add(&m_top, ANb);
        TSize y = Kokkos::atomic_fetch_add(&m_top_R2N, 4 * (TSize) ANb);
        for (auto i = x; i < x + ANb; i++) {
                m_first(i) = y;
                m_nb_nodes(i) = 4;
                y = y + 4;
        }
        return x;
}
/*----------------------------------------------------------------------------*/
TCellID
RContainer::addHexahedra(const TInt32 ANb)
{
        TSize x = Kokkos::atomic_fetch_add(&m_top, ANb);
        TSize y = Kokkos::atomic_fetch_add(&m_top_R2N, 8 * (TSize) ANb);
        for (auto i = x; i < x + ANb; i++) {
                m_first(i) = y;
                m_nb_nodes(i) = 8;
                y = y + 8;
        }
        return x;
}
/*----------------------------------------------------------------------------*/
TCellID
RContainer::addPyramids(const TInt32 ANb)
{
        TSize x = Kokkos::atomic_fetch_add(&m_top, ANb);
        TSize y = Kokkos::atomic_fetch_add(&m_top_R2N, 5 * (TSize) ANb);
        for (auto i = x; i < x + ANb; i++) {
                m_first(i) = y;
                m_nb_nodes(i) = 5;
                y = y + 5;
        }
        return x;
}
/*----------------------------------------------------------------------------*/
TCellID
RContainer::addPrism3s(const TInt32 ANb)
{
        TSize x = Kokkos::atomic_fetch_add(&m_top, ANb);
        TSize y = Kokkos::atomic_fetch_add(&m_top_R2N, 6 * (TSize) ANb);
        for (auto i = x; i < x + ANb; i++) {
                m_first(i) = y;
                m_nb_nodes(i) = 6;
                y = y + 6;
        }
        return x;
}
/*----------------------------------------------------------------------------*/
TCellID
RContainer::newTetrahedron(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4)
{
        TSize i = addTetrahedron();
        TSize first = m_first(i);
        m_R2N(first) = AN1;
        m_R2N(first + 1) = AN2;
        m_R2N(first + 2) = AN3;
        m_R2N(first + 3) = AN4;
        return i;
}
/*----------------------------------------------------------------------------*/
TCellID
RContainer::newHexahedron(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                        const TCellID AN5, const TCellID AN6, const TCellID AN7, const TCellID AN8)
{
        TSize i = addHexahedron();
        TSize first = m_first(i);
        m_R2N(first) = AN1;
        m_R2N(first + 1) = AN2;
        m_R2N(first + 2) = AN3;
        m_R2N(first + 3) = AN4;
        m_R2N(first + 4) = AN5;
        m_R2N(first + 5) = AN6;
        m_R2N(first + 6) = AN7;
        m_R2N(first + 7) = AN8;
        return i;
}
/*----------------------------------------------------------------------------*/
TCellID
RContainer::newPyramid(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4, const TCellID AN5)
{
        TSize i = addPyramid();
        TSize first = m_first(i);
        m_R2N(first) = AN1;
        m_R2N(first + 1) = AN2;
        m_R2N(first + 2) = AN3;
        m_R2N(first + 3) = AN4;
        m_R2N(first + 4) = AN5;
        return i;
}
/*----------------------------------------------------------------------------*/
TCellID
RContainer::newPrism3(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                          const TCellID AN5, const TCellID AN6)
{
        TSize i = addPrism3();
        TSize first = m_first(i);
        m_R2N(first) = AN1;
        m_R2N(first + 1) = AN2;
        m_R2N(first + 2) = AN3;
        m_R2N(first + 3) = AN4;
        m_R2N(first + 4) = AN5;
        m_R2N(first + 5) = AN6;
        return i;
}
/*----------------------------------------------------------------------------*/
void
RContainer::remove(const TCellID AId)
{
        m_first(AId) = NullTSize;
}
/*----------------------------------------------------------------------------*/
void
RContainer::get(const TCellID AI, Kokkos::View<TCellID*>& AN) const
{
        TSize first = m_first(AI);
        int nb = m_nb_nodes(AI);
        AN = Kokkos::subview(m_R2N, std::make_pair(first, first + nb));
}
/*----------------------------------------------------------------------------*/
void
RContainer::get(const TCellID AI, TCellID AIds[12], int* ASize) const
{
        TSize first = m_first(AI);
        int nb = m_nb_nodes(AI);

        *ASize = nb;
        for(int i=0; i<nb; i++) {
            AIds[i] = m_R2N(first+i);
        }
}
/*----------------------------------------------------------------------------*/
TInt32
RContainer::getNbNodes(const TCellID AI) const
{
        return m_nb_nodes(AI);
}
/*----------------------------------------------------------------------------*/
void
RContainer::set(const TCellID AI, const Kokkos::View<TCellID*>& AN)
{
        TSize first = m_first(AI);
        int nb = m_nb_nodes(AI);
        if (AN.size() <= nb) {
                setUnsafe(AI, AN);
        } else {
                // We don't have enough space to put the new node ids
                // face id doesn't change but first and nb must point on new
                // memory locations in R2N
                // it means so that we let some new holes in the m_R2N RContainer
                nb = AN.size();
                // we get a new first item in R2N
                first = Kokkos::atomic_fetch_add(&m_top_R2N, nb);
                // we update data of face AI now
                m_first(AI) = first;
                m_nb_nodes(AI) = nb;
                for (auto i = 0; i < nb; i++) {
                        m_R2N(first + i) = AN(i);
                }
        }
}
/*----------------------------------------------------------------------------*/
void
RContainer::setUnsafe(const TCellID AI, const Kokkos::View<TCellID*>& AN)
{
        TSize first = m_first(AI);
        int nb = AN.size();
        // this allocation is useful when we replace a hex by a tet
        // for instance
        m_nb_nodes(AI) = nb;
        for (auto i = 0; i < nb; i++) {
                m_R2N(first + i) = AN(i);
        }
}
/*----------------------------------------------------------------------------*/
void
RContainer::setUnsafe(const TCellID AI, TCellID* AIds, int ASize)
{
    TSize first = m_first(AI);
    int nb = ASize;
    // this allocation is useful when we replace a hex by a tet
    // for instance
    m_nb_nodes(AI) = nb;
    for (auto i = 0; i < nb; i++) {
        m_R2N(first + i) = AIds[i];
    }
}
/*----------------------------------------------------------------------------*/
TSize
RContainer::capacity() const
{
        return m_first.extent(0);
}
/*----------------------------------------------------------------------------*/
TSize
RContainer::top() const
{
        return m_top;
}
/*----------------------------------------------------------------------------*/
TSize
RContainer::nbCells() const
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
RContainer::nbTetrahedra() const
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
RContainer::nbHexahedra() const
{
        TSize n = m_top;
        TSize s = 0;
        Kokkos::parallel_reduce(n,
                                KOKKOS_LAMBDA(const TSize i, TSize& ls) {
                                        if (this->has(i) && m_nb_nodes(i) == 8) {
                                                ls = ls + 1;
                                        }
                                },
                                s);
        return s;
}
/*----------------------------------------------------------------------------*/
TSize
RContainer::nbPyramids() const
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
TSize
RContainer::nbPrism3s() const
{
        TSize n = m_top;
        TSize s = 0;
        Kokkos::parallel_reduce(n,
                                KOKKOS_LAMBDA(const TSize i, TSize& ls) {
                                    if (this->has(i) && m_nb_nodes(i) == 6) {
                                            ls = ls + 1;
                                    }
                                },
                                s);
        return s;
}
/*----------------------------------------------------------------------------*/
void
RContainer::resize(const TInt32 ASize)
{
        TSize prev_capacity = capacity();

        if(prev_capacity < ASize) {

                Kokkos::resize(m_first, ASize);
                Kokkos::resize(m_nb_nodes, ASize);
                Kokkos::resize(m_R2N, (TSize) ASize * 8);

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
RContainer::removeAll()
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
