/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The KMDS library is a computer program whose purpose is to provide a set of
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
 * RContainer.h
 *
 *  Created on: 06/22/2017
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_RCONTAINER_H_
#define KMDS_RCONTAINER_H_
/*----------------------------------------------------------------------------*/
// STL headers
#include <iostream>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// KMDS headers
#include <KM/Utils/KTypes.h>
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
class Mesh;
/*----------------------------------------------------------------------------*/
class EXPORT_KMDS RContainer
{
 public:
        /*------------------------------------------------------------------------*/
        /** \brief Constructor
                *
                * \param[in] AOwner the mesh that owns this cell container
                * \param[in] ACapacity initial capacity of the container
                                */
        RContainer(Mesh* AOwner, const TInt32 ACapacity = 16);

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        virtual ~RContainer();

        /*------------------------------------------------------------------------*/
        /** \brief Indicate if this container contains a cell of id AID
         *
         *  \param AID a cell id
         */
        bool has(const TCellID AID) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate \p ANb new cells
         *
         * \param[in] ANb the number of cells to be created
         * \return    The index of the first created cell. Others are contiguous righ behind.
         */
        TCellID addHexahedra(const TInt32 ANb);
        TCellID addTetrahedra(const TInt32 ANb);
        TCellID addPyramids(const TInt32 ANb);
        TCellID addPrism3s(const TInt32 ANb);
        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate one
         *
         * \return    The index of the created cell
         */
        TCellID addHexahedron();
        TCellID addTetrahedron();
        TCellID addPyramid();
        TCellID addPrism3();

        TCellID newHexahedron(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                    const TCellID AN5, const TCellID AN6, const TCellID AN7, const TCellID AN8);
        TCellID newTetrahedron(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4);
        TCellID newPyramid(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                            const TCellID AN5);
        TCellID newPrism3(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                       const TCellID AN5, const TCellID AN6);

        /*------------------------------------------------------------------------*/
        /** \brief  Remove/Desallocate cell of id \p AId.
                 *
         * \param[in] AId The index of the cell to be removed.
         */
        void remove(const TCellID AId);

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the capacity in terms of array size (including holes)
         */
        TSize capacity() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the top id (including holes)
         */
        TSize top() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of cells really stored
         */
        TSize nbCells() const;
        TSize nbTetrahedra() const;
        TSize nbHexahedra() const;
        TSize nbPyramids() const;
        TSize nbPrism3s() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the nodes of cell \p AI
         *
         *  \param[in]  AI cell id we want to get the nodes
         *  \param[out] AN obtained nodes
         */
        void get(const TCellID AI, std::vector<TCellID>& AN) const;
        void get(const TCellID AI, Kokkos::View<TCellID*>& AV) const;
        void get(const TCellID AI, TCellID AIds[12], int* ASize) const;
        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of nodes of cell \p AI
         *
         *  \param[in]  AI cell id we want to get the nodes
         */
        TInt32 getNbNodes(const TCellID AI) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Set the nodes of cell \p AI
         *
         *  \param[in] AI cell id we want to get the nodes
         *  \param[in] AN new nodes of the cell
         */
        void set(const TCellID AI, const Kokkos::View<TCellID*>& AN);

        /*------------------------------------------------------------------------*/
        /** \brief  Set the nodes of cell \p AI. Unlike set(...), this unsafe
         *          version, just replace values but don't check if you have the
         *          right number of nodes.
         *
         *  \param[in] AI cell id we want to get the nodes
         *  \param[in] AN new nodes of the cell
         */
        void setUnsafe(const TCellID AI, const Kokkos::View<TCellID*>& AN);
        void setUnsafe(const TCellID AI, TCellID* AIds, int ASize);

        /*------------------------------------------------------------------------*/
        /** \brief  Resize the container
         */
        void resize(const TInt32 ASize);

        /*------------------------------------------------------------------------*/
        /** \brief  Empty the container
         */
        void removeAll();

 private:
        /** mesh owner of this cell container */
        Mesh* m_mesh;

        /** First node, if = NullID means empty*/
        Kokkos::View<TSize*> m_first;
        /** Nb nodes*/
        Kokkos::View<int*> m_nb_nodes;

        Kokkos::View<TCellID*> m_R2N;

        /** top of the heap including holes so */
        TSize m_top;
        TSize m_top_R2N;
};
/*----------------------------------------------------------------------------*/
}  // namespace KMDS
/*----------------------------------------------------------------------------*/
#endif /* KMDS_RCONTAINER_H_ */
/*----------------------------------------------------------------------------*/
