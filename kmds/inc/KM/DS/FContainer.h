/*----------------------------------------------------------------------------*/
/*
 * FContainer.h
 *
 *  Created on: 03/08/2017
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_FCONTAINER_H_
#define KMDS_FCONTAINER_H_
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
class EXPORT_KMDS FContainer
{
 public:
        /*------------------------------------------------------------------------*/
        /** \brief Constructor
                *
                * \param[in] AOwner the mesh that owns this cell container
                * \param[in] ACapacity initial capacity of the container
                                */
        FContainer(Mesh* AOwner, const TInt32 ACapacity = 16);

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        virtual ~FContainer();

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
        TCellID addTriangles(const TInt32 ANb);
        TCellID addQuads(const TInt32 ANb);
        TCellID addPentagons(const TInt32 ANb);
        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate one
         *
         * \return    The index of the created cell
         */
        TCellID addQuad();
        TCellID addTriangle();
        TCellID addPentagon();

        TCellID newQuad(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4);
        TCellID newTriangle(const TCellID AN1, const TCellID AN2, const TCellID AN3);
        TCellID newPentagon(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4,
                            const TCellID AN5);

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
        TSize nbTriangles() const;
        TSize nbQuads() const;
        TSize nbPentagons() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the nodes of cell \p AI
         *
         *  \param[in]  AI cell id we want to get the nodes
         *  \param[out] AN obtained nodes
         */
        void get(const TCellID AI, std::vector<TCellID>& AN) const;
        void get(const TCellID AI, Kokkos::View<TCellID*>& AV) const;
        void get(const TCellID AI, TCellID AIds[5], int* ASize) const;
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

        Kokkos::View<TCellID*> m_F2N;
        /** top of the heap including holes so */
        TSize m_top;
        TSize m_top_F2N;
};
/*----------------------------------------------------------------------------*/
}  // namespace KMDS
/*----------------------------------------------------------------------------*/
#endif /* KMDS_FCONTAINER_H_ */
/*----------------------------------------------------------------------------*/
