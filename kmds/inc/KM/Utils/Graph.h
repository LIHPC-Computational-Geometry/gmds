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
/** \file    Graph.h
 *  \author  legoff
 *  \date    03/06/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_GRAPH_H_
#define KMDS_GRAPH_H_
/*----------------------------------------------------------------------------*/
// Kokkos File Headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
#include <set>
#include <vector>
/*----------------------------------------------------------------------------*/
// KMDS File Headers
#include "KM/Utils/GrowingView.h"
#include "KM/Utils/KTypes.h"
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
/** \class Graph
 *  \brief
 */
/*----------------------------------------------------------------------------*/
class Graph
{
 public:
        /*------------------------------------------------------------------------*/
        /** \brief Default constructor.
         *
         *  \param[in] AName view name for Kokkos debug purpose
         *  \param[in] ANbVec number of vertices in the graph
         *  \param[in] ANbVec max number of neighbors per vertices
         */
        Graph(const std::string AName, const TSize ANbVec, const TSize AMaxNeighbors)
         : m_nbvec(ANbVec)
         , m_maxneighbors(AMaxNeighbors)
         , m_nbneighbors(AName, ANbVec)
         , m_neighbors(AName, ANbVec, AMaxNeighbors)
        {
                ;
        }

        /*------------------------------------------------------------------------*/
        /** \brief Getter on the number of vertices
         *
         *  \return the number of vertices
         */
        TCellID getNbVec() const;

    /*------------------------------------------------------------------------*/
    /** \brief Getter on the number of edges
     *
     *  \return the number of edges
     */
    TCellID getNbEdges() const;

        /*------------------------------------------------------------------------*/
        /** \brief Set the neighbors of a vertex
         *
         *  \param[in]  AId the vertex id
         *  \param[in]  AN the indices of the neighbors
         *
         */
        void setNeighbors(const kmds::TCellID AId, const std::set<kmds::TCellID>& AN);

        /*------------------------------------------------------------------------*/
        /** \brief Set the neighbors of a vertex
         *
         *  \param[in]  AId the vertex id
         *  \param[in]  AN the indices of the neighbors
         *
         */
        void setNeighbors(const kmds::TCellID AId, const std::vector<kmds::TCellID>& AN);

        /*------------------------------------------------------------------------*/
        /** \brief Get the neighbors of a vertex
         *
         *  \param[in]  AId the vertex id
         *  \param[out] AN a subview containing the neighbors
         *
         *  \return the number of neighbors
         */
        TCellID getNeighbors(const TCellID AId, Kokkos::View<TCellID*>& AN) const;

        /*------------------------------------------------------------------------*/
        /** \brief Check whether AId1 is a neighbor of AId0 (directed edge AID0 -> AId1)
        *
        *  \param[in]  AId0 the vertex of id
        *  \param[in]  AId1 the vertex of id
        *
        *  \return whether AId1 is a neighbor of AId0
        */
        bool isNeighbor(const TCellID AId0, const TCellID AId1) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns a container with ids of graph vertices forming an independent set.
         *          Will not work if the graph contains self-neighboring nodes
         *
         *          // WARNING works only on non-oriented graphs
         *
         *  \param[out] ASelection a container with independent vertices
         */
        void getIndependentSet(kmds::GrowingView<kmds::TCellID>* ASelection) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns whether the selection of vertices forms an independent set
         *
         *          // WARNING works only on non-oriented graphs
         *
         *  \param[in] ASelection a container of vertices to check
         *
         *  \return true if the set is independent, else otherwise
         */
        bool checkIndependentSet(kmds::GrowingView<kmds::TCellID>* ASelection) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns a container with ids of graph vertices forming an independent set.
         *          Will not work if the graph contains self-neighboring nodes
         *
         *          // WARNING works only on non-oriented graphs
         *
         *  \param[out] ASelection a container with independent vertices
         */
        void buildColoring(const int nbMaxColors);


        kmds::GrowingView<kmds::TCellID>* getColoring(const int AColor);

        /*------------------------------------------------------------------------*/
        /** \brief  Populates the graph
         *
         *  \param[in]  AMeanNbNeighbors the average number of neighbors per vertex
         */
        void populateNonOriented(TFloat AMeanNbNeighbors) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Prints the graph
         *
         */
        void print() const;

        // WARNING hard-coded value of 20 for MaxNbColors
        const static int MaxNbColors = 20;

        TCellID m_nbvec;
        TCellID m_maxneighbors;
        Kokkos::View<TCellID*> m_nbneighbors;
        Kokkos::View<TCellID**> m_neighbors;

        std::vector<kmds::GrowingView<kmds::TCellID> > m_colorsets;
};
/*----------------------------------------------------------------------------*/
}  // end namespace kmds
/*----------------------------------------------------------------------------*/
#endif /* KMDS_GRAPH_H_ */
/*----------------------------------------------------------------------------*/
