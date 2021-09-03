/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * This software is a computer program whose purpose is to provide a set of
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
 * data to be ensured and, more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    MaterialAssignment.h
 *  \author  legoff
 *  \date    30/11/2017
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_MATERIALASSIGNMENT_H_
#define ELG3D_MATERIALASSIGNMENT_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <map>
#include <set>
#include <string>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/Utils/KTypes.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/
/** \class MaterialAssignment
 *  \brief This is a dummy class.
 */
/*----------------------------------------------------------------------------*/
    class MaterialAssignment
    {
    public:
        /*------------------------------------------------------------------------*/
        /** \brief  Constructor
         */
        MaterialAssignment(kmds::TCellID ACapacity = 16);

        /*------------------------------------------------------------------------*/
        /** \brief Copy constructor
         */
        MaterialAssignment(const MaterialAssignment&) =delete;

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        ~MaterialAssignment() =default;

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator=
         */
        MaterialAssignment& operator=(const MaterialAssignment&) =delete;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of materials
         *
         *  \return the number of materials
         */
        int getNbMaterials() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the id of a material
         *
         *  \param[in]  AName the name of the material
         *
         *  \return the id of the material
         */
        int getMaterialID(std::string AName) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the name of a material
         *
         *  \param[in]  AId the id of the material
         *
         *  \return the name of the material
         */
        std::string getMaterialName(int AId) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the names of the materials
         *
         *  \return the name of the material
         */
        std::map<int, std::string> getMaterialList() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Adds a material
         *
         *  \param[in]  AName the name of the material
         *
         *  \return the id of the added material
         */
        int createMaterial(std::string AName);

        /*------------------------------------------------------------------------*/
        /** \brief  Adds a material
         *
         *  \param[in]  AMatList the list of the materials names and their ids
         */
// TODO create with specified ID
        void createMaterials(std::map<int, std::string> AMatList);

        /*------------------------------------------------------------------------*/
        /** \brief  Assigns a cell to a material
         *
         *  \param[in]  AMat the id of the material
         *  \param[in]  AId the cell id
         */
// TODO threadsafety
        void setMaterial(int AMat, kmds::TCellID AId);

        /*------------------------------------------------------------------------*/
        /** \brief  Gets the MaterialAssignment of a cell
         *
         *  \param[in]  AId the cell id
         *
         *  \return  the MaterialAssignment of cell AId
         */
        int getMaterial(kmds::TCellID AId) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Gets the cells assigned to a material
         *
         *  \param[in]  AMat the id of the material
         *
         *  \return  a container of the IDs of the cells assigned to AMat
         */
        std::set<kmds::TCellID> getCells(int AMat) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Changes the assignment of a cell
         *
         *  \param[in]  AMat the id of the material
         *  \param[in]  AId the cell id
         */
// TODO threadsafety
        void changeMaterial(int AMat, kmds::TCellID AId);

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the capacity in terms of array size (including holes)
         */
        kmds::TSize capacity() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Resize the container
         */
        void updateCapacity(const kmds::TSize ASize);

    private:

        /* mapping between material ids and their names */
        std::map<int, std::string> m_matId2Name;

        /* mapping between material names and their ids */
        std::map<std::string, int> m_matName2Id;

        /* mapping between cells and their material assignment, for each material */
        //std::map<int, std::set<kmds::TCellID> > m_mat2Cells;

        /* mapping between cells and their material assignment, for each material */
        //std::map<kmds::TCellID, int> m_cells2Mat;
        Kokkos::View<int *> m_cells2Mat;


    };
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_MATERIALASSIGNMENT_H_ */
