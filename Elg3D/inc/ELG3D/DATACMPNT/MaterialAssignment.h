/*----------------------------------------------------------------------------*/

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
