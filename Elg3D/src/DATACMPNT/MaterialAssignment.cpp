/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    MaterialAssignment.cpp
 *  \author  legoff
 *  \date    04/12/2017
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/
    MaterialAssignment::MaterialAssignment(kmds::TCellID ACapacity)
    : m_cells2Mat("CELLS2MAT", ACapacity)
    {

    }
/*----------------------------------------------------------------------------*/
//MaterialAssignment::MaterialAssignment(const MaterialAssignment& AMaterialAssignment)
//{
//}
///*----------------------------------------------------------------------------*/
//MaterialAssignment::~MaterialAssignment()
//{
//}
///*----------------------------------------------------------------------------*/
//MaterialAssignment&
//MaterialAssignment::operator=(const MaterialAssignment& AMaterialAssignment)
//{
//        return *this;
//}
/*----------------------------------------------------------------------------*/
    int
    MaterialAssignment::getNbMaterials() const
    {
        return m_matId2Name.size();
    }
/*----------------------------------------------------------------------------*/
    int
    MaterialAssignment::getMaterialID(std::string AName) const
    {
        return m_matName2Id.at(AName);
    }
/*----------------------------------------------------------------------------*/
    std::string
    MaterialAssignment::getMaterialName(int AId) const
    {
        return m_matId2Name.at(AId);
    }
/*----------------------------------------------------------------------------*/
    std::map<int, std::string>
    MaterialAssignment::getMaterialList() const
    {
        return m_matId2Name;
    }
/*----------------------------------------------------------------------------*/
    int
    MaterialAssignment::createMaterial(std::string AName)
    {
        const int id = m_matId2Name.size();
        m_matName2Id.emplace(AName, id);
        m_matId2Name.emplace(id, AName);

        return id;
    }
/*----------------------------------------------------------------------------*/
    void
    MaterialAssignment::createMaterials(std::map<int, std::string> AMatList)
    {
        for(auto name: AMatList) {
            this->createMaterial(name.second);
        }
    }
/*----------------------------------------------------------------------------*/
    void
    MaterialAssignment::setMaterial(int AMat, kmds::TCellID AId)
    {
        m_cells2Mat[AId] = AMat;
    }
/*----------------------------------------------------------------------------*/
    int
    MaterialAssignment::getMaterial(kmds::TCellID AId) const
    {
        return m_cells2Mat[AId];
    }
/*----------------------------------------------------------------------------*/
//    std::set<kmds::TCellID>
//    MaterialAssignment::getCells(int AMat) const
//    {
//        //return m_mat2Cells.at(AMat);
//    }
    /*----------------------------------------------------------------------------*/
    void
    MaterialAssignment::changeMaterial(int AMat, kmds::TCellID AId)
    {
        m_cells2Mat[AId] = AMat;
    }
    /*----------------------------------------------------------------------------*/
    kmds::TSize
    MaterialAssignment::capacity() const
    {
        return m_cells2Mat.extent(0);
    }
    /*----------------------------------------------------------------------------*/
    void
    MaterialAssignment::updateCapacity(const kmds::TSize ASize)
    {
        int prev_capacity = capacity();

        if(prev_capacity < ASize) {

            Kokkos::resize(m_cells2Mat, ASize);
        }
    }
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
