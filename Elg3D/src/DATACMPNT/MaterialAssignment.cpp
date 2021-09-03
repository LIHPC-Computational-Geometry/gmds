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
