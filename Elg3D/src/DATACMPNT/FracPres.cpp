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
/** \file    FracPres.cpp
 *  \author  legoff
 *  \date    20/11/2017
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/FracPres.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <set>
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
    FracPres::FracPres()
    {
    }
/*----------------------------------------------------------------------------*/
//FracPres::FracPres(const FracPres& AFracPres)
//{
//}
///*----------------------------------------------------------------------------*/
//FracPres::~FracPres()
//{
//}
///*----------------------------------------------------------------------------*/
FracPres&
FracPres::operator=(const FracPres& AFracPres)
{
    this->m_matId2Name = AFracPres.m_matId2Name;
    this->m_matName2Id = AFracPres.m_matName2Id;
    this->m_mat2CellFracPres = AFracPres.m_mat2CellFracPres;
    return *this;
}
/*----------------------------------------------------------------------------*/
    void
    FracPres::clear()
    {
        m_matId2Name.clear();
        m_matName2Id.clear();
        m_mat2CellFracPres.clear();
    }
/*----------------------------------------------------------------------------*/
    int
    FracPres::getNbMaterials() const
    {
        return m_matId2Name.size();
    }
/*----------------------------------------------------------------------------*/
    int
    FracPres::getMaterialID(std::string AName) const
    {
        return m_matName2Id.at(AName);
    }
/*----------------------------------------------------------------------------*/
    std::string
    FracPres::getMaterialName(int AId) const
    {
        return m_matId2Name.at(AId);
    }
/*----------------------------------------------------------------------------*/
    std::map<int, std::string>
    FracPres::getMaterialList() const
    {
        return m_matId2Name;
    }
    /*----------------------------------------------------------------------------*/
    void
    FracPres::setMaterialList(const std::map<int, std::string> AMatList)
    {
        for(auto mat: AMatList) {
            m_matName2Id.emplace(mat.second, mat.first);
            m_matId2Name.emplace(mat.first, mat.second);
            m_mat2CellFracPres.emplace(mat.first, std::map<kmds::TCellID, double> () );
        }
    }
    /*----------------------------------------------------------------------------*/
    bool
    FracPres::materialExists(std::string AName) const
    {
        return (m_matName2Id.find(AName) != m_matName2Id.end());
    }
/*----------------------------------------------------------------------------*/
    int
    FracPres::createMaterial(std::string AName)
    {
        const int id = m_matId2Name.size();
        m_matName2Id.emplace(AName, id);
        m_matId2Name.emplace(id, AName);
        m_mat2CellFracPres.emplace(id, std::map<kmds::TCellID, double> () );

        return id;
    }
/*----------------------------------------------------------------------------*/
    void
    FracPres::setFracPres(int AMat, kmds::TCellID AId, double AFracPres)
    {
        m_mat2CellFracPres[AMat][AId] = AFracPres;
    }
/*----------------------------------------------------------------------------*/
    double
    FracPres::getFracPres(int AMat, kmds::TCellID AId) const
    {
        if(m_mat2CellFracPres.at(AMat).find(AId) == m_mat2CellFracPres.at(AMat).end()) {
            return 0.;
        } else {
            return m_mat2CellFracPres.at(AMat).at(AId);
        }
    }
/*----------------------------------------------------------------------------*/
    bool
    FracPres::isMixedCell(kmds::TCellID AId) const
    {
        for(int imat=0; imat<this->getNbMaterials(); imat++) {
            double matfp = this->getFracPres(imat, AId);

            if((matfp != 0.) && (matfp != 1.)) {
//            if((matfp > FracPres_TRACE_AMOUNT) && (matfp < FracPres_TRACE_AMOUNT)) {
                return true;
            }
        }

        return false;
    }
/*----------------------------------------------------------------------------*/
    int
    FracPres::getMaxMatFracPresIndex(kmds::TCellID AId) const
    {
        int matIndex = -1;
        double maxfracpres = -HUGE_VALF;

        for(int imat=0; imat<this->getNbMaterials(); imat++) {
            double matfp = this->getFracPres(imat, AId);
            if((matfp != 0.) && (matfp > maxfracpres)) {
                maxfracpres = matfp;
                matIndex = imat;
            }
        }

        return matIndex;
    }
/*----------------------------------------------------------------------------*/
    double
    FracPres::getMaxMatFracPresValue(kmds::TCellID AId) const
    {
        double maxfracpres = -HUGE_VALF;

        for(int imat=0; imat<this->getNbMaterials(); imat++) {
            double matfp = this->getFracPres(imat, AId);
            if((matfp != 0.) && (matfp > maxfracpres)) {
                maxfracpres = matfp;
            }
        }

        return maxfracpres;
    }
/*----------------------------------------------------------------------------*/
    std::map<kmds::TCellID, double>
    FracPres::getMatFracPres(int AMat) const
    {
        return m_mat2CellFracPres.at(AMat);
    }
/*----------------------------------------------------------------------------*/
    void
    FracPres::setMatFracPres(int AMat, const std::map<kmds::TCellID, double> AMatFracPres)
    {
        m_mat2CellFracPres.at(AMat) = AMatFracPres;
    }
/*----------------------------------------------------------------------------*/
    void
    FracPres::normalize(kmds::TCellID AId)
    {
        double fp_tot = 0.;
        const int nbMat = this->getNbMaterials();

        for(int imat=0; imat<nbMat; imat++) {
            double matfp = this->getFracPres(imat, AId);
            fp_tot += matfp;
        }

        for(int imat=0; imat<nbMat; imat++) {
            double matfp = this->getFracPres(imat, AId);
            if(matfp > 0.) {
                matfp *= 1. / fp_tot;
                this->setFracPres(imat, AId, matfp);
            }
        }

    }
/*----------------------------------------------------------------------------*/
bool
FracPres::normalize()
    {
        bool wasModified = false;

        // first gather all the cells
        std::set<kmds::TCellID> cellIDs;

        const std::map<int, std::string> mats = this->getMaterialList();

        for(auto m: mats) {
            const int matid = m.first;

            const std::map<kmds::TCellID, double> cellsfp = this->getMatFracPres(matid);

            for(auto cfp: cellsfp) {
                cellIDs.insert(cfp.first);
            }
        }

        for(auto cid: cellIDs) {
            double fp_tot = 0.;

            for(auto m: mats) {
                double matfp = this->getFracPres(m.first, cid);
                fp_tot += matfp;
            }

            if(std::abs(1. - fp_tot) > 10e-12) {
                wasModified = true;
                this->normalize(cid);
            }
        }

        return wasModified;
    }
/*----------------------------------------------------------------------------*/
bool
FracPres::checkNormalized() const
    {
        // first gather all the cells
        std::set<kmds::TCellID> cellIDs;

        const std::map<int, std::string> mats = this->getMaterialList();

        for(auto m: mats) {
            const int matid = m.first;

            const std::map<kmds::TCellID, double> cellsfp = this->getMatFracPres(matid);

            for(auto cfp: cellsfp) {
                cellIDs.insert(cfp.first);
            }
        }

        for(auto cid: cellIDs) {
            double fp_tot = 0.;

            for(auto m: mats) {
                double matfp = this->getFracPres(m.first, cid);
                fp_tot += matfp;
            }

            if(std::abs(1. - fp_tot) > 10e-12) {
                std::cerr<<"FracPres::checkNormalized false at cellid "<<cid<<" with fp_tot "<<fp_tot<<std::endl;
                return false;
            }
        }

        return true;
    }
/*----------------------------------------------------------------------------*/
bool
FracPres::removeTraceAmounts()
{
    bool hasChanged = false;

    // first gather all the cells
    std::set<kmds::TCellID> cellIDs;

    const std::map<int, std::string> mats = this->getMaterialList();

    for(auto m: mats) {
        const int matid = m.first;

        const std::map<kmds::TCellID, double> cellsfp = this->getMatFracPres(matid);

        for(auto cfp: cellsfp) {
            cellIDs.insert(cfp.first);
        }
    }

    for(auto cid: cellIDs) {
        double fp_tot = 0.;
        int nbMat_loc = 0;
        int mat_index = -1;

        for(auto m: mats) {
            double matfp = this->getFracPres(m.first, cid);

            // remove trace materials from the cell ( fp < 10e-10)
            if(matfp <= FracPres_TRACE_AMOUNT) {
                this->m_mat2CellFracPres[m.first].erase(cid);
                hasChanged = true;
            } else {
                fp_tot += matfp;
                nbMat_loc++;
                mat_index = m.first;
            }
        }

        // force fp to 1. when only one material present
        if(nbMat_loc == 1) {
            const double fp_prev = this->getFracPres(mat_index, cid);
            this->setFracPres(mat_index, cid, 1.);
            if(fp_prev != 1.) {
                hasChanged = true;
            }
        } else {
            if (std::abs(1. - fp_tot) >= FracPres_TRACE_AMOUNT) {
                this->normalize(cid);
                hasChanged = true;
            }
        }
    }

    return hasChanged;
}
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
