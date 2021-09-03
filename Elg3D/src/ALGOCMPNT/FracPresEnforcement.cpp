/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    Tools.cpp
 *  \author  legoff
 *  \date    01/14/2020
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/FracPresEnforcement.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <map>
#include <set>
#include <string>
#include <vector>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
    double FracPresEnforcement_maintain_xD(const kmds::GrowingView<kmds::TCellID>* cellIDs,
                                           const elg3d::FracPres* Afp_ref,
                                           elg3d::FracPres* Afp_new,
                                           const std::vector<std::string> AMatNames)
    {
        Afp_new->clear();

        const int nbMat = Afp_ref->getNbMaterials();

        // copy the frac pres data
        *Afp_new = *Afp_ref;

        std::vector<int> matId2remove_vec;
        for(auto name: AMatNames) {
            matId2remove_vec.push_back(Afp_ref->getMaterialID(name));
        }
        const int nbMat2remove = matId2remove_vec.size();

        const kmds::TCellID nbCells = cellIDs->getNbElems();

        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID cid = cellIDs->get(i);

                                 for(int imat=0; imat<nbMat2remove; imat++) {
                                     double fp = Afp_ref->getFracPres(matId2remove_vec[imat], cid);

                                     if(fp > 0.) {
                                         if(Afp_ref->getMaxMatFracPresIndex(cid) == matId2remove_vec[imat]) {
                                             break;
                                         } else {
                                             double fp_max = Afp_ref->getMaxMatFracPresValue(cid);
                                             double fp_modif = fp_max + 0.01;
                                             Afp_new->setFracPres(matId2remove_vec[imat], cid, fp_modif);

                                             Afp_new->normalize(cid);
                                             break;
                                         }
                                     }

                                 }
                             });

    }
    /*----------------------------------------------------------------------------*/
    double FracPresEnforcement_fuse_xD(const kmds::GrowingView<kmds::TCellID>* cellIDs,
                                       const elg3d::FracPres* Afp_ref,
                                       elg3d::FracPres* Afp_new,
                                       const std::vector<std::string> AMatNames,
                                       const std::string AMatName,
                                       const bool AKeepNumbering)
    {
        Afp_new->clear();

        const int nbMat_ref = Afp_ref->getNbMaterials();
        const std::map<int, std::string> matList_ref = Afp_ref->getMaterialList();

        std::set<int> matId2remove;
        std::vector<int> matId2remove_vec;
        for(auto name: AMatNames) {
            matId2remove.insert(Afp_ref->getMaterialID(name));
            matId2remove_vec.push_back(Afp_ref->getMaterialID(name));
        }
        const int nbMat2remove = matId2remove.size();

        // build the new material list
        std::map<int, std::string> matList_new;

        if(AKeepNumbering) {
            matList_new = matList_ref;
        } else {
            for(auto mat: matList_ref) {
                if(matId2remove.end() == matId2remove.find(mat.first)) {
                    matList_new.emplace(mat);
                }
            }
        }

        Afp_new->setMaterialList(matList_new);

        // set the unmodified material fracpres
        for(auto mat: matList_ref) {
            if(matId2remove.end() == matId2remove.find(mat.first)) {
                Afp_new->setMatFracPres(Afp_new->getMaterialID(mat.second), Afp_ref->getMatFracPres(mat.first));
            }
        }

        // create the new material and compute its fracpres
        const int matID = Afp_new->createMaterial(AMatName);
        const kmds::TCellID nbCells = cellIDs->getNbElems();

        Kokkos::View<double*> fracpres_newMat("FRACPRES_NEWMAT", nbCells);

        Kokkos::parallel_for(nbCells,
                             KOKKOS_LAMBDA(const kmds::TCellID i) {
                                 const kmds::TCellID cid = cellIDs->get(i);
                                 double fp_sum = 0.;

                                 for(int imat=0; imat<nbMat2remove; imat++) {
                                     fp_sum += Afp_ref->getFracPres(matId2remove_vec[imat], cid);
                                 }
                                 fracpres_newMat[i] = fp_sum;
                             });

        // TODO: currently cannot be done in parallel because of the FracPres non-threadsafe structure
        for(kmds::TCellID i=0; i<nbCells; i++) {
            Afp_new->setFracPres(matID, cellIDs->get(i), fracpres_newMat[i]);
        }

    }
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
