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
/** \file    InitData.cpp
 *  \author  legoff
 *  \date    22/11/2017
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/InitData.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <sstream>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
// exodusII
#ifdef ELG3D_WITH_EXODUSII
#include <exodusII.h>
#endif
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
// KMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/Utils/InitTools.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/Tools.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/
    void
    initData_3x3_2D(kmds::Mesh* AMesh, elg3d::FracPres* Afp)
    {
        const int nb_x = 4;
        const int nb_y = 4;
        const int nb_z = 1;
        const int nb_i = nb_x - 1;
        const int nb_j = nb_y - 1;

        AMesh->updateNodeCapacity(nb_x*nb_y);
        AMesh->updateFaceCapacity(nb_i*nb_j);

        int indexNode = AMesh->addNodes(nb_x*nb_y*nb_z);
        for(int i=0; i<nb_x; i++) {
            for(int j=0; j<nb_y; j++) {
                    AMesh->setNodeLocation(indexNode,i,j,0.);
                    indexNode++;
            }
        }

        kmds::TCellID indexCell = AMesh->addQuads(nb_i*nb_j);
        for(int i=0; i<nb_i; i++) {
            for(int j=0; j<nb_j; j++) {
                kmds::Face q = AMesh->getFace(indexCell);
                kmds::TCellID v[4];
                v[0] = i*nb_y + j;
                v[1] = (i+1)*nb_y + j;
                v[2] = (i+1)*nb_y + j+1;
                v[3] = i*nb_y + j+1;
                q.setNodes(v,4);
                indexCell++;
            }
        }

        int aId = Afp->createMaterial("matA");
        int bId = Afp->createMaterial("matB");
        int cId = Afp->createMaterial("matC");

        Afp->setFracPres(aId, 0, 1.);
        Afp->setFracPres(aId, 1, 0.8);
        Afp->setFracPres(aId, 3, 0.9);
        Afp->setFracPres(aId, 4, 0.72);
        Afp->setFracPres(bId, 1, 0.2);
        Afp->setFracPres(bId, 2, 1.);
        Afp->setFracPres(bId, 4, 0.18);
        Afp->setFracPres(bId, 5, 0.9);
        Afp->setFracPres(cId, 3, 0.1);
        Afp->setFracPres(cId, 4, 0.1);
        Afp->setFracPres(cId, 5, 0.1);
        Afp->setFracPres(cId, 6, 1.);
        Afp->setFracPres(cId, 7, 1.);
        Afp->setFracPres(cId, 8, 1.);

    }
/*----------------------------------------------------------------------------*/
    void
    initData_2x2_2D_nonmanifold(kmds::Mesh* AMesh, elg3d::FracPres* Afp)
    {
        const int nb_x = 3;
        const int nb_y = 3;
        const int nb_z = 1;
        const int nb_i = nb_x - 1;
        const int nb_j = nb_y - 1;

        AMesh->updateNodeCapacity(nb_x*nb_y);
        AMesh->updateFaceCapacity(nb_i*nb_j);

        int indexNode = AMesh->addNodes(nb_x*nb_y*nb_z);
        for(int i=0; i<nb_x; i++) {
            for(int j=0; j<nb_y; j++) {
                AMesh->setNodeLocation(indexNode,i,j,0.);
                indexNode++;
            }
        }

        kmds::TCellID indexCell = AMesh->addQuads(nb_i*nb_j);
        for(int i=0; i<nb_i; i++) {
            for(int j=0; j<nb_j; j++) {
                kmds::Face q = AMesh->getFace(indexCell);
                kmds::TCellID v[4];
                v[0] = i*nb_y + j;
                v[1] = (i+1)*nb_y + j;
                v[2] = (i+1)*nb_y + j+1;
                v[3] = i*nb_y + j+1;
                q.setNodes(v,4);
                indexCell++;
            }
        }

        int aId = Afp->createMaterial("matA");
        int bId = Afp->createMaterial("matB");

        Afp->setFracPres(aId, 0, 1.);
        Afp->setFracPres(aId, 1, 0.2);
        Afp->setFracPres(aId, 2, 0.3);
        Afp->setFracPres(aId, 3, 0.9);
        Afp->setFracPres(bId, 1, 0.8);
        Afp->setFracPres(bId, 2, 0.7);
        Afp->setFracPres(bId, 3, 0.1);
    }
/*----------------------------------------------------------------------------*/
    void
    initData_3x3x3_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp)
    {
        const int nb_x = 4;
        const int nb_y = 4;
        const int nb_z = 4;
        const int nb_i = nb_x - 1;
        const int nb_j = nb_y - 1;
        const int nb_k = nb_z - 1;

        AMesh->updateNodeCapacity(nb_x*nb_y*nb_z);
        AMesh->updateRegionCapacity(nb_i*nb_j*nb_k);

        int indexNode = AMesh->addNodes(nb_x*nb_y*nb_z);
        for(int i=0; i<nb_x; i++) {
            for(int j=0; j<nb_y; j++) {
                for(int k=0; k<nb_z; k++) {
                    AMesh->setNodeLocation(indexNode, i, j, k);
                    indexNode++;
                }
            }
        }



        kmds::TCellID indexCell = AMesh->addHexahedra(nb_i*nb_j*nb_k);
        for(int i=0; i<nb_i; i++) {
            for(int j=0; j<nb_j; j++) {
                for(int k=0; k<nb_k; k++) {
                    kmds::Region h = AMesh->getRegion(indexCell);
                    kmds::TCellID v[8];
                    v[0] = i * nb_y * nb_z + j * nb_z + k;
                    v[1] = (i + 1) * nb_y * nb_z + j * nb_z + k;
                    v[2] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k;
                    v[3] = i * nb_y * nb_z + (j + 1) * nb_z + k;
                    v[4] = i * nb_y * nb_z + j * nb_z + k +1;
                    v[5] = (i + 1) * nb_y * nb_z + j * nb_z + k +1;
                    v[6] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k +1;
                    v[7] = i * nb_y * nb_z + (j + 1) * nb_z + k + 1;
                    h.setNodes(v, 8);
                    indexCell++;
                }
            }
        }

        int aId = Afp->createMaterial("matA");
        int bId = Afp->createMaterial("matB");
        int cId = Afp->createMaterial("matC");

        Afp->setFracPres(aId, 0, 1.);
        Afp->setFracPres(aId, 3, 1.);
        Afp->setFracPres(aId, 6, 1.);
        Afp->setFracPres(aId, 1, 0.8);
        Afp->setFracPres(aId, 4, 0.8);
        Afp->setFracPres(aId, 7, 0.8);
        Afp->setFracPres(aId, 9, 0.9);
        Afp->setFracPres(aId, 12, 0.9);
        Afp->setFracPres(aId, 15, 0.9);
        Afp->setFracPres(aId, 10, 0.72);
        Afp->setFracPres(aId, 13, 0.72);
        Afp->setFracPres(aId, 16, 0.72);
        Afp->setFracPres(bId, 1, 0.2);
        Afp->setFracPres(bId, 4, 0.2);
        Afp->setFracPres(bId, 7, 0.2);
        Afp->setFracPres(bId, 2, 1.);
        Afp->setFracPres(bId, 5, 1.);
        Afp->setFracPres(bId, 8, 1.);
        Afp->setFracPres(bId, 10, 0.18);
        Afp->setFracPres(bId, 13, 0.18);
        Afp->setFracPres(bId, 16, 0.18);
        Afp->setFracPres(bId, 11, 0.9);
        Afp->setFracPres(bId, 14, 0.9);
        Afp->setFracPres(bId, 17, 0.9);
        Afp->setFracPres(cId, 9, 0.1);
        Afp->setFracPres(cId, 12, 0.1);
        Afp->setFracPres(cId, 15, 0.1);
        Afp->setFracPres(cId, 10, 0.1);
        Afp->setFracPres(cId, 13, 0.1);
        Afp->setFracPres(cId, 16, 0.1);
        Afp->setFracPres(cId, 11, 0.1);
        Afp->setFracPres(cId, 14, 0.1);
        Afp->setFracPres(cId, 17, 0.1);
        Afp->setFracPres(cId, 18, 1.);
        Afp->setFracPres(cId, 21, 1.);
        Afp->setFracPres(cId, 24, 1.);
        Afp->setFracPres(cId, 19, 1.);
        Afp->setFracPres(cId, 22, 1.);
        Afp->setFracPres(cId, 25, 1.);
        Afp->setFracPres(cId, 20, 1.);
        Afp->setFracPres(cId, 23, 1.);
        Afp->setFracPres(cId, 26, 1.);

    }
    /*----------------------------------------------------------------------------*/
    void
    initData_2x1x2_3D_nonmanifold(kmds::Mesh* AMesh, elg3d::FracPres* Afp)
    {
        const int nb_x = 3;
        const int nb_y = 2;
        const int nb_z = 3;
        const int nb_i = nb_x - 1;
        const int nb_j = nb_y - 1;
        const int nb_k = nb_z - 1;

        AMesh->updateNodeCapacity(nb_x*nb_y*nb_z);
        AMesh->updateRegionCapacity(nb_i*nb_j*nb_k);

        int indexNode = AMesh->addNodes(nb_x*nb_y*nb_z);
        for(int i=0; i<nb_x; i++) {
            for(int j=0; j<nb_y; j++) {
                for(int k=0; k<nb_z; k++) {
                    AMesh->setNodeLocation(indexNode, i, j, k);
                    indexNode++;
                }
            }
        }

        kmds::TCellID indexCell = AMesh->addHexahedra(nb_i*nb_j*nb_k);
        for(int i=0; i<nb_i; i++) {
            for(int j=0; j<nb_j; j++) {
                for(int k=0; k<nb_k; k++) {
                    kmds::Region h = AMesh->getRegion(indexCell);
                    kmds::TCellID v[8];
                    v[0] = i * nb_y * nb_z + j * nb_z + k;
                    v[1] = (i + 1) * nb_y * nb_z + j * nb_z + k;
                    v[2] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k;
                    v[3] = i * nb_y * nb_z + (j + 1) * nb_z + k;
                    v[4] = i * nb_y * nb_z + j * nb_z + k +1;
                    v[5] = (i + 1) * nb_y * nb_z + j * nb_z + k +1;
                    v[6] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k +1;
                    v[7] = i * nb_y * nb_z + (j + 1) * nb_z + k + 1;
                    h.setNodes(v, 8);
                    indexCell++;
                }
            }
        }

        int aId = Afp->createMaterial("matA");
        int bId = Afp->createMaterial("matB");

        Afp->setFracPres(aId, 0, 1.);
        Afp->setFracPres(aId, 1, 0.2);
        Afp->setFracPres(aId, 2, 0.3);
        Afp->setFracPres(aId, 3, 0.9);
        Afp->setFracPres(bId, 1, 0.8);
        Afp->setFracPres(bId, 2, 0.7);
        Afp->setFracPres(bId, 3, 0.1);
    }
/*----------------------------------------------------------------------------*/
    void
    initData_I_2D_nonmanifold(kmds::Mesh* AMesh, elg3d::FracPres* Afp, int nb_i)
    {
        const int nb_j = nb_i;

        const int nb_x = nb_i+1;
        const int nb_y = nb_j+1;

        AMesh->updateNodeCapacity(nb_x*nb_y);
        AMesh->updateFaceCapacity(nb_i*nb_j);

        int indexNode = AMesh->addNodes(nb_x*nb_y);
        for(int i=0; i<nb_x; i++) {
            for(int j=0; j<nb_y; j++) {
                AMesh->setNodeLocation(indexNode, i, j, 0);
                indexNode++;
            }
        }

        kmds::TCellID indexCell = AMesh->addQuads(nb_i*nb_j);
        for(int i=0; i<nb_i; i++) {
            for(int j=0; j<nb_j; j++) {
                kmds::Face f = AMesh->getFace(indexCell);
                kmds::TCellID v[4];
                v[0] = i * nb_y + j;
                v[1] = (i + 1) * nb_y + j;
                v[2] = (i + 1) * nb_y + (j + 1);
                v[3] = i * nb_y + (j + 1);
                f.setNodes(v, 4);
                indexCell++;
            }
        }

        int aId = Afp->createMaterial("matA");
        int bId = Afp->createMaterial("matB");
        int cId = Afp->createMaterial("matC");

        indexCell = 0;
        for(int i=0; i<nb_i; i++) {
            for (int j = 0; j < nb_j; j++) {

                // matA
                if((nb_i-1-i)+j <= nb_i-2) {
                    Afp->setFracPres(aId, indexCell, 1.);
                }
                if((nb_j-1-j)+i <= nb_j-2) {
                    Afp->setFracPres(aId, indexCell, 1.);
                }
                if((nb_i-1-i)+j == nb_i-3) {
                    Afp->setFracPres(aId, indexCell, 3./4.);
                }
                if((nb_j-1-j)+i == nb_j-3) {
                    Afp->setFracPres(aId, indexCell, 2./3.);
                }

                // matB
                if((i==j) && i<=nb_i/2.) {
                    Afp->setFracPres(bId, indexCell, 1.);
                }
                if(((nb_i-1-i)+j == nb_i-3) && (i<=nb_i/2.)) {
                    Afp->setFracPres(bId, indexCell, 1./4.);
                }
                if(((nb_j-1-j)+i == nb_j-3) && (i<=nb_i/2.)) {
                    Afp->setFracPres(bId, indexCell, 1./3.);
                }

                // matC
                if((i==j) && i>nb_i/2.) {
                    Afp->setFracPres(cId, indexCell, 1.);
                }
                if(((nb_i-1-i)+j == nb_i-3) && (i>nb_i/2.)) {
                    Afp->setFracPres(cId, indexCell, 1./4.);
                }
                if(((nb_j-1-j)+i == nb_j-3) && (i>nb_i/2.)) {
                    Afp->setFracPres(cId, indexCell, 1./3.);
                }

                indexCell++;
            }
        }

    }
/*----------------------------------------------------------------------------*/
    void
    initData_IxJxI_3D_nonmanifold(kmds::Mesh* AMesh, elg3d::FracPres* Afp, int nb_i, int nb_j)
    {
        const int nb_k = nb_i;

        const int nb_x = nb_i+1;
        const int nb_y = nb_j+1;
        const int nb_z = nb_k+1;

        AMesh->updateNodeCapacity(nb_x*nb_y*nb_z);
        AMesh->updateRegionCapacity(nb_i*nb_j*nb_k);

        int indexNode = AMesh->addNodes(nb_x*nb_y*nb_z);
        for(int i=0; i<nb_x; i++) {
            for(int j=0; j<nb_y; j++) {
                for(int k=0; k<nb_z; k++) {
                    AMesh->setNodeLocation(indexNode, i, j, k);
                    indexNode++;
                }
            }
        }

        kmds::TCellID indexCell = AMesh->addHexahedra(nb_i*nb_j*nb_k);
        for(int i=0; i<nb_i; i++) {
            for(int j=0; j<nb_j; j++) {
                for(int k=0; k<nb_k; k++) {
                    kmds::Region h = AMesh->getRegion(indexCell);
                    kmds::TCellID v[8];
                    v[0] = i * nb_y * nb_z + j * nb_z + k;
                    v[1] = (i + 1) * nb_y * nb_z + j * nb_z + k;
                    v[2] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k;
                    v[3] = i * nb_y * nb_z + (j + 1) * nb_z + k;
                    v[4] = i * nb_y * nb_z + j * nb_z + k +1;
                    v[5] = (i + 1) * nb_y * nb_z + j * nb_z + k +1;
                    v[6] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k +1;
                    v[7] = i * nb_y * nb_z + (j + 1) * nb_z + k + 1;
                    h.setNodes(v, 8);
                    indexCell++;
                }
            }
        }

        int aId = Afp->createMaterial("matA");
        int bId = Afp->createMaterial("matB");
        int cId = Afp->createMaterial("matC");

        indexCell = 0;
        for(int i=0; i<nb_i; i++) {
            for (int j = 0; j < nb_j; j++) {
                for (int k = 0; k < nb_k; k++) {

                    // matA
                    if((nb_i-1-i)+k <= nb_i-2) {
                        Afp->setFracPres(aId, indexCell, 1.);
                    }
                    if((nb_k-1-k)+i <= nb_k-2) {
                        Afp->setFracPres(aId, indexCell, 1.);
                    }
                    if((nb_i-1-i)+k == nb_i-3) {
                        Afp->setFracPres(aId, indexCell, 3./4.);
                    }
                    if((nb_k-1-k)+i == nb_k-3) {
                        Afp->setFracPres(aId, indexCell, 2./3.);
                    }

                    // matB
                    if((i==k) && i<=nb_i/2.) {
                        Afp->setFracPres(bId, indexCell, 1.);
                    }
                    if(((nb_i-1-i)+k == nb_i-3) && (i<=nb_i/2.)) {
                        Afp->setFracPres(bId, indexCell, 1./4.);
                    }
                    if(((nb_k-1-k)+i == nb_k-3) && (i<=nb_i/2.)) {
                        Afp->setFracPres(bId, indexCell, 1./3.);
                    }

                    // matC
                    if((i==k) && i>nb_i/2.) {
                        Afp->setFracPres(cId, indexCell, 1.);
                    }
                    if(((nb_i-1-i)+k == nb_i-3) && (i>nb_i/2.)) {
                        Afp->setFracPres(cId, indexCell, 1./4.);
                    }
                    if(((nb_k-1-k)+i == nb_k-3) && (i>nb_i/2.)) {
                        Afp->setFracPres(cId, indexCell, 1./3.);
                    }

                    indexCell++;
                }
            }
        }

    }
    /*----------------------------------------------------------------------------*/
    void
    initData_rainbow_3D_nonmanifold(kmds::Mesh* AMesh, elg3d::FracPres* Afp)
    {
        const int nb_x = 3;
        const int nb_y = 3;
        const int nb_z = 3;
        const int nb_i = nb_x - 1;
        const int nb_j = nb_y - 1;
        const int nb_k = nb_z - 1;

        AMesh->updateNodeCapacity(nb_x * nb_y * nb_z);
        AMesh->updateRegionCapacity(nb_i * nb_j * nb_k);

        int indexNode = AMesh->addNodes(nb_x * nb_y * nb_z);
        for (int i = 0; i < nb_x; i++) {
            for (int j = 0; j < nb_y; j++) {
                for (int k = 0; k < nb_z; k++) {
                    AMesh->setNodeLocation(indexNode, i, j, k);
                    indexNode++;
                }
            }
        }

        kmds::TCellID indexCell = AMesh->addHexahedra(nb_i * nb_j * nb_k);
        for (int i = 0; i < nb_i; i++) {
            for (int j = 0; j < nb_j; j++) {
                for (int k = 0; k < nb_k; k++) {
                    kmds::Region h = AMesh->getRegion(indexCell);
                    kmds::TCellID v[8];
                    v[0] = i * nb_y * nb_z + j * nb_z + k;
                    v[1] = (i + 1) * nb_y * nb_z + j * nb_z + k;
                    v[2] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k;
                    v[3] = i * nb_y * nb_z + (j + 1) * nb_z + k;
                    v[4] = i * nb_y * nb_z + j * nb_z + k + 1;
                    v[5] = (i + 1) * nb_y * nb_z + j * nb_z + k + 1;
                    v[6] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k + 1;
                    v[7] = i * nb_y * nb_z + (j + 1) * nb_z + k + 1;
                    h.setNodes(v, 8);
                    indexCell++;
                }
            }
        }

        int aId = Afp->createMaterial("matA");
        int bId = Afp->createMaterial("matB");
        int cId = Afp->createMaterial("matC");
        int dId = Afp->createMaterial("matD");
        int eId = Afp->createMaterial("matE");
        int fId = Afp->createMaterial("matF");
        int gId = Afp->createMaterial("matG");


        Afp->setFracPres(aId, 0, 1.);
        Afp->setFracPres(bId, 1, 1.);
        Afp->setFracPres(cId, 2, 1.);
        Afp->setFracPres(dId, 3, 1.);
        Afp->setFracPres(eId, 4, 1.);
        Afp->setFracPres(bId, 5, 1.);
        Afp->setFracPres(bId, 6, 1.);
        Afp->setFracPres(aId, 7, 1.);
    }
    /*----------------------------------------------------------------------------*/
    void
    initData_internalFormat(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFilename) {

        std::ifstream input(AFilename.c_str(), std::ios::in);
        if (!input)
            throw kmds::KException("initData_internalFormat Impossible to read this fracpres file");

        int dim = -1;
        input >> dim;

        if(dim == 2) {

            int nb_i = -1;
            int nb_j = -1;
            input >> nb_i >> nb_j;

            const int nb_x = nb_i + 1;
            const int nb_y = nb_j + 1;

            AMesh->updateNodeCapacity(nb_x * nb_y);
            AMesh->updateFaceCapacity(nb_i * nb_j);

            int indexNode = AMesh->addNodes(nb_x * nb_y);
            for (int i = 0; i < nb_x; i++) {
                for (int j = 0; j < nb_y; j++) {
                    AMesh->setNodeLocation(indexNode, i, j, 0.);
                    indexNode++;
                }
            }

            kmds::TCellID indexCell = AMesh->addQuads(nb_i * nb_j);
            for (int i = 0; i < nb_i; i++) {
                for (int j = 0; j < nb_j; j++) {
                    kmds::Face q = AMesh->getFace(indexCell);
                    kmds::TCellID v[4];
                    v[0] = i * nb_y + j;
                    v[1] = (i + 1) * nb_y + j;
                    v[2] = (i + 1) * nb_y + (j + 1);
                    v[3] = i * nb_y + (j + 1);
                    q.setNodes(v, 4);
                    indexCell++;
                }
            }

        } else {

            int nb_i = -1;
            int nb_j = -1;
            int nb_k = -1;
            input >> nb_i >> nb_j >> nb_k;

            const int nb_x = nb_i + 1;
            const int nb_y = nb_j + 1;
            const int nb_z = nb_k + 1;


            AMesh->updateNodeCapacity(nb_x * nb_y * nb_z);
            AMesh->updateRegionCapacity(nb_i * nb_j * nb_k);

            int indexNode = AMesh->addNodes(nb_x * nb_y * nb_z);
            for (int i = 0; i < nb_x; i++) {
                for (int j = 0; j < nb_y; j++) {
                    for (int k = 0; k < nb_z; k++) {
                        AMesh->setNodeLocation(indexNode, i, j, k);
                        indexNode++;
                    }
                }
            }

            kmds::TCellID indexCell = AMesh->addHexahedra(nb_i * nb_j * nb_k);
            for (int k = 0; k < nb_k; k++) {
                for (int j = 0; j < nb_j; j++) {
                    for (int i = 0; i < nb_i; i++) {
                        kmds::Region h = AMesh->getRegion(indexCell);
                        kmds::TCellID v[8];
                        v[0] = i * nb_y * nb_z + j * nb_z + k;
                        v[1] = (i + 1) * nb_y * nb_z + j * nb_z + k;
                        v[2] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k;
                        v[3] = i * nb_y * nb_z + (j + 1) * nb_z + k;
                        v[4] = i * nb_y * nb_z + j * nb_z + k + 1;
                        v[5] = (i + 1) * nb_y * nb_z + j * nb_z + k + 1;
                        v[6] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k + 1;
                        v[7] = i * nb_y * nb_z + (j + 1) * nb_z + k + 1;
                        h.setNodes(v, 8);
                        indexCell++;
                    }
                }
            }
        }

        int nbBlocks = -1;
        input >> nbBlocks;

        std::cout<<"nbBlocks "<<nbBlocks<<std::endl;

        for(int b = 0; b<nbBlocks; b++) {

            int nbMat = -1;
            input >> nbMat;

            std::cout<<"block "<<b<<" "<<nbMat<<std::endl;

            for(int mat=0; mat<nbMat; mat++) {

                int dummyIndex = -1;
                std::string matname("");
                input >> dummyIndex >> matname;

                if (!Afp->materialExists(matname)) {
                    Afp->createMaterial(matname);
                }

                int matIndex = Afp->getMaterialID(matname);

                int nbCells_in_mat = -1;
                input >> nbCells_in_mat;

                std::cout<<"block "<<b<<" "<<matIndex<<" "<<matname<<" "<<nbCells_in_mat<<std::endl;

                for(int c=0; c<nbCells_in_mat; c++) {
                    int index = -1;
                    double fp = -HUGE_VALF;
                    input >> index >> fp;

                    Afp->setFracPres(matIndex, index, fp);
                }
            }

        }

        input.close();
    }
    /*----------------------------------------------------------------------------*/
    void
    initData_unstructFormat(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFilename) {

        std::ifstream input(AFilename.c_str(), std::ios::in);
        if (!input) {
            throw kmds::KException("initData_unstructFormat Impossible to read this fracpres file");
        }

        int dim = -1;
        input >> dim;
        std::cout<<"dim "<<dim<<std::endl;

        int nbNodes;
        input >> nbNodes;
        std::cout<<"nbNodes "<<nbNodes<<std::endl;

        AMesh->updateNodeCapacity(nbNodes);

        int indexNode = AMesh->addNodes(nbNodes);
        for (int i = 0; i < nbNodes; i++) {
            double x, y, z;
            input >> x >> y >> z;

            AMesh->setNodeLocation(indexNode, x, y, z);
            indexNode++;

        }

        if(dim == 2) {
            int nbFaces;
            input >> nbFaces;

            AMesh->updateFaceCapacity(nbFaces);

            for (int i = 0; i < nbFaces; i++) {

                int nids;
                input >> nids;
                std::vector<kmds::TCellID> ids(nids);

                for (int j = 0; j < nids; j++) {
                    int val;
                    input >> val;
                    ids[j] = val;
                }

                switch (nids) {
                    case 3 : {
                        kmds::TCellID index = AMesh->addTriangle();
                        kmds::Face f = AMesh->getFace(index);
                        f.setNodes(ids.data(), nids);
                    }
                        break;
                    case 4 : {
                        kmds::TCellID index = AMesh->addQuad();
                        kmds::Face f = AMesh->getFace(index);
                        f.setNodes(ids.data(), nids);
                    }
                        break;
                    default:
                        throw kmds::KException("initData_unstructFormat unknown face type.");
                }

            }
        } else {
            int nbRegions;
            input >> nbRegions;

            AMesh->updateRegionCapacity(nbRegions);

            for (int i = 0; i < nbRegions; i++) {

                int nids;
                input >> nids;
                std::vector<kmds::TCellID> ids(nids);

                for (int j = 0; j < nids; j++) {
                    int val;
                    input >> val;
                    ids[j] = val;
                }

                switch (nids) {
                    case 4 : {
                        kmds::TCellID index = AMesh->addTetrahedron();
                        kmds::Region r = AMesh->getRegion(index);
                        r.setNodes(ids.data(), nids);
                    }
                        break;
                    case 8 : {
                        kmds::TCellID index = AMesh->addHexahedron();
                        kmds::Region r = AMesh->getRegion(index);
                        r.setNodes(ids.data(), nids);
                    }
                        break;
                    default:
                        throw kmds::KException("initData_unstructFormat unknown region type.");
                }

            }
        }

        int nbMat = -1;
        input >> nbMat;

        std::cout<<"nbMat "<<nbMat<<std::endl;

        for(int mat=0; mat<nbMat; mat++) {

            int dummyIndex = -1;
            std::string matname("");
            input >> dummyIndex >> matname;

            if (!Afp->materialExists(matname)) {
                Afp->createMaterial(matname);
            }

            int matIndex = Afp->getMaterialID(matname);

            int nbCells_in_mat = -1;
            input >> nbCells_in_mat;

            std::cout<<"mat "<<matIndex<<" "<<matname<<" "<<nbCells_in_mat<<std::endl;

            for(int c=0; c<nbCells_in_mat; c++) {
                int index = -1;
                double fp = -HUGE_VALF;
                input >> index >> fp;

                Afp->setFracPres(matIndex, index, fp);
            }
        }

        input.close();
    }
    /*----------------------------------------------------------------------------*/
    void
    initData_2x2x2_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp)
    {
        const int nb_x = 3;
        const int nb_y = 3;
        const int nb_z = 3;
        const int nb_i = nb_x - 1;
        const int nb_j = nb_y - 1;
        const int nb_k = nb_z - 1;

        AMesh->updateNodeCapacity(nb_x*nb_y*nb_z);
        AMesh->updateRegionCapacity(nb_i*nb_j*nb_k);

        int indexNode = AMesh->addNodes(nb_x*nb_y*nb_z);
        for(int i=0; i<nb_x; i++) {
            for(int j=0; j<nb_y; j++) {
                for(int k=0; k<nb_z; k++) {
                    AMesh->setNodeLocation(indexNode, i, j, k);
                    indexNode++;
                }
            }
        }



        kmds::TCellID indexCell = AMesh->addHexahedra(nb_i*nb_j*nb_k);
        for(int i=0; i<nb_i; i++) {
            for(int j=0; j<nb_j; j++) {
                for(int k=0; k<nb_k; k++) {
                    kmds::Region h = AMesh->getRegion(indexCell);
                    kmds::TCellID v[8];
                    v[0] = i * nb_y * nb_z + j * nb_z + k;
                    v[1] = (i + 1) * nb_y * nb_z + j * nb_z + k;
                    v[2] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k;
                    v[3] = i * nb_y * nb_z + (j + 1) * nb_z + k;
                    v[4] = i * nb_y * nb_z + j * nb_z + k +1;
                    v[5] = (i + 1) * nb_y * nb_z + j * nb_z + k +1;
                    v[6] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k +1;
                    v[7] = i * nb_y * nb_z + (j + 1) * nb_z + k + 1;
                    h.setNodes(v, 8);
                    indexCell++;
                }
            }
        }

        int aId = Afp->createMaterial("matA");
        int bId = Afp->createMaterial("matB");

        Afp->setFracPres(aId, 0, 1.);
        Afp->setFracPres(bId, 1, 1.);
        Afp->setFracPres(aId, 2, 1.);
        Afp->setFracPres(bId, 3, 1.);
        Afp->setFracPres(bId, 4, 1.);
        Afp->setFracPres(aId, 5, 1.);
        Afp->setFracPres(bId, 6, 1.);
        Afp->setFracPres(aId, 7, 1.);

    }

    /*----------------------------------------------------------------------------*/
    void
    initData_2x2x2_badpillow_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp)
    {
        const int nb_x = 3;
        const int nb_y = 3;
        const int nb_z = 3;
        const int nb_i = nb_x - 1;
        const int nb_j = nb_y - 1;
        const int nb_k = nb_z - 1;

        AMesh->updateNodeCapacity(nb_x*nb_y*nb_z);
        AMesh->updateRegionCapacity(nb_i*nb_j*nb_k);

        int indexNode = AMesh->addNodes(nb_x*nb_y*nb_z);
        for(int i=0; i<nb_x; i++) {
            for(int j=0; j<nb_y; j++) {
                for(int k=0; k<nb_z; k++) {
                    AMesh->setNodeLocation(indexNode, i, j, k);
                    indexNode++;
                }
            }
        }



        kmds::TCellID indexCell = AMesh->addHexahedra(nb_i*nb_j*nb_k);
        for(int i=0; i<nb_i; i++) {
            for(int j=0; j<nb_j; j++) {
                for(int k=0; k<nb_k; k++) {
                    kmds::Region h = AMesh->getRegion(indexCell);
                    kmds::TCellID v[8];
                    v[0] = i * nb_y * nb_z + j * nb_z + k;
                    v[1] = (i + 1) * nb_y * nb_z + j * nb_z + k;
                    v[2] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k;
                    v[3] = i * nb_y * nb_z + (j + 1) * nb_z + k;
                    v[4] = i * nb_y * nb_z + j * nb_z + k +1;
                    v[5] = (i + 1) * nb_y * nb_z + j * nb_z + k +1;
                    v[6] = (i + 1) * nb_y * nb_z + (j + 1) * nb_z + k +1;
                    v[7] = i * nb_y * nb_z + (j + 1) * nb_z + k + 1;
                    h.setNodes(v, 8);
                    indexCell++;
                }
            }
        }

        int aId = Afp->createMaterial("matA");
        int bId = Afp->createMaterial("matB");
        int cId = Afp->createMaterial("matC");
        int dId = Afp->createMaterial("matD");

        Afp->setFracPres(aId, 0, 1.);
        Afp->setFracPres(bId, 1, 1.);
        Afp->setFracPres(aId, 2, 1.);
        Afp->setFracPres(aId, 3, .8);
        Afp->setFracPres(bId, 3, .2);
//        Afp->setFracPres(aId, 3, 1.);
        Afp->setFracPres(bId, 4, 1.);
        Afp->setFracPres(bId, 5, 1.);
        Afp->setFracPres(bId, 6, 1.);
        Afp->setFracPres(aId, 7, 0.9);
        Afp->setFracPres(bId, 7, 0.08);
        Afp->setFracPres(cId, 7, 0.01);
        Afp->setFracPres(dId, 7, 0.01);

    }

/*----------------------------------------------------------------------------*/
    void
    initData_exodusReader_2D(gmds::Mesh* AMesh, std::string AFileName)
    {
#ifdef ELG3D_WITH_EXODUSII

        int CPU_word_size = 8, IO_word_size = 8;
        float version;
        int exoid = ex_open( AFileName.c_str(),
                             EX_READ,
                             &CPU_word_size,
                             &IO_word_size,
                             &version );
        if ( exoid<0 )
        {
            std::cerr<<"ERROR: Cannot open Exo File "<<AFileName<<std::endl;
            throw kmds::KException("initData_fromExodus_2D failed to open file.");
        }

        char title[MAX_LINE_LENGTH+1];
        int num_dim,
                num_nodes,
                num_elem,
                num_elem_blk,
                num_node_sets,
                num_side_sets;
        ex_get_init( exoid,
                     title,
                     &num_dim,
                     &num_nodes,
                     &num_elem,
                     &num_elem_blk,
                     &num_node_sets,
                     &num_side_sets );

        std::cout<<"reading "<<num_nodes<<" nodes, "<<num_elem<<" elements, "<<num_elem_blk<<" blocks"<<std::endl;

        double *x = new double [num_nodes];
        double *y = new double [num_nodes];
        double *z = new double [num_nodes];
        ex_get_coord( exoid, x, y, z );

        std::vector<gmds::Node> all_nodes(num_nodes);

        int *node_id_map = new int[num_nodes];
        ex_get_node_num_map(exoid, node_id_map);

        // get nodes max id
        /*
        gmds::TCellID maxID = 0;
        for ( int inode=0; inode<num_nodes; inode++ ) {
                if(node_id_map[inode] > maxID) {
                        maxID = node_id_map[inode];
                }
        }
        m_mesh.clearAndResizeNodeIDContainer(maxID);
      */
        for ( int inode=0; inode<num_nodes; inode++ ) {
            int id = node_id_map[inode];
            //gmds::Node node = m_mesh.newNodeWithID(x[inode], y[inode], z[inode],id);
            gmds::Node node = AMesh->newNode(x[inode], y[inode], z[inode]);

            all_nodes[inode] = node;
        }
        delete [] node_id_map;
        node_id_map = NULL;

        int *ids = new int[num_elem_blk];
        ex_get_elem_blk_ids( exoid, ids );

        int num_hexes = 0;
        int num_tets = 0;
        int iblock;
        for ( iblock = 0; iblock < num_elem_blk; iblock++ )
        {
            char elem_type[MAX_STR_LENGTH+1];
            int num_elems_in_block;
            int num_nodes_per_elem;
            int num_attr;
            ex_get_elem_block( exoid,
                               ids[iblock],
                               elem_type,
                               &num_elems_in_block,
                               &num_nodes_per_elem,
                               &num_attr );

            // for now assume hexes only.
            std::cout<<"reading block "<<iblock<<" of id "<<ids[iblock]<<" with "<<num_elems_in_block<<" elements of "<<num_nodes_per_elem<<" per element"<<std::endl;
            if (num_elems_in_block > 0) {
                assert( (num_nodes_per_elem == 3) || (num_nodes_per_elem == 4));
            }

            int *connect = new int [num_elems_in_block*num_nodes_per_elem];
            ex_get_elem_conn( exoid, ids[iblock], connect );

            std::string block_prefix("block_");
            std::string block_name = block_prefix + std::to_string(iblock);
            gmds::CellGroup<gmds::Face>* mat = AMesh->newGroup<gmds::Face>(block_name);

            for(int ielem=0; ielem<num_elems_in_block; ielem++) {

                gmds::Face f;

                switch(num_nodes_per_elem) {
                    case 3:
                        f = AMesh->newTriangle(
                                all_nodes[connect[num_nodes_per_elem*ielem+0]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+1]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+2]-1]
                        );
                        break;
                    case 4:
                        f = AMesh->newQuad(
                                all_nodes[connect[num_nodes_per_elem*ielem+0]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+1]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+2]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+3]-1]
                        );
                        break;
                    default:
                        throw kmds::KException("initData_fromExodus_2D unknown cell type.");
                        break;
                }

                mat->add(f);
            }

            delete[] connect;
        }


        ex_close( exoid );

        delete [] ids;
        delete [] x;
        delete [] y;
        delete [] z;

#else  // ELG3D_WITH_EXODUSII

        throw kmds::KException("initData_fromExodus_2D: exodusii not activated.");

#endif  // ELG3D_WITH_EXODUSII
    }

    /*----------------------------------------------------------------------------*/
    void
    initData_exodusReader_3D(gmds::Mesh* AMesh, std::string AFileName)
    {
#ifdef ELG3D_WITH_EXODUSII

        int CPU_word_size = 8, IO_word_size = 8;
        float version;
        int exoid = ex_open( AFileName.c_str(),
                             EX_READ,
                             &CPU_word_size,
                             &IO_word_size,
                             &version );
        if ( exoid<0 )
        {
            std::cerr<<"ERROR: Cannot open Exo File "<<AFileName<<std::endl;
            throw kmds::KException("initData_fromExodus_3D failed to open file.");
        }

        char title[MAX_LINE_LENGTH+1];
        int num_dim,
                num_nodes,
                num_elem,
                num_elem_blk,
                num_node_sets,
                num_side_sets;
        ex_get_init( exoid,
                     title,
                     &num_dim,
                     &num_nodes,
                     &num_elem,
                     &num_elem_blk,
                     &num_node_sets,
                     &num_side_sets );

        std::cout<<"reading "<<num_nodes<<" nodes, "<<num_elem<<" elements, "<<num_elem_blk<<" blocks"<<std::endl;

        double *x = new double [num_nodes];
        double *y = new double [num_nodes];
        double *z = new double [num_nodes];
        ex_get_coord( exoid, x, y, z );

        std::vector<gmds::Node> all_nodes(num_nodes);

        int *node_id_map = new int[num_nodes];
        ex_get_node_num_map(exoid, node_id_map);

        // get nodes max id
        /*
        gmds::TCellID maxID = 0;
        for ( int inode=0; inode<num_nodes; inode++ ) {
                if(node_id_map[inode] > maxID) {
                        maxID = node_id_map[inode];
                }
        }
        m_mesh.clearAndResizeNodeIDContainer(maxID);
      */
        for ( int inode=0; inode<num_nodes; inode++ ) {
            int id = node_id_map[inode];
            //gmds::Node node = m_mesh.newNodeWithID(x[inode], y[inode], z[inode],id);
            gmds::Node node = AMesh->newNode(x[inode], y[inode], z[inode]);

            all_nodes[inode] = node;
        }
        delete [] node_id_map;
        node_id_map = NULL;

        int *ids = new int[num_elem_blk];
        ex_get_elem_blk_ids( exoid, ids );

        int num_hexes = 0;
        int num_tets = 0;
        int iblock;
        for ( iblock = 0; iblock < num_elem_blk; iblock++ )
        {
            char elem_type[MAX_STR_LENGTH+1];
            int num_elems_in_block;
            int num_nodes_per_elem;
            int num_attr;
            ex_get_elem_block( exoid,
                               ids[iblock],
                               elem_type,
                               &num_elems_in_block,
                               &num_nodes_per_elem,
                               &num_attr );

            // for now assume hexes only.
            std::cout<<"reading block "<<iblock<<" of id "<<ids[iblock]<<" with "<<num_elems_in_block<<" elements of "<<num_nodes_per_elem<<" per element"<<std::endl;
            if (num_elems_in_block > 0) {
                assert( (num_nodes_per_elem == 4) || (num_nodes_per_elem == 8));
            }

            int *connect = new int [num_elems_in_block*num_nodes_per_elem];
            ex_get_elem_conn( exoid, ids[iblock], connect );

            std::string block_prefix("block_");
            std::string block_name = block_prefix + std::to_string(iblock);
            gmds::CellGroup<gmds::Region>* mat = AMesh->newGroup<gmds::Region>(block_name);

            for(int ielem=0; ielem<num_elems_in_block; ielem++) {

                gmds::Region r;

                switch(num_nodes_per_elem) {
                    case 4:
                        r = AMesh->newTet(
                                all_nodes[connect[num_nodes_per_elem*ielem+0]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+1]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+2]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+3]-1]
                        );
                        break;
                    case 8:
                        r = AMesh->newHex(
                                all_nodes[connect[num_nodes_per_elem*ielem+0]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+1]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+2]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+3]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+4]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+5]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+6]-1],
                                all_nodes[connect[num_nodes_per_elem*ielem+7]-1]
                        );
                        break;
                    default:
                        throw kmds::KException("initData_fromExodus_3D unknown cell type.");
                        break;
                }

                mat->add(r);
            }

            delete[] connect;
        }


        ex_close( exoid );

        delete [] ids;
        delete [] x;
        delete [] y;
        delete [] z;

#else  // ELG3D_WITH_EXODUSII

        throw kmds::KException("initData_fromExodus_3D: exodusii not activated.");

#endif  // ELG3D_WITH_EXODUSII
    }

/*----------------------------------------------------------------------------*/
    void
    initData_fromExodus_2D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFilename, const double Axyz_min[3], const double Axyz_max[3], int ANi, int ANj, bool ADeduceVoid)
    {
        // build the grid
        kmds::InitTools_createGrid_2D(AMesh, Axyz_min, Axyz_max, ANi, ANj);

        // read the exodusii file
        gmds::Mesh gmdsMesh(gmds::DIM2|gmds::N|gmds::F|gmds::F2N);
        initData_exodusReader_2D(&gmdsMesh, AFilename);

        // fill in the volume fractions
        elg3d::Tools_compute_fracpres_2D(AMesh, &gmdsMesh, ADeduceVoid, Afp);

    }

/*----------------------------------------------------------------------------*/
    void
    initData_fromExodus_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFilename, const double Axyz_min[3], const double Axyz_max[3], int ANi, int ANj, int ANk, bool ADeduceVoid)
    {
        // build the grid
        kmds::InitTools_createGrid_3D(AMesh, Axyz_min, Axyz_max, ANi, ANj, ANk);

        std::cout<<"aaaa"<<std::endl;

        // read the exodusii file
        gmds::Mesh gmdsMesh(gmds::DIM3|gmds::N|gmds::R|gmds::R2N);
        initData_exodusReader_3D(&gmdsMesh, AFilename);

        std::cout<<"bbbb"<<std::endl;

        // fill in the volume fractions
        elg3d::Tools_compute_fracpres_3D(AMesh, &gmdsMesh, ADeduceVoid, Afp);

    }

    /*----------------------------------------------------------------------------*/
    void
    initData_fromExodus_vf_3D(int AVoidID, kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFileName)
    {
#ifdef ELG3D_WITH_EXODUSII

        int CPU_word_size = 8, IO_word_size = 8;
        float version;
        int exoid = ex_open( AFileName.c_str(),
                             EX_READ,
                             &CPU_word_size,
                             &IO_word_size,
                             &version );
        if ( exoid<0 )
        {
            std::cerr<<"initData_fromExodus_vf_3D Cannot open Exo File "<<AFileName<<std::endl;
            throw kmds::KException("initData_fromExodus_vf_3D Cannot open Exo File");
        }

        char title[MAX_LINE_LENGTH+1];
        int num_dim,
                num_nodes,
                num_elem,
                num_elem_blk,
                num_node_sets,
                num_side_sets;
        ex_get_init( exoid,
                     title,
                     &num_dim,
                     &num_nodes,
                     &num_elem,
                     &num_elem_blk,
                     &num_node_sets,
                     &num_side_sets );

        std::cout<<"reading "<<num_nodes<<" nodes, "<<num_elem<<" elements, "<<num_elem_blk<<" blocks"<<std::endl;
/*
        if(1 != num_elem_blk) {
                std::cerr<<"MeshFracPres::readFromFileExodus number of blocks is not one!"<<AFileName<<std::endl;
                throw GMDSException("MeshFracPres::readFromFileExodus number of blocks is not one!");
        }
*/

        // get element variables
        int num_ele_vars = 0;
        ex_get_var_param( exoid, "e", &num_ele_vars);
        std::cout<<"number of variables "<< num_ele_vars<<std::endl;
        if( num_ele_vars <= 0 ) {
            std::cerr<<"number of variable is zero"<<std::endl;
            exit(-1);
        }
        for(int imat=0; imat<num_ele_vars; imat++) {
            // TODO : what to do when the void material is explicitly stored ?
            if(imat == AVoidID ) continue; // we do not read the void
            std::string mat_name = std::string("mat_") + std::to_string(imat);
            Afp->createMaterial(mat_name);
        }

        double *x = new double [num_nodes];
        double *y = new double [num_nodes];
        double *z = new double [num_nodes];
        ex_get_coord( exoid, x, y, z );

        std::vector<kmds::TCellID > all_nodes(num_nodes);

        int *node_id_map = new int[num_nodes];
        ex_get_node_num_map(exoid, node_id_map);

        // get nodes max id
        /*
        gmds::TCellID maxID = 0;
        for ( int inode=0; inode<num_nodes; inode++ ) {
                if(node_id_map[inode] > maxID) {
                        maxID = node_id_map[inode];
                }
        }
        m_mesh.clearAndResizeNodeIDContainer(maxID);
      */

        AMesh->updateNodeCapacity(num_nodes);
        AMesh->addNodes(num_nodes);

        for ( int inode=0; inode<num_nodes; inode++ ) {
            int id = node_id_map[inode];
            //gmds::Node node = m_mesh.newNodeWithID(x[inode], y[inode], z[inode],id);
            AMesh->setNodeLocation(inode, x[inode], y[inode], z[inode]);

            all_nodes[inode] = inode;
        }
        delete [] node_id_map;
        delete [] x;
        delete [] y;
        delete [] z;
        node_id_map = NULL;

        int *ids = new int[num_elem_blk];
        ex_get_elem_blk_ids( exoid, ids );

        int num_hexes = 0;
        int num_tets = 0;
        int* num_elems_in_block = new int[num_elem_blk];
        int iblock;
        for ( iblock = 0; iblock < num_elem_blk; iblock++ )
        {
            char elem_type[MAX_STR_LENGTH+1];
            //int num_elems_in_block;
            int num_nodes_per_elem;
            int num_attr;
            ex_get_elem_block( exoid,
                               ids[iblock],
                               elem_type,
                               &num_elems_in_block[iblock],
                               &num_nodes_per_elem,
                               &num_attr );

            // for now assume hexes only.
            std::cout<<"reading block "<<iblock<<" of id "<<ids[iblock]<<" with "<<num_elems_in_block[iblock]<<" elements of "<<num_nodes_per_elem<<" per element"<<" and "<<num_attr<<" attributes"<<std::endl;

            if (num_elems_in_block[iblock] > 0) {
                assert( num_nodes_per_elem == 8);
            }

            int *connect = new int [num_elems_in_block[iblock]*num_nodes_per_elem];
            ex_get_elem_conn( exoid, ids[iblock], connect );

            AMesh->updateRegionCapacity(AMesh->getNbRegions() + num_elems_in_block[iblock]);

            for(int ihex=0; ihex<num_elems_in_block[iblock]; ihex++) {

                AMesh->newHexahedron(
                        all_nodes[connect[8*ihex+0]-1],
                        all_nodes[connect[8*ihex+1]-1],
                        all_nodes[connect[8*ihex+2]-1],
                        all_nodes[connect[8*ihex+3]-1],
                        all_nodes[connect[8*ihex+4]-1],
                        all_nodes[connect[8*ihex+5]-1],
                        all_nodes[connect[8*ihex+6]-1],
                        all_nodes[connect[8*ihex+7]-1]
                );

            }

            double** elem_var_vals;
            elem_var_vals = (double **)calloc(num_ele_vars, sizeof (double *));

            for (int imat=0; imat < num_ele_vars; imat++) {

                elem_var_vals[imat] = (double *)calloc(num_elems_in_block[iblock], sizeof(double));

                if (ex_get_var(exoid, 1, EX_ELEM_BLOCK, imat+1, ids[iblock], num_elems_in_block[iblock],

                               elem_var_vals[imat]) < 0) {

                }

            }

            for(int imat=0; imat<num_ele_vars; imat++) {
                // TODO : void materials will not be stored
                if(imat == AVoidID) continue;
                int numMat = (imat<AVoidID) ? imat : imat-1;
                for(int ihex=0; ihex<num_elems_in_block[iblock]; ihex++) {
                    // TODO: here region id is ihex but it could be different
                    Afp->setFracPres(imat-1, ihex, elem_var_vals[imat][ihex]);
                }
            }

            delete[] connect;
        }


        ex_close( exoid );

        delete [] num_elems_in_block;

        delete [] ids;
        delete [] x;
        delete [] y;
        delete [] z;

#else  // ELG3D_WITH_EXODUSII

        throw kmds::KException("initData_fromExodus_3D: exodusii not activated.");

#endif  // ELG3D_WITH_EXODUSII
    }

    /*----------------------------------------------------------------------------*/
    void
    initData_TEClike_2D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFilename) {

        std::ifstream input(AFilename.c_str(), std::ios::in);
        if (!input)
            throw kmds::KException("initData_TEClike_2D Impossible to read this fracpres file");


        int nb_i = -1;
        int nb_j = -1;
//        input >> nb_i >> nb_j;

        std::string line;
        std::getline(input, line);

        std::cout<<"line " <<line<<std::endl;

        std::size_t found_i = line.find("nb_i=");
        if (found_i!=std::string::npos) {
            std::size_t found_end = line.find(",", found_i+5);

            if (found_end!=std::string::npos) {
                nb_i = std::stoi(line.substr(found_i+5, found_end));
                std::cout<<"nb_i "<<nb_i<<std::endl;
            } else {
                std::cout<<"coma not found"<<std::endl;
                exit(-1);
            }
        } else {
            std::cout<<"nb_i not found"<<std::endl;
            exit(-1);
        }

        std::size_t found_j = line.find("nb_j=");
        if (found_j!=std::string::npos) {

                nb_j = std::stoi(line.substr(found_j+5, line.length()));
                std::cout<<"nb_j "<<nb_j<<std::endl;
        } else {
            std::cout<<"nb_j not found"<<std::endl;
            exit(-1);
        }


//        nb_i = 750;
//        nb_j = 750;

        const int nb_x = nb_i + 1;
        const int nb_y = nb_j + 1;

        AMesh->updateNodeCapacity(nb_x * nb_y);
        AMesh->updateFaceCapacity(nb_i * nb_j);


        double xyz_min[3] = {0., 0., 0.};
        double xyz_max[3] = {nb_x, nb_y, 0.};
        kmds::InitTools_createGrid_2D(AMesh, xyz_min, xyz_max, nb_i, nb_j);

        const int nbCells = nb_i * nb_j;

        const int imatA = Afp->createMaterial("matA");
        const int imatB = Afp->createMaterial("matB");

        for(int icell=0; icell<nbCells; icell++) {
            double x;
            double y;
            double vf;

            input >> x >> y >> vf;

//            int index_i = std::rint(x / .002);
//            int index_j = std::rint(y / .002);
            int index_i = std::rint(x / .0005);
            int index_j = std::rint(y / .0005);


            int index = index_i * nb_j + index_j;

//            std::cout<<"index "<<index_i<<" "<<index_j<<" "<<index<<std::endl;

            Afp->setFracPres(imatA, index, vf);
            Afp->setFracPres(imatB, index, 1. - vf);
        }

        input.close();

        std::cout<<"nbMat "<<Afp->getNbMaterials()<<std::endl;
}

/*----------------------------------------------------------------------------*/
    void
    initData_TEClike_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFilename) {

        std::ifstream input(AFilename.c_str(), std::ios::in);
        if (!input)
            throw kmds::KException("initData_TEClike_3D Impossible to read this fracpres file");


        int nb_i = -1;
        int nb_j = -1;
        int nb_k = -1;
//        input >> nb_i >> nb_j;

        std::string line;
        std::getline(input, line);

        std::cout<<"line " <<line<<std::endl;

        std::size_t found_i = line.find("nb_i=");
        if (found_i!=std::string::npos) {
            std::size_t found_end = line.find(",", found_i+5);

            if (found_end!=std::string::npos) {
                nb_i = std::stoi(line.substr(found_i+5, found_end));
                std::cout<<"nb_i "<<nb_i<<std::endl;
            } else {
                std::cout<<"coma not found"<<std::endl;
                exit(-1);
            }
        } else {
            std::cout<<"nb_i not found"<<std::endl;
            exit(-1);
        }

        std::size_t found_j = line.find("nb_j=");
        if (found_j!=std::string::npos) {

            nb_j = std::stoi(line.substr(found_j+5, line.length()));
            std::cout<<"nb_j "<<nb_j<<std::endl;
        } else {
            std::cout<<"nb_j not found"<<std::endl;
            exit(-1);
        }

        nb_k = nb_i;
        std::cout<<"nb_k "<<nb_k<<std::endl;

//        nb_i = 750;
//        nb_j = 750;

        const int nb_x = nb_i + 1;
        const int nb_y = nb_j + 1;
        const int nb_z = nb_k + 1;

        AMesh->updateNodeCapacity(nb_x * nb_y * nb_z);
        AMesh->updateRegionCapacity(nb_i * nb_j * nb_k);


        double xyz_min[3] = {0., 0., 0.};
        double xyz_max[3] = {nb_x, nb_y, nb_z};
        kmds::InitTools_createGrid_3D(AMesh, xyz_min, xyz_max, nb_i, nb_j, nb_k);

        const int nbCells = nb_i * nb_j * nb_k;

        const int imatA = Afp->createMaterial("matA");
        const int imatB = Afp->createMaterial("matB");

        for(int icell=0; icell<nbCells; icell++) {
            double x;
            double y;
            double z;
            double vf;

            input >> x >> y >> z >> vf;

//            int index_i = std::rint(x / .002);
//            int index_j = std::rint(y / .002);
//            int index_k = std::rint(z / .002);
            int index_i = std::rint(x / .008);
            int index_j = std::rint(y / .008);
            int index_k = std::rint(z / .008);


            int index = index_i * nb_j * nb_k + index_j * nb_k + index_k;

            if(icell % 100000 == 0) {
                std::cout << "index " << index_i << " " << index_j << " " << index_k << " " << index << std::endl;
            }

            Afp->setFracPres(imatA, index, vf);
            Afp->setFracPres(imatB, index, 1. - vf);
        }

        input.close();

        std::cout<<"nbMat "<<Afp->getNbMaterials()<<std::endl;
    }
/*----------------------------------------------------------------------------*/
    void
    initData_fromTXT_3D(kmds::Mesh* AMesh, elg3d::FracPres* Afp, std::string AFilename, const double Axyz_min[3], const double Axyz_max[3], int ANi, int ANj, int ANk, int ANbMat)
    {
        // build the grid
        kmds::InitTools_createGrid_3D(AMesh, Axyz_min, Axyz_max, ANi, ANj, ANk);

        for(int imat=0; imat<ANbMat; imat++) {
            std::string mat_name = std::string("mat") + std::to_string(imat);
            Afp->createMaterial(mat_name);
        }

        const int nbCells = AMesh->getNbRegions();

        // open the file
        std::ifstream myfile;
        myfile.open (AFilename);
        for(int i=0; i<nbCells; i++) {

            for(int imat=0; imat<ANbMat; imat++) {
                double scalar;
                myfile >> scalar;
                if(scalar > 0.) {
                    Afp->setFracPres(imat, i, scalar);
                }
            }
        }
        myfile.close();

        // fill in the volume fractions


    }
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
