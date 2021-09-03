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
/*
 * InitTools.cpp
 *
 *  Created on: 15 mars 2019
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include "KM/Utils/InitTools.h"
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
    void
    InitTools_createGrid_2D(kmds::Mesh* AMesh, const double xyz_min[3], const double xyz_max[3], int ni, int nj)
    {
        const int nb_x = ni + 1;
        const int nb_y = nj + 1;

        AMesh->updateNodeCapacity(nb_x * nb_y);
        AMesh->updateFaceCapacity(ni * nj);

        int indexNode = AMesh->addNodes(nb_x * nb_y);
        for (int i = 0; i < nb_x; i++) {
            for (int j = 0; j < nb_y; j++) {

                double x = xyz_min[0] + ( (double) i / (double) ni ) * (xyz_max[0] - xyz_min[0]);
                double y = xyz_min[1] + ( (double) j / (double) nj ) * (xyz_max[1] - xyz_min[1]);

                AMesh->setNodeLocation(indexNode, x, y, 0.);
                indexNode++;
            }
        }


        kmds::TCellID indexCell = AMesh->addQuads(ni * nj);
        for (int i = 0; i < ni; i++) {
            for (int j = 0; j < nj; j++) {
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

    }

/*----------------------------------------------------------------------------*/
    void
    InitTools_createGrid_3D(kmds::Mesh* AMesh, const double xyz_min[3], const double xyz_max[3], int ni, int nj, int nk)
    {
        const int nb_x = ni + 1;
        const int nb_y = nj + 1;
        const int nb_z = nk + 1;

        AMesh->updateNodeCapacity(nb_x * nb_y * nb_z);
        AMesh->updateRegionCapacity(ni * nj * nk);

        int indexNode = AMesh->addNodes(nb_x * nb_y * nb_z);
        for (int i = 0; i < nb_x; i++) {
            for (int j = 0; j < nb_y; j++) {
                for (int k = 0; k < nb_z; k++) {

                    double x = xyz_min[0] + ((double) i / (double)
                    ni) *(xyz_max[0] - xyz_min[0]);
                    double y = xyz_min[1] + ((double) j / (double)
                    nj) *(xyz_max[1] - xyz_min[1]);
                    double z = xyz_min[2] + ((double) k / (double)
                    nk) *(xyz_max[2] - xyz_min[2]);

                    AMesh->setNodeLocation(indexNode, x, y, z);
                    indexNode++;
                }
            }
        }

        kmds::TCellID indexCell = AMesh->addHexahedra(ni * nj * nk);
        for (int i = 0; i < ni; i++) {
            for (int j = 0; j < nj; j++) {
                for (int k = 0; k < nk; k++) {
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

/*----------------------------------------------------------------------------*/
}  // namespace kmds
/*----------------------------------------------------------------------------*/
