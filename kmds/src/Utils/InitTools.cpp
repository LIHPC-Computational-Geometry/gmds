/*----------------------------------------------------------------------------*/
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
