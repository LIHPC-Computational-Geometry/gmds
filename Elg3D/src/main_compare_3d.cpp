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
/** \file    main_pixelVF3d.cpp
 *  \author  legoff
 *  \date    11/21/2018
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/IO/VTKWriter.h>

#include <gmds/cad/FACManager.h>
#include <KM/DS/Mesh.h>

#include <gmds/ig/Mesh.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/IGMeshIOService.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/ALGOCMPNT/ExtractGeomModel.h>
#include <ELG3D/ALGOCMPNT/InitData.h>
#include <ELG3D/ALGOCMPNT/InterfaceNodesPosSmoothVF.h>
#include <ELG3D/ALGOCMPNT/Tools.h>
#include <ELG3D/DATACMPNT/FracPres.h>
#include <ELG3D/DATACMPNT/MaterialAssignment.h>


#include <ELG3D/ALGOCMPNT/AssignCells.h>
#include <ELG3D/ALGOCMPNT/BoundingBoxGeomAssociation.h>
#include <ELG3D/ALGOCMPNT/InterfaceNodesPos.h>
#include <ELG3D/ALGOCMPNT/ManifoldDetection.h>
#include <ELG3D/ALGOCMPNT/MaterialGradientComputation.h>
#include <ELG3D/ALGOCMPNT/MaterialInterfaces.h>
#include <ELG3D/ALGOCMPNT/MeshExtractor.h>
#include <ELG3D/ALGOCMPNT/MoveToNewPos.h>
#include <ELG3D/ALGOCMPNT/Pillow.h>
#include <ELG3D/ALGOCMPNT/SmartLaplacian.h>
#include <ELG3D/DATACMPNT/FacetedCurveGeomServices.h>
#include <ELG3D/DATACMPNT/FacetedSurfaceGeomServices.h>
#include <ELG3D/DATACMPNT/Parameters.h>
/*----------------------------------------------------------------------------*/


int main(int argc, char* argv[]) {
    int num_threads = -1;

    // nbthreads writeoutput fromfile filename deducevoid

    std::istringstream iss(argv[1]);
    int val;

    if (iss >> val) {
        num_threads = val;
    } else {
        std::cerr << "could not convert number of threads argument." << std::endl;
        exit(-1);
    }

    bool writeOutputFiles = false;
    std::istringstream iss_2(argv[2]);
    iss_2 >> writeOutputFiles;

    bool isFromFile = false;
    std::istringstream iss_3(argv[3]);
    iss_3 >> isFromFile;

    std::string case_filename_fp;
    std::string case_filename_ma;
    if (isFromFile) {
        std::istringstream iss_4(argv[4]);
        iss_4 >> case_filename_fp;

        std::istringstream iss_5(argv[5]);
        iss_5 >> case_filename_ma;
    }

    Kokkos::InitArguments kargs;
    kargs.num_threads = num_threads;
    Kokkos::initialize(kargs);

    std::cout<<"num_threads "<<num_threads<<std::endl;

    Kokkos::Timer timer_tot;
    timer_tot.reset();

    Kokkos::Timer timer;
    timer.reset();


    elg3d::Parameters::dimension = 3;
    elg3d::Parameters::celltype = kmds::KMDS_REGION;

    kmds::Mesh mesh_fp;
    elg3d::FracPres fp_ref;

    if (isFromFile) {
//        elg3d::initData_internalFormat(&mesh_source, &fp_source, case_filename);
        elg3d::Tools_read_fracpres_vtk_3D(&mesh_fp, &fp_ref, case_filename_fp, false);
//        elg3d::initData_TEClike_3D(&mesh_source, &fp_source, case_filename);

    } else {
        const double xyz_min[3] = {0., 0., 0.};
        const double xyz_max[3] = {16., 8., 8.};
        elg3d::initData_fromTXT_3D(&mesh_fp, &fp_ref, "volfrac.txt", xyz_min, xyz_max, 120, 60, 60, 2);
//        std::cerr<<"needs to be read from a file."<<std::endl;
    }

    std::cout<<"nbcells "<<mesh_fp.getNbRegions()<<std::endl;
    std::cout<<"fp.nbMat "<<fp_ref.getNbMaterials()<<std::endl;


    kmds::Mesh mesh_ma;
    elg3d::MaterialAssignment ma;
    ma.createMaterials(fp_ref.getMaterialList());

    elg3d::Tools_read_ma_vtk_3D(&mesh_ma, &ma, case_filename_ma);

    elg3d::FracPres fp_new;
//    elg3d::Tools_compute_fracpres_3D(&mesh_fp, &mesh_ma, &ma, &fp_new);
    elg3d::Tools_compute_fracpres_3D(&mesh_ma, &mesh_fp, &ma, &fp_new);

    kmds::Variable<double>* varDiscrepancy = mesh_fp.createVariable<double>(-HUGE_VALF, kmds::KMDS_REGION, "discrepancy");
    elg3d::Tools_compute_discrepancy_3D(&mesh_fp, &fp_ref, &fp_new, varDiscrepancy);

    double minJacobian = elg3d::Tools_computeScaledJacobian_3D(&mesh_ma);

    const kmds::TCellID nbCells_ref = mesh_fp.getNbRegions();
    double discrepancy = 0;
    Kokkos::parallel_reduce(nbCells_ref,
                            KOKKOS_LAMBDA(const int i, double &sum) {
                                sum += (*varDiscrepancy)[i];

                            },
                            discrepancy);

    std::cout<<"sumdiscrepancy "<<discrepancy<<std::endl;
    std::cout<<"minJacobian "<<minJacobian<<std::endl;

    elg3d::Tools_write_3D(&mesh_fp, &fp_ref, &ma, "discrepancy");
    //==========================================================================


    Kokkos::finalize();
    return 0;
}