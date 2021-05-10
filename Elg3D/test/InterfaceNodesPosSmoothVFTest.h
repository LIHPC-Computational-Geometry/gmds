/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/InterfaceNodesPosSmoothVF.h"

#include <ELG3D/ALGOCMPNT/Tools.h>

#include "KM/Utils/InitTools.h"
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
class InterfaceNodesPosSmoothVFTest : public ::testing::Test
{
protected:
    InterfaceNodesPosSmoothVFTest()
    {
        ;
    }
    virtual ~InterfaceNodesPosSmoothVFTest()
    {
        ;
    }

    static void
    SetUpTestCase()
    {
        // Kokkos::Serial::initialize();
        // Kokkos::Threads::initialize();
        Kokkos::InitArguments kargs;
        kargs.num_threads = 1;
//                int num_threads = 4;
//                int use_numa = 1;
//                int use_core = 1;
//                Kokkos::OpenMP::initialize(num_threads, use_numa, use_core);
        Kokkos::initialize(kargs);
    }

    static void
    TearDownTestCase()
    {
        // Kokkos::Serial::finalize();
        // Kokkos::Threads::finalize();
//                Kokkos::OpenMP::finalize();
        Kokkos::finalize();
    }
};
/*----------------------------------------------------------------------------*/
TEST_F(InterfaceNodesPosSmoothVFTest, dummy)
{
    EXPECT_EQ(1, 1);
}
/*----------------------------------------------------------------------------*/
TEST_F(InterfaceNodesPosSmoothVFTest, smallGraph)
{
    const int nsub = 10;

    const int nbVert = nsub * nsub + 3 * nsub;
    const int nbColor = 2;

    const double fpA = 0.7;
    const double fpB = 0.3;

    // init (i,j) -> index
    int indices[nsub+2][nsub+2];

    for(int i=0; i<nsub+2; i++) {
        for(int j=0; j<nsub+2; j++) {
            indices[i][j] = -1;
        }
    }
    int index = 0;
    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {
            indices[i][j] = index;
            index++;
        }
    }
    // left
    for(int j=1; j<=nsub; j++) {
        indices[0][j] = index;
        index++;
    }
    //right
    for(int j=1; j<=nsub; j++) {
        indices[nsub+1][j] = index;
        index++;
    }
    //top
    for(int i=1; i<=nsub; i++) {
        indices[i][nsub+1] = index;
        index++;
    }

    // graph builders
    const int maxNbNeighbors = 8;
    kmds::Graph gr("graph", nbVert, maxNbNeighbors);

    for(int i=0; i<nsub+2; i++) {
        for(int j=0; j<nsub+2; j++) {
            if(indices[i][j] != -1) {
                std::set<kmds::TCellID> nids;

                // list neighbors
                for(int di=-1; di<=1; di++) {
                    for(int dj=-1; dj<=1; dj++) {
                        if((di != 0) && (dj !=0)) {
                            if((0 <=i + di ) && (i + di <= nsub+1)) {
                                if((0 <=j + dj ) && (j + dj <= nsub+1)) {
                                    if(indices[i+di][j+dj] != -1) {
                                        nids.insert(indices[i+di][j+dj]);
                                    }
                                }
                            }
                        }
                    }
                }

                gr.setNeighbors(indices[i][j], nids);
            }
        }
    }


    // initialize vf
    int nbFixed = 0;
    std::vector<bool> isfixed(nbVert, 0);
    std::vector<int> mat(nbVert, -1);
    std::vector<gmds::math::Vector> vf_current(nbVert, gmds::math::Vector(fpA,fpB,0.));
    std::vector<gmds::math::Vector> vf_next(nbVert, gmds::math::Vector(0.,0.,0.));

    // left
    for(int j=1; j<=nsub; j++) {
        int index = indices[0][j];
        vf_current[index] = gmds::math::Vector(1.,0.,0.);
    }
    //right
    for(int j=1; j<=nsub; j++) {
        int index = indices[nsub+1][j];
        vf_current[index] = gmds::math::Vector(0.,1.,0.);
    }
    //top
    for(int i=1; i<=nsub; i++) {
        int index = indices[i][nsub+1];
        vf_current[index] = gmds::math::Vector(1.,0.,0.);
    }

    double threshold = 1.;
    int nbOuterIter = 0;
    while(nbFixed != nbVert) {

        std::cout<<"nbFixed "<<nbFixed<<" threshold "<<threshold<<std::endl;

        // fix pixels above threshold
        int nbFixed_toadd = 0;
        for(int p=0; p<nbVert; p++) {
            if(!isfixed[p]) {
                if(vf_current[p].get(0) >= threshold) {
                    isfixed[p] = true;
                    mat[p] = 0;
                    nbFixed_toadd++;
                    vf_current[p].set(0,1.);
                    vf_current[p].set(1,0.);
                } else {
                    if(vf_current[p].get(1) >= threshold) {
                        isfixed[p] = true;
                        mat[p] = 1;
                        nbFixed_toadd++;
                        vf_current[p].set(1,1.);
                        vf_current[p].set(0,0.);
                    }
                }
            }
        }
        nbFixed += nbFixed_toadd;

        // update threshold if no pixel was fixed this round
        if(nbFixed_toadd == 0) {
            threshold -= 0.01;
        }


        // update vf considering already fixed pixels
        int nblabel_A = 0;
        int nblabel_B = 0;

        for(int i=1; i<=nsub; i++) {
            for(int j=1; j<=nsub; j++) {
                int p = indices[i][j];

                if (isfixed[p]) {
                    if(mat[p] == 0) {
                        nblabel_A++;
                    } else {
                        nblabel_B++;
                    }
                }
            }
        }

        // all subpixels are marked
        if(nsub*nsub - nblabel_A - nblabel_B == 0) {
            break;
        }

        double target_A = nsub*nsub*fpA - nblabel_A;
        if(target_A < 0.) {
            target_A = 0.;
        }

        double fpA_next = target_A / (nsub*nsub - nblabel_A - nblabel_B);
        double fpB_next = 1. - fpA_next;

        std::cout<<"nblabel_A "<<nblabel_A<<" nblabel_B "<<nblabel_B<<std::endl;
        std::cout<<fpA_next<<" "<<fpB_next<<std::endl;

        for(int i=1; i<=nsub; i++) {
            for(int j=1; j<=nsub; j++) {
                int p = indices[i][j];

                if (!isfixed[p]) {
                    vf_current[p] = gmds::math::Vector (fpA_next, fpB_next, 0.);
                }
            }
        }

        // smooth vf
        for(int iter=0; iter<10; iter++) {
            for (int p = 0; p < nbVert; p++) {
                if (!isfixed[p]) {

                    gmds::math::Vector v(vf_current[p]);

                    // do not smooth if one component is higher than threshold
                    if((v.get(0) >= threshold) || (v.get(1) >= threshold)) {
                        vf_next[p] = v;
                        break;
                    }

                    Kokkos::View<kmds::TCellID*> nids;
                    const int nbNeighbors = gr.getNeighbors(p, nids);

                    for(int i=0; i<nbNeighbors; i++) {
                        v = v + vf_current[nids[i]];
                    }

                    v = 1. / (nbNeighbors + 1) * v;

                    // no need for normalization, the norm1(average) = 1
//                    double norm1 = v.get(0) + v.get(1);
//                    vf_next[p] = 1./norm1 * v;
                    vf_next[p] = v;
                }
            }

            for (int p = 0; p < nbVert; p++) {
                if (!isfixed[p]) {
                    vf_current[p] = vf_next[p];
                }
            }

        }

        nbOuterIter++;
    }

    std::cout<<"=============="<<std::endl;

    for(int j=nsub+1; j>=0; j--) {
        for(int i=0; i<nsub+2; i++) {

            if(indices[i][j] == -1) {
                std::cout<<" ";
            } else {
                std::cout<<mat[indices[i][j]];
            }
        }
        std::cout<<std::endl;
    }
    std::cout<<"=============="<<std::endl;

    // cost fp-wise
    int nblabel_A = 0;
    int nblabel_B = 0;

    for(int i=1; i<=nsub; i++) {
        for(int j=1; j<=nsub; j++) {

            if(mat[indices[i][j]] == 0) {
                nblabel_A++;
            }
            if(mat[indices[i][j]] == 1) {
                nblabel_B++;
            }
        }
    }
    std::cout<<"fpA "<<(double)nblabel_A/(nsub*nsub)<<" fpB "<<(double)nblabel_B/(nsub*nsub)<<std::endl;

}
/*----------------------------------------------------------------------------*/
TEST_F(InterfaceNodesPosSmoothVFTest, heuristic_5x5_2D) {

    kmds::Mesh mesh_source;
    elg3d::FracPres fp_source;

    elg3d::initData_internalFormat(&mesh_source, &fp_source, "/home/legoffn/travail/GMDS/gmds_build_clion/Elg3D/test/Samples/internalFormat_5x5_2D.txt");


    kmds::VTKWriter<kmds::Mesh> w(mesh_source);
    w.write("heuristic_5x5_source", kmds::F);

    kmds::Mesh mesh_target;
    elg3d::MaterialAssignment ma_target;

    kmds::Connectivity* c_N2F = mesh_source.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh_source);
    ch.buildN2F();

    kmds::Connectivity* c_E2F = mesh_source.createConnectivity(kmds::E2F);
    ch.buildEandE2F_2D_variant_1();


    kmds::Variable<kmds::TCellID> *varSubCells2oldCells =
            mesh_target.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "subCells2oldCells");
    kmds::Variable<kmds::TCellID> *varOldCells2firstSubCells =
            mesh_source.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "oldCells2firstSubCells");

    kmds::Variable<bool> *varMixedCells_source =
            mesh_source.createVariable<bool>(false, kmds::KMDS_FACE, "mixedCells_source");
    kmds::Variable<bool> *varMixedCells_target =
            mesh_target.createVariable<bool>(false, kmds::KMDS_FACE, "mixedCells_target");

    kmds::Variable<double> *varSurfvol_source =
            mesh_source.createVariable<double>(false, kmds::KMDS_FACE, "surfvol_source");
    kmds::Variable<double> *varSurfvol_target =
            mesh_target.createVariable<double>(false, kmds::KMDS_FACE, "surfvol_target");


    elg3d::InterfaceNodesPosSmoothVF_buildSubMesh_2D(&mesh_source,
                                                     &mesh_target,
                                                     c_N2F,
                                                     c_E2F,
                                                     &fp_source,
                                                     &ma_target,
                                                     varSubCells2oldCells,
                                                     varOldCells2firstSubCells,
                                                     varMixedCells_source,
                                                     varMixedCells_target,
                                                     varSurfvol_source,
                                                     varSurfvol_target
    );

    elg3d::Tools_write_lite_2D(&mesh_target, &ma_target, "heuristic_5x5_before");

    kmds::Connectivity* c_N2F_target = mesh_target.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch_target(&mesh_target);
    ch_target.buildN2F();
    kmds::Connectivity* c_F2F_byN_target = mesh_target.createConnectivity(kmds::F2F_byN);
    ch_target.buildF2F_byN(c_N2F_target);

    elg3d::FracPres fp_target;


    // build graph
    kmds::Variable<kmds::TCellID> *varCells2vertices =
            mesh_target.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "cell2vertices");

    kmds::GrowingView<kmds::TCellID> cellsIDs_target("CELLS", mesh_target.getNbFaces());
    mesh_target.getFaceIDs(&cellsIDs_target);

    kmds::GrowingView<kmds::TCellID> selectionCells("CELLS", mesh_target.getNbFaces());

    elg3d::InterfaceNodesPosSmoothVF_selectCells_xD(&cellsIDs_target,
                                                    c_F2F_byN_target,
                                                    varMixedCells_target,
                                                    varCells2vertices,
                                                    &selectionCells);

    kmds::Graph graph("CELLS_SUBMESH_GRAPH", selectionCells.getNbElems(), 20);
    elg3d::InterfaceNodesPosSmoothVF_buildGraph_xD(&selectionCells,
                                                   c_F2F_byN_target,
                                                   varCells2vertices,
                                                   &graph);


    kmds::GrowingView<kmds::TCellID> cellsIDs_source("CELLS", mesh_source.getNbFaces());
    mesh_source.getFaceIDs(&cellsIDs_source);


    const int nbSubPixels = elg3d::InterfaceNodesPosSmoothVF_NBSUB2;


    elg3d::InterfaceNodesPosSmoothVF_assignHeuristic_xD(&cellsIDs_source,
                                                        &cellsIDs_target,
                                                        &fp_source,
                                                        &ma_target,
                                                        varSubCells2oldCells,
                                                        varOldCells2firstSubCells,
                                                        varMixedCells_source,
                                                        varCells2vertices,
                                                        varSurfvol_source,
                                                        varSurfvol_target,
                                                        &graph,
                                                        nbSubPixels
    );


    elg3d::Tools_write_2D(&mesh_target, &fp_target, &ma_target, "heuristic_5x5_after");
}
/*----------------------------------------------------------------------------*/
TEST_F(InterfaceNodesPosSmoothVFTest, heuristic_3x1x1_3D) {

    kmds::Mesh mesh_source;
    elg3d::FracPres fp_source;

    const int ni = 3;
    const int nj = 1;
    const int nk = 1;

    double xyz_min[3] = {0.,0.,0.};
    double xyz_max[3] = {3.,1.,1.};
    kmds::InitTools_createGrid_3D(&mesh_source, xyz_min, xyz_max, ni, nj, nk);

    kmds::Region r = mesh_source.getRegion(1);
    Kokkos::View<kmds::TCellID *> nids;
    r.nodeIds(nids);
    for(int i=0; i<8; i++) {
        std::cout<<"node "<<nids[i]<<std::endl;
    }
    nids[0] = 6;
    nids[1] = 4;
    nids[2] = 8;
    nids[3] = 10;
    nids[4] = 7;
    nids[5] = 5;
    nids[6] = 9;
    nids[7] = 11;

    fp_source.createMaterial("matA");
    fp_source.createMaterial("matB");

    fp_source.setFracPres(0, 0, 1.);
    fp_source.setFracPres(0, 1, 0.5);
    fp_source.setFracPres(1, 1, 0.5);
    fp_source.setFracPres(1, 2, 1.);

    elg3d::MaterialAssignment ma_source;
    ma_source.createMaterial("matA");
    ma_source.createMaterial("matB");
    ma_source.setMaterial(0,0);
    ma_source.setMaterial(0,1);
    ma_source.setMaterial(0,2);
    elg3d::Tools_write_3D(&mesh_source, &fp_source, &ma_source, "heuristic_3x1x1_source");

//    kmds::VTKWriter<kmds::Mesh> w(mesh_source);
//    w.write("heuristic_3x1x1_source", kmds::R);

    kmds::Mesh mesh_target;
    elg3d::MaterialAssignment ma_target;

    kmds::Connectivity* c_N2R = mesh_source.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch(&mesh_source);
    ch.buildN2R();

    kmds::Connectivity* c_E2R = mesh_source.createConnectivity(kmds::E2R);
    ch.buildEandE2R();

    kmds::Connectivity* c_F2R = mesh_source.createConnectivity(kmds::F2R);
    ch.buildFandF2R_variant_0();


    kmds::Variable<kmds::TCellID> *varSubCells2oldCells =
            mesh_target.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_REGION, "subCells2oldCells");
    kmds::Variable<kmds::TCellID> *varOldCells2firstSubCells =
            mesh_source.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_REGION, "oldCells2firstSubCells");

    kmds::Variable<bool> *varMixedCells_source =
            mesh_source.createVariable<bool>(false, kmds::KMDS_REGION, "mixedCells_source");
    kmds::Variable<bool> *varMixedCells_target =
            mesh_target.createVariable<bool>(false, kmds::KMDS_REGION, "mixedCells_target");

    kmds::Variable<double> *varSurfvol_source =
            mesh_source.createVariable<double>(false, kmds::KMDS_REGION, "surfvol_source");
    kmds::Variable<double> *varSurfvol_target =
            mesh_target.createVariable<double>(false, kmds::KMDS_REGION, "surfvol_target");


    elg3d::InterfaceNodesPosSmoothVF_buildSubMesh_3D(&mesh_source,
                                                     &mesh_target,
                                                     c_N2R,
                                                     c_E2R,
                                                     c_F2R,
                                                     &fp_source,
                                                     &ma_target,
                                                     varSubCells2oldCells,
                                                     varOldCells2firstSubCells,
                                                     varMixedCells_source,
                                                     varMixedCells_target,
                                                     varSurfvol_source,
                                                     varSurfvol_target
    );

//    elg3d::Tools_write_lite_2D(&mesh_target, &ma_target, "heuristic_5x5_before");

    kmds::Connectivity* c_N2R_target = mesh_target.createConnectivity(kmds::N2R);
    kmds::ConnectivityHelper ch_target(&mesh_target);
    ch_target.buildN2R();
    kmds::Connectivity* c_R2R_byN_target = mesh_target.createConnectivity(kmds::R2R_byN);
    ch_target.buildR2R_byN(c_N2R_target);

    elg3d::FracPres fp_target;


    // build graph
    kmds::Variable<kmds::TCellID> *varCells2vertices =
            mesh_target.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_REGION, "cell2vertices");

    kmds::GrowingView<kmds::TCellID> cellsIDs_target("CELLS", mesh_target.getNbRegions());
    mesh_target.getRegionIDs(&cellsIDs_target);

    kmds::GrowingView<kmds::TCellID> selectionCells("CELLS", mesh_target.getNbRegions());

    elg3d::InterfaceNodesPosSmoothVF_selectCells_xD(&cellsIDs_target,
                                                    c_R2R_byN_target,
                                                    varMixedCells_target,
                                                    varCells2vertices,
                                                    &selectionCells);

    kmds::Graph graph("CELLS_SUBMESH_GRAPH", selectionCells.getNbElems(), 20);
    elg3d::InterfaceNodesPosSmoothVF_buildGraph_xD(&selectionCells,
                                                   c_R2R_byN_target,
                                                   varCells2vertices,
                                                   &graph);


    kmds::GrowingView<kmds::TCellID> cellsIDs_source("CELLS", mesh_source.getNbRegions());
    mesh_source.getRegionIDs(&cellsIDs_source);


    const int nbSubPixels = elg3d::InterfaceNodesPosSmoothVF_NBSUB3;


    elg3d::InterfaceNodesPosSmoothVF_assignHeuristic_xD(&cellsIDs_source,
                                                        &cellsIDs_target,
                                                        &fp_source,
                                                        &ma_target,
                                                        varSubCells2oldCells,
                                                        varOldCells2firstSubCells,
                                                        varMixedCells_source,
                                                        varCells2vertices,
                                                        varSurfvol_source,
                                                        varSurfvol_target,
                                                        &graph,
                                                        nbSubPixels
    );


    elg3d::Tools_write_3D(&mesh_target, &fp_target, &ma_target, "heuristic_3x1x1_after");
}
/*----------------------------------------------------------------------------*/
TEST_F(InterfaceNodesPosSmoothVFTest, lp_5x5_2D) {

    kmds::Mesh mesh_source;
    elg3d::FracPres fp_source;

    elg3d::initData_internalFormat(&mesh_source, &fp_source, "Elg3D/test/Samples/internalFormat_5x5_2D.txt");


    kmds::VTKWriter<kmds::Mesh> w(mesh_source);
    w.write("lp_5x5_source", kmds::F);

    kmds::Mesh mesh_target;
    elg3d::MaterialAssignment ma_target;

    kmds::Connectivity* c_N2F = mesh_source.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh_source);
    ch.buildN2F();

    kmds::Connectivity* c_E2F = mesh_source.createConnectivity(kmds::E2F);
    ch.buildEandE2F_2D_variant_1();


    kmds::Variable<kmds::TCellID> *varSubCells2oldCells =
            mesh_target.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "subCells2oldCells");
    kmds::Variable<kmds::TCellID> *varOldCells2firstSubCells =
            mesh_source.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "oldCells2firstSubCells");

    kmds::Variable<bool> *varMixedCells_source =
            mesh_source.createVariable<bool>(false, kmds::KMDS_FACE, "mixedCells_source");
    kmds::Variable<bool> *varMixedCells_target =
            mesh_target.createVariable<bool>(false, kmds::KMDS_FACE, "mixedCells_target");

    kmds::Variable<double> *varSurfvol_source =
            mesh_source.createVariable<double>(false, kmds::KMDS_FACE, "surfvol_source");
    kmds::Variable<double> *varSurfvol_target =
            mesh_target.createVariable<double>(false, kmds::KMDS_FACE, "surfvol_target");


    elg3d::InterfaceNodesPosSmoothVF_buildSubMesh_2D(&mesh_source,
                                                     &mesh_target,
                                                     c_N2F,
                                                     c_E2F,
                                                     &fp_source,
                                                     &ma_target,
                                                     varSubCells2oldCells,
                                                     varOldCells2firstSubCells,
                                                     varMixedCells_source,
                                                     varMixedCells_target,
                                                     varSurfvol_source,
                                                     varSurfvol_target
    );

    elg3d::Tools_write_lite_2D(&mesh_target, &ma_target, "lp_5x5_before");

    kmds::Connectivity* c_N2F_target = mesh_target.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch_target(&mesh_target);
    ch_target.buildN2F();
    kmds::Connectivity* c_F2F_byN_target = mesh_target.createConnectivity(kmds::F2F_byN);
    ch_target.buildF2F_byN(c_N2F_target);

    elg3d::FracPres fp_target;


    // build graph
    kmds::Variable<kmds::TCellID> *varCells2vertices =
            mesh_target.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "cell2vertices");

    kmds::GrowingView<kmds::TCellID> cellsIDs_target("CELLS", mesh_target.getNbFaces());
    mesh_target.getFaceIDs(&cellsIDs_target);

    kmds::GrowingView<kmds::TCellID> selectionCells("CELLS", mesh_target.getNbFaces());

    elg3d::InterfaceNodesPosSmoothVF_selectCells_xD(&cellsIDs_target,
                                                    c_F2F_byN_target,
                                                    varMixedCells_target,
                                                    varCells2vertices,
                                                    &selectionCells);

    kmds::Graph graph("CELLS_SUBMESH_GRAPH", selectionCells.getNbElems(), 20);
    elg3d::InterfaceNodesPosSmoothVF_buildGraph_xD(&selectionCells,
                                                   c_F2F_byN_target,
                                                   varCells2vertices,
                                                   &graph);


    kmds::GrowingView<kmds::TCellID> cellsIDs_source("CELLS", mesh_source.getNbFaces());
    mesh_source.getFaceIDs(&cellsIDs_source);


    const int nbSubPixels = elg3d::InterfaceNodesPosSmoothVF_NBSUB2;


    elg3d::InterfaceNodesPosSmoothVF_assignLP_xD(&cellsIDs_source,
                                                 &cellsIDs_target,
                                                 &fp_source,
                                                 &fp_target,
                                                 &ma_target,
                                                 varSubCells2oldCells,
                                                 varOldCells2firstSubCells,
                                                 varMixedCells_source,
                                                 varCells2vertices,
                                                 &graph,
                                                 nbSubPixels
    );


    elg3d::Tools_write_2D(&mesh_target, &fp_target, &ma_target, "lp_5x5_after");
}
/*----------------------------------------------------------------------------*/
TEST_F(InterfaceNodesPosSmoothVFTest, graphcut_5x5_2D) {

    kmds::Mesh mesh_source;
    elg3d::FracPres fp_source;

    elg3d::initData_internalFormat(&mesh_source, &fp_source, "Elg3D/test/Samples/internalFormat_5x5_2D.txt");


    kmds::VTKWriter<kmds::Mesh> w(mesh_source);
    w.write("graphcut_5x5_source", kmds::F);

    kmds::Mesh mesh_target;
    elg3d::MaterialAssignment ma_target;

    kmds::Connectivity* c_N2F = mesh_source.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch(&mesh_source);
    ch.buildN2F();

    kmds::Connectivity* c_E2F = mesh_source.createConnectivity(kmds::E2F);
    ch.buildEandE2F_2D_variant_1();


    kmds::Variable<kmds::TCellID> *varSubCells2oldCells =
            mesh_target.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "subCells2oldCells");
    kmds::Variable<kmds::TCellID> *varOldCells2firstSubCells =
            mesh_source.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "oldCells2firstSubCells");

    kmds::Variable<bool> *varMixedCells_source =
            mesh_source.createVariable<bool>(false, kmds::KMDS_FACE, "mixedCells_source");
    kmds::Variable<bool> *varMixedCells_target =
            mesh_target.createVariable<bool>(false, kmds::KMDS_FACE, "mixedCells_target");

    kmds::Variable<double> *varSurfvol_source =
            mesh_source.createVariable<double>(false, kmds::KMDS_FACE, "surfvol_source");
    kmds::Variable<double> *varSurfvol_target =
            mesh_target.createVariable<double>(false, kmds::KMDS_FACE, "surfvol_target");


    elg3d::InterfaceNodesPosSmoothVF_buildSubMesh_2D(&mesh_source,
                                                     &mesh_target,
                                                     c_N2F,
                                                     c_E2F,
                                                     &fp_source,
                                                     &ma_target,
                                                     varSubCells2oldCells,
                                                     varOldCells2firstSubCells,
                                                     varMixedCells_source,
                                                     varMixedCells_target,
                                                     varSurfvol_source,
                                                     varSurfvol_target
    );

    elg3d::Tools_write_lite_2D(&mesh_target, &ma_target, "graphcut_5x5_before");

    kmds::Connectivity* c_N2F_target = mesh_target.createConnectivity(kmds::N2F);
    kmds::ConnectivityHelper ch_target(&mesh_target);
    ch_target.buildN2F();
    kmds::Connectivity* c_F2F_byN_target = mesh_target.createConnectivity(kmds::F2F_byN);
    ch_target.buildF2F_byN(c_N2F_target);

    elg3d::FracPres fp_target;


    // build graph
    kmds::Variable<kmds::TCellID> *varCells2vertices =
            mesh_target.createVariable<kmds::TCellID>(kmds::NullID, kmds::KMDS_FACE, "cell2vertices");

    kmds::GrowingView<kmds::TCellID> cellsIDs_target("CELLS", mesh_target.getNbFaces());
    mesh_target.getFaceIDs(&cellsIDs_target);

    kmds::GrowingView<kmds::TCellID> selectionCells("CELLS", mesh_target.getNbFaces());

    elg3d::InterfaceNodesPosSmoothVF_selectCells_xD(&cellsIDs_target,
                                                    c_F2F_byN_target,
                                                    varMixedCells_target,
                                                    varCells2vertices,
                                                    &selectionCells);

    kmds::Graph graph("CELLS_SUBMESH_GRAPH", selectionCells.getNbElems(), 20);
    elg3d::InterfaceNodesPosSmoothVF_buildGraph_xD(&selectionCells,
                                                   c_F2F_byN_target,
                                                   varCells2vertices,
                                                   &graph);


    kmds::GrowingView<kmds::TCellID> cellsIDs_source("CELLS", mesh_source.getNbFaces());
    mesh_source.getFaceIDs(&cellsIDs_source);


    const int nbSubPixels = elg3d::InterfaceNodesPosSmoothVF_NBSUB2;


    kmds::Variable<gmds::math::Point> *varMidPoints_target =
            mesh_target.createVariable<gmds::math::Point>(kmds::NullID, kmds::KMDS_FACE, "midpoints_target");

    elg3d::Tools_computeMidPointVar_2D(&cellsIDs_target,
                                       &mesh_target,
                                       varMidPoints_target);

    elg3d::InterfaceNodesPosSmoothVF_assignGraphCut_xD(&cellsIDs_source,
                                                       &cellsIDs_target,
                                                       &fp_source,
                                                       &ma_target,
                                                       varSubCells2oldCells,
                                                       varOldCells2firstSubCells,
                                                       varMixedCells_source,
                                                       varCells2vertices,
                                                       &graph,
                                                       varMidPoints_target,
                                                       nbSubPixels
    );

    mesh_target.deleteVariable(kmds::KMDS_FACE, "midpoints_target");


    elg3d::Tools_write_2D(&mesh_target, &fp_target, &ma_target, "graphcut_5x5_after");
}
/*----------------------------------------------------------------------------*/