/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    main_elg2d.cpp
 *  \author  legoff
 *  \date    07/04/2018
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
#include <GMDS/CAD/FacetedGeomManager.h>
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <ELG3D/ALGOCMPNT/AssignCells.h>
#include <ELG3D/ALGOCMPNT/BoundingBoxGeomAssociation.h>
#include <ELG3D/ALGOCMPNT/InitData.h>
#include <ELG3D/ALGOCMPNT/InterfaceNodesPos.h>
#include <ELG3D/ALGOCMPNT/ManifoldDetection.h>
#include <ELG3D/ALGOCMPNT/MaterialGradientComputation.h>
#include <ELG3D/ALGOCMPNT/MaterialInterfaces.h>
#include <ELG3D/ALGOCMPNT/MoveToNewPos.h>
#include <ELG3D/ALGOCMPNT/Pillow.h>
#include <ELG3D/ALGOCMPNT/SmartLaplacian.h>
#include <ELG3D/ALGOCMPNT/Tools.h>
#include <ELG3D/DATACMPNT/FracPres.h>
#include <ELG3D/DATACMPNT/MaterialAssignment.h>
/*----------------------------------------------------------------------------*/

struct fillgv
{
    kmds::GrowingView<kmds::TCellID>* gv;

    fillgv(
            kmds::GrowingView<kmds::TCellID>* gv_)
            : gv(gv_)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void
    operator()(int i) const
    {
//        if(i%2 == 0) {
            gv->push_back(i);
        //}
    }
};


int main(int argc, char* argv[])
{
    int num_threads = -1;

    std::istringstream iss( argv[1] );
    int val;

    if (iss >> val)
    {
        num_threads = val;
    } else {
        std::cerr<<"could not convert number of threads argument."<<std::endl;
        exit(-1);
    }

    Kokkos::InitArguments kargs;
    kargs.num_threads = num_threads;
    Kokkos::initialize(kargs);



    Kokkos::Timer timer;
    timer.reset();

    const int nbelems = 100000000;
    kmds::GrowingView<kmds::TCellID> gv("GV", nbelems);

    timer.reset();
    Kokkos::parallel_for(nbelems,
                         fillgv(&gv));
    Kokkos::fence();

    std::cout<<"timer gv "<<timer.seconds()<<std::endl;

    std::cout<<"gv.getNbElems "<<gv.getNbElems()<<std::endl;

//    Kokkos::UnorderedMap<kmds::TCellID, void> kmap(nbelems*2);
//    timer.reset();
//    Kokkos::parallel_for(nbelems,
//                         KOKKOS_LAMBDA(const int i) {
//                             kmap.insert(i);}
//                             );
//    Kokkos::fence();
//
//    std::cout<<"timer kmap "<<timer.seconds()<<std::endl;
//
//    std::cout<<"kmap.size() "<<kmap.size()<<std::endl;

    Kokkos::View<kmds::TCellID *> nb_tmp("nb_tmp", 7+1);
    nb_tmp(0) = 2;
    nb_tmp(1) = 2;
    nb_tmp(2) = 3;
    nb_tmp(3) = 7;
    nb_tmp(4) = 0;
    nb_tmp(5) = 1;
    nb_tmp(6) = 2;
    Kokkos::parallel_scan(7+1,
                          KOKKOS_LAMBDA(const int& i, kmds::TCellID& upd, const bool& final)
                          {
                              std::cout<<"aaa "<<i<<" "<<nb_tmp(i)<<std::endl;
                              const kmds::TCellID val_i = nb_tmp(i);
                              if(final) {
                                  nb_tmp(i) = upd;
                                  std::cout<<"bbb "<<i<<" "<<nb_tmp(i)<<" "<<upd<<std::endl;
                              }

                              upd += val_i;
                          }
    );

    for(int i=0; i<7+1; i++) {
        std::cout << "nb_tmp("<<i<<") "<<nb_tmp(i)<<std::endl;
    }

    Kokkos::finalize();
    return 0;
}