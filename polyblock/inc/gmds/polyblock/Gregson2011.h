/*----------------------------------------------------------------------------*/
#ifndef POLYCUBE_GREGSON2011_H
#define POLYCUBE_GREGSON2011_H
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <gmds/ig/Mesh.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
/*----------------------------------------------------------------------------*/
using namespace std;
/*----------------------------------------------------------------------------*/
typedef Eigen::Matrix<double, 2, 1> vect2D;
/*----------------------------------------------------------------------------*/
gmds::math::Vector eigen2math(vect2D);
/*----------------------------------------------------------------------------*/
namespace gmds {

    class Gregson2011_2D{
    private:
        Mesh *mesh;
        int mark_smooth_bnd_nodes;
        int mark_bnd_nodes;
        int mark_bnd_egdes;
    public:

        Gregson2011_2D(Mesh *mesh);
        void transformMeshToPolycube();
        Mesh getPolycube();
        Mesh getMesh();

    private:
        void markOnBoundary();
        void updateBoundary();
        vect2D normal(Edge);
        vect2D getClosestCanonical(vect2D);
        int getCanonicalId(vect2D v);
        void projectOnPolycube(Variable<vect2D>* edge_chosenNormal);

    };

}
/*----------------------------------------------------------------------------*/
#endif
/*----------------------------------------------------------------------------*/
