/*------------------------------------------------------------------------*/
#include <gmds/math/TransfiniteInterpolation.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::math;
/*------------------------------------------------------------------------*/
bool TransfiniteInterpolation::
compute(std::vector<std::vector<Point>> &AG) {
    //We first check that the dimensions of AG are correct
    auto M = AG.size();
    if(M==0)
        return false;
    auto N = AG[0].size();
    for(auto t:AG)
        if(t.size()!=N)
            return false;

    for(auto i=1;i<M-1;i++){

        double wi_0 = (double)(M-1-i)/(double)(M-1);
        double wi_1 = (double)i/(double)(M-1);

        for(auto j=1;j<N-1;j++){

            double wj_0 = (double)(N-1-j)/(double)(N-1);
            double wj_1 = (double)j/(double)(N-1);

            AG[i][j] =   wi_0*AG[0][j] + wi_1*AG[M-1][j]
                         + wj_0*AG[i][0] + wj_1*AG[i][N-1]
                         - wi_0*wj_0*AG[0  ][0  ]
                         - wi_0*wj_1*AG[0  ][N-1]
                         - wi_1*wj_0*AG[M-1][0  ]
                         - wi_1*wj_1*AG[M-1][N-1];
        }

    }
    return true;
}
/*------------------------------------------------------------------------*/
bool TransfiniteInterpolation::
computeQuad(Array2D<Point> &AG) {
    //We first check that the dimensions of AG are correct
    auto M = AG.nbLines();
    if(M==0)
        return false;
    auto N = AG.nbColumns();

    for(auto i=1;i<M-1;i++){

        double wi_0 = (double)(M-1-i)/(double)(M-1);
        double wi_1 = (double)i/(double)(M-1);

        for(auto j=1;j<N-1;j++){

            double wj_0 = (double)(N-1-j)/(double)(N-1);
            double wj_1 = (double)j/(double)(N-1);

            AG(i,j) =   wi_0*AG(0,j) + wi_1*AG(M-1,j)
                         + wj_0*AG(i,0) + wj_1*AG(i,N-1)
                         - wi_0*wj_0*AG(0  ,0  )
                         - wi_0*wj_1*AG(0  ,N-1)
                         - wi_1*wj_0*AG(M-1,0  )
                         - wi_1*wj_1*AG(M-1,N-1);
        }
    }
    return true;
}
/*------------------------------------------------------------------------*/
/* s4 +-----c3-----+ s3
      |            |
      |            |
     c4            c2
      |            |
      |            |
   s1 +-----c1-----+ s2
   For the triangular case, s1=s4=c4,
   and we have the transfinite interpolation defined as
        tfi_triangle(c1, c2, c3, s1, s2, s3, u, v) =
        u *c2 + (1. - v) * c1 + v *c3 - (u * (1. - v) * s2 + u * v * s3)
*/
/*------------------------------------------------------------------------*/
bool TransfiniteInterpolation::computeTri(Array2D<Point> &AGrid) {
    // the transfinite interpolation is given by formulation in
    // https://www.ljll.math.upmc.fr/perronnet/transfini/transfini.html
    // we also check how gmsh does the transfinite interpolation
    //first side = first line
    int dim = AGrid.nbLines();
    int dim_2 = AGrid.nbColumns();

    //transfinite interpolation is only available for square data
    if(dim!=dim_2)
        return false;

    math::Point s2 = AGrid(0,0);
    math::Point s1 = AGrid(dim-1,0);
    math::Point s3 = AGrid(0,dim-1);

    for(auto i=1; i<dim-1; i++){
        for(auto j=1; j<dim-1-i; j++){
            TCoord u = (double)j/(double)(dim-1);
            TCoord v = (double)i/(double)(dim-1);
            math::Point c2 = AGrid(0,j);
            math::Point c1 = AGrid(i,0);
            math::Point c3 = AGrid(i,dim-1-i);
            AGrid(i,j)= u *c2 + (1. - v) * c1 + v *c3
                        - (u * (1. - v) * s2 + u * v * s3);
        }
    }
    return true;
}

