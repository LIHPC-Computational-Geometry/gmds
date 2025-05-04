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
                         +(- wi_0*wj_0*AG[0  ][0  ])
                         +(- wi_0*wj_1*AG[0  ][N-1])
                         +(- wi_1*wj_0*AG[M-1][0  ])
                         +(- wi_1*wj_1*AG[M-1][N-1]);
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
                         +(- wi_0*wj_0*AG(0  ,0  ))
                         +(- wi_0*wj_1*AG(0  ,N-1))
                         +(- wi_1*wj_0*AG(M-1,0  ))
                         +(- wi_1*wj_1*AG(M-1,N-1));
        }
    }
    return true;
}
/*------------------------------------------------------------------------*/
bool TransfiniteInterpolation::computeTri(TriArray<Point> &AGrid) {
    // the transfinite interpolation is given by formulation in
    // https://www.ljll.math.upmc.fr/perronnet/transfini/transfini.html
    // we also check how gmsh does the transfinite interpolation
    //first side = first line
    int dim = AGrid.size();
    int n = dim-1;

    math::Point p100 = AGrid(n,0,0);
    math::Point p010 = AGrid(0,n,0);
    math::Point p001 = AGrid(0,0,n);

    for(auto i=1; i<n-1; i++){
        for(auto j=1; j<n-i; j++){
            //i and j are parametric coordinate in a discrete way
            //corresponding barycentric weights are
            TCoord s = (double)i/(double)(n);
            TCoord t = (double)j/(double)(n);
            TCoord a1 = s;
            TCoord a2 = t;
            TCoord a3 = 1-s-t;
            auto b1 = i;
            auto b2 = j;
            auto b3 = n-b1-b2;
            auto k=n-i-j;
            AGrid(i,j,k)=
                    a1*((AGrid(n-b2,b2,0)+AGrid(n-b3,0,b3))+math::Point(-p100.X(),-p100.Y(),-p100.Z()))+
                    a2*((AGrid(0,n-b3,b3)+AGrid(b1,n-b1,0))+math::Point(-p010.X(),-p010.Y(),-p010.Z()))+
                    a3*((AGrid(b1,0,n-b1)+AGrid(0,b2,n-b2))+(math::Point(-p001.X(),-p001.Y(),-p001.Z())));
        }
    }
    return true;
}

