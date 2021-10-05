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
