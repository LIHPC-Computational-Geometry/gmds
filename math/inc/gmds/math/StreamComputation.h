/*----------------------------------------------------------------------------*/
// Created by fledoux on 2019-02-12.
/*----------------------------------------------------------------------------*/
#ifndef GMDS_STREAM_COMPUTATION_H
#define GMDS_STREAM_COMPUTATION_H
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
    namespace math {
/*----------------------------------------------------------------------------*/
        class StreamComputation {
        public:
            /**
             *  \brief Method to compute the output point @APout and the output
             *  vector @AVOut. This is an implementation of the RK4 algorithm
            *  where the dt value is given by boundary intersection only.
             *  @AOutCellDim gives the dimension of the cell we go out in the
             *  triangle (0 in a point, 1 in an edge). For dimension 0,
             *  @AOutCellId is 0 for @AP0, 1 for @AP1, 2 for @AP3. For dimension
             *  1, @AOutCellId is 0 for edge [@AP0,@AP1], 1 for edge [@AP1,@AP2]
             *  and 2 for edge [@AP2,@AP0].
             *
             *  @AVin and @AVOut are normal vectors
             *
             * @param[in]  AP   Triangle points
             * @param[in]  AV   Vectors at triangle points
             * @param[in]  APIn  Input point
             * @param[in]  AVIn  Vector at the input point @APIn
             * @param[in]  AInCellDim dimension of the cell where @APIn is.
             * @param[in]  AInCellId  local id of the cell where @APin is.
             * @param[out] APOut Output point
             * @param[out] AVOut Ouput vector at @APout
             * @param[out] AOutCellDim dimension of the cell where @APOut is.
             * @param[out] AOutCellId  local id of the cell where @APOut is.
             *
             * @return true if we truly propagate into the triangle, false otherwise
             */
            static bool RK4(const Point AP[3], const Vector3d AV[3],
                            const Point& APIn, const Vector3d &AVIn,
                            const int& AInCellDim,
                            const int& AInCellId,
                            Point& APOut, Vector3d& AVOut,
                            int &AOutCellDim,
                            int &AOutCellId);
        private:

            static bool RK4FromNode(const Point AP[3], const Vector3d AV[3],
                                    const Point& APIn, const Vector3d &AVIn,
                                    const int& AInCellId,
                                    Point& APOut, Vector3d& AVOut,
                                    int &AOutCellDim,
                                    int &AOutCellId);

            static bool RK4FromEdge(const Point AP[3], const Vector3d AV[3],
                                    const Point& APIn, const Vector3d &AVIn,
                                    const int& AInCellId,
                                    Point& APOut, Vector3d& AVOut,
                                    int &AOutCellDim,
                                    int &AOutCellId);
        };

    }
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_STREAMLINE_H
/*----------------------------------------------------------------------------*/
