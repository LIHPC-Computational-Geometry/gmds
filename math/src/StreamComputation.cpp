/*----------------------------------------------------------------------------*/
// Created by fledoux on 2019-02-12.
/*----------------------------------------------------------------------------*/
#include <gmds/math/StreamComputation.h>
#include <gmds/math/Segment.h>
#include <gmds/math/Plane.h>
/*----------------------------------------------------------------------------*/
using namespace gmds::math;
/*----------------------------------------------------------------------------*/
const int local_id[3][3]= {{0,1,2}, {1,2,0}, {2,0,1} };
/*----------------------------------------------------------------------------*/
bool StreamComputation::RK4(const Point AP[3], const Vector3d AV[3],
                            const Point& APIn, const Vector3d &AVIn,
                            const int& AInCellDim,
                            const int& AInCellId,
                            Point& APOut, Vector3d& AVOut,
                            int &AOutCellDim,
                            int &AOutCellId)
{

    if(AInCellId<0 && AInCellId>2){
        std::string mess = "StreamComputation::RK4) Wrong local input "+std::to_string(AInCellId);
        throw GMDSException(mess);
    }
    if(AInCellDim==0){
        return RK4FromNode(AP,AV,  APIn,AVIn,AInCellId,
                    APOut,AVOut,AOutCellDim,AOutCellId);
    }
    else{
        return RK4FromEdge(AP,AV,APIn,AVIn,AInCellId,
                    APOut,AVOut,AOutCellDim,AOutCellId);
    }
}
/*----------------------------------------------------------------------------*/
bool StreamComputation::RK4FromNode(const Point AP[3], const Vector3d AV[3],
                                    const Point &APIn, const Vector3d &AVIn,
                                    const int &AInCellId,
                                    Point &APOut, Vector3d &AVOut,
                                    int &AOutCellDim, int &AOutCellId)
{
    //We start from a node and so we only try and intersect the opposite edge
    Point p1 = AP[local_id[AInCellId][1]];
    Point p2 = AP[local_id[AInCellId][2]];

    Vector3d v1 = AV[local_id[AInCellId][1]];
    Vector3d v2 = AV[local_id[AInCellId][2]];

    Segment seg(p1,p2);
    //=============================================================
    // 1ST LAUNCH
    //=============================================================
    Plane pl(APIn,AVIn);
    Point pi; //first intersection point
    double w0=0, w1=0;
    Plane::IntersectionType  i_type =pl.intersect(seg, pi, w0, w1);
    if(i_type==Plane::NO_INTERSECTION)
        return false;

    math::Vector3d v_result_1 = w0*v1 + w1*v2;

    //=============================================================
    // 2ND LAUNCH
    //=============================================================
    //We go back to the origin with a new vector
    math::Vector3d v_launched_2 = 0.5*(AVIn+v_result_1);

    Plane pl2(APIn,v_launched_2);
    i_type =pl.intersect(seg, pi, w0, w1);
    if(i_type==Plane::NO_INTERSECTION)
        return false;

    math::Vector3d v_result_2 = w0*v1 + w1*v2;

    //Output now
    AVOut =v_result_2;
    APOut=pi;
    if(i_type==Plane::SEGMENT_FIRST_PNT){
        AOutCellDim=0;
    }
    else if(i_type==Plane::SEGMENT_SECOND_PNT){
        AOutCellDim=0;

    }
    else{ //intersect in the edge so
        AOutCellDim=1;

    }
    return true;
}
/*----------------------------------------------------------------------------*/
bool StreamComputation::RK4FromEdge(const gmds::math::Point *AP, const gmds::math::Vector3d *AV,
                                    const gmds::math::Point &APIn, const gmds::math::Vector3d &AVIn,
                                    const int &AInCellId, gmds::math::Point &APOut, gmds::math::Vector3d &AVOut,
                                    int &AOutCellDim, int &AOutCellId)
{
    std::string mess("StreamComputation::RK4FromEdge) not implemented.");
    throw GMDSException(mess);
}
/*----------------------------------------------------------------------------*/
