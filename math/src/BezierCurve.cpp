
/*----------------------------------------------------------------------------*/
/*
 * BezieCurve.cpp
 *
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/BezierCurve.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*--------------------------------------------------------------------------*/
    namespace math{
        /*------------------------------------------------------------------------*/
        BezierCurve::BezierCurve(const Point& AP1, const Point& AP2,
                                 const Point& AP3) {
            m_control_points.resize(3);
            m_control_points[0] = AP1;
            m_control_points[1] = AP2;
            m_control_points[2] = AP3;
        }
        /*------------------------------------------------------------------------*/
        BezierCurve::BezierCurve(const Point& AP1, const Point& AP2,
                                 const Point& AP3, const Point& AP4) {
            m_control_points.resize(4);
            m_control_points[0] = AP1;
            m_control_points[1] = AP2;
            m_control_points[2] = AP3;
            m_control_points[3] = AP4;
        }
        /*------------------------------------------------------------------------*/
        BezierCurve:: BezierCurve(const Point& AP1, const Vector3d& AV1,
                                  const Point& AP2, const Vector3d& AV2){
            m_control_points.resize(4);
            m_control_points[0] = AP1;

            m_control_points[1] = Point(AP1.X()+(1.0/3.0)*AV1.X(),
                                        AP1.Y()+(1.0/3.0)*AV1.Y(),
                                        AP1.Z()+(1.0/3.0)*AV1.Z());
            m_control_points[2] = Point(AP2.X()-(1.0/3.0)*AV2.X(),
                                        AP2.Y()-(1.0/3.0)*AV2.Y(),
                                        AP2.Z()-(1.0/3.0)*AV2.Z());
            m_control_points[3] = AP2;
        }

        /*------------------------------------------------------------------------*/
        BezierCurve::BezierCurve(const std::vector<Point>& APts) {
            m_control_points = APts;
        }
        /*------------------------------------------------------------------------*/
        Point BezierCurve::operator()(const double& AT) const {
            if(AT>1.0 || AT<0.0)
                throw GMDSException("BezierCurve Query out of the range [0,1]");


            std::vector<Point> current_control = m_control_points;
            while(current_control.size()>1){

                std::vector<Point> next_control;
                const int deg = current_control.size();

                next_control.resize(deg-1);
                for(int i=0; i<deg-1; i++){
                    next_control[i] = (1-AT)*current_control[i]+ AT*current_control[i+1];
                }

                current_control = next_control;
            }

            return current_control[0];
        }
        /*------------------------------------------------------------------------*/
        std::vector<Point> BezierCurve:: getDiscretization(const int ANb) const {
            if(ANb<1)
                throw GMDSException("BezierCurve discretization impossible with this parameter");


            std::vector<Point> points;
            points.resize(ANb+1);
            double step = 1.0/ANb;
            double val =0;
            for(int i=0; i<=ANb; i++){
                points[i] = this->operator()(val);
                val+=step;
            }
            return points;
        }

        /*----------------------------------------------------------------------------*/
    } // namespace math
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
