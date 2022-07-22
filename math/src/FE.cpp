// GMDS header files
#include <gmds/math/Constants.h>
#include <gmds/math/Segment.h>
#include <gmds/math/FE.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*------------------------------------------------------------------------*/
    namespace math{
        /*--------------------------------------------------------------------*/
        Vector2d TriangleP1::grad_ref[3] = {
                Vector2d({-1.0,-1.0}),
                Vector2d({ 1.0, 0.0}),
                Vector2d({ 0.0, 1.0})
        };
        /*--------------------------------------------------------------------*/
        Vector3d TetrahedronP1::grad_ref[4] = {
                Vector3d({-1.0,-1.0,-1.0}),
                Vector3d({ 1.0, 0.0, 0.0}),
                Vector3d({ 0.0, 1.0, 0.0}),
                Vector3d({0.0, 0.0, 1.0})
        };
        /*--------------------------------------------------------------------*/
        Matrix<2,2,double> TriangleP1::B(const Point& AP1,
                                         const Point& AP2,
                                         const Point& AP3)
        {
            math::Matrix<2,2,double> s;
            double x1 = AP1.X(), y1 = AP1.Y();
            double x2 = AP2.X(), y2 = AP2.Y();
            double x3 = AP3.X(), y3 = AP3.Y();

            //2 x area(Triangle P1,P2,P3)
            double area2 = x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 -x3*y2;

            double inv_2area =1.0/(area2);
            // First line
            s[0][0]= (x2-x1)*inv_2area;
            s.set(0, 1, (x3-x1)*inv_2area);
            // Second line
            s.set(1, 0, (y2-y1)*inv_2area);
            s.set(1, 1, (y3-y1)*inv_2area);

            return s;
        }
        /*----------------------------------------------------------------------------*/
        Matrix<3,3,double> TriangleP1::stiffnessMatrix(const Point& AP1,
                                                       const Point& AP2,
                                                       const Point& AP3)
        {
            math::Matrix<3,3,double> s;

            Matrix<2,2,double> b = B(AP1,AP2,AP3);

            //We build each gradient using the reference element
            Vector2d grad_real[3] ={
                    b*grad_ref[0],
                    b*grad_ref[1],
                    b*grad_ref[2]
            };
            for(auto i=0;i<4;i++){
                for(auto j=0;j<4;i++){
                    s.set(i, j, grad_real[i].dot(grad_real[j]));

                }
            }
            return s;
        }
        /*----------------------------------------------------------------------------*/
        Matrix<4,4,double> TetrahedronP1::stiffnessMatrix(const Point& AP1,
                                                          const Point& AP2,
                                                          const Point& AP3,
                                                          const Point& AP4)
        {
            math::Matrix<4,4,double> s;
            double x1 = AP1.X(), y1 = AP1.Y();
            double x2 = AP2.X(), y2 = AP2.Y();
            double x3 = AP3.X(), y3 = AP3.Y();

            double x23 = x2-x3;
            double x31 = x3-x1;
            double x12 = x1-x2;
            double y23 = y2-y3;
            double y31 = y3-y1;
            double y12 = y1-y2;

            Vector3d v12=AP2-AP1;
            Vector3d v13=AP3-AP1;
            double area2 = x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 -x3*y2;
            //      double area2 = v12.cross(v13).norm(); //get Area Triangle(AP1,AP2,AP3)
            double area4 = 4*area2*area2;
            double inv_area4 = 1.0/area4;
            // First row
            s.set(0, 0, (y23*y23 + x23*x23)*inv_area4);
            s.set(0, 1, (y23*y31 + x23*x31)*inv_area4);
            s.set(0, 2, (y23*y12 + x23*x12)*inv_area4);
            // Second row
            s.set(1, 0, (y23*y31 + x23*x31)*inv_area4);
            s.set(1, 1, (y31*y31 + x31*x31)*inv_area4);
            s.set(1, 2, (y31*y12 + x31*x12)*inv_area4);
            // Second row
            s.set(2, 0, (y23*y12 + x23*x12)*inv_area4);
            s.set(2, 1, (y31*y12 + x31*x12)*inv_area4);
            s.set(2, 2, (y12*y12 + x12*x12)*inv_area4);

            return s;
        }

        /*----------------------------------------------------------------------------*/
    } // namespace math
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/

