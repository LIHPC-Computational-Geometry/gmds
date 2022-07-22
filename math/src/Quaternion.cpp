/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
 * Quaternion.cpp
 *
 *  Created on: 03/01/2015
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#include "gmds/math/Quaternion.h"
#include "gmds/math/Chart.h"
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <cmath>
/*----------------------------------------------------------------------------*/
#define a 1.0/sqrt(2.0)
#define b 1.0/2.0
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*------------------------------------------------------------------------*/
    namespace math {
        /*--------------------------------------------------------------------*/
        const Quaternion Quaternion::img[24] =
        {
            Quaternion(1.0, 0.0, 0.0, 0.0),
            Quaternion(a, a, 0.0, 0.0),
            Quaternion(0.0, 1.0, 0.0, 0.0),
            Quaternion(a, -a, 0.0, 0.0),
            Quaternion(0.0, a, a, 0.0),
            Quaternion(b, -b, -b, b),
            Quaternion(a, 0.0, 0.0, a),
            Quaternion(b, b, b, b),
            Quaternion(a, 0.0, -a, 0.0),
            Quaternion(b, b, -b, b),
            Quaternion(0.0, a, 0.0, a),
            Quaternion(b, -b, -b, -b),
            Quaternion(0.0, 0.0, 1.0, 0.0),
            Quaternion(0.0, 0.0, a, -a),
            Quaternion(0.0, 0.0, 0.0, 1.0),
            Quaternion(0.0, 0.0, a, a),
            Quaternion(a, 0.0, 0.0, -a),
            Quaternion(b, b, -b, -b),
            Quaternion(0.0, a, -a, 0.0),
            Quaternion(b, -b, b, -b),
            Quaternion(b, b, b, -b),
            Quaternion(0.0, a, 0.0, -a),
            Quaternion(b, -b, b, b),
            Quaternion(a, 0.0, a, 0.0)
        };
        //------------------------------------------------------//
        //--------------------CONSTRUCTORS----------------------//
        //------------------------------------------------------//
        
        Quaternion::Quaternion(const Chart & t)
        {
            int choice = 0; // default choice
            Vector3d U = t.X();
            Vector3d V = t.Y();
            Vector3d W = t.Z();
            
            //we first determine the greates diagonal value
            if (((U[0] * U[0]) < (V[1] * V[1])) &&
                ((W[2] * W[2]) < (V[1] * V[1])))
            {
                choice = 1;
            }
            
            if ((U[0] * U[0] < W[2] * W[2]) &&
                (W[2] * W[2] > V[1] * V[1]))
            {
                choice = 2;
            }
            double Quu, Quv, Quw, Qvu, Qvv, Qvw, Qwu, Qwv, Qww;
            //Underlying idea : doing the affectation in the test loops
            //without using choice variable;
            switch (choice) {
                case 0:
                    Quu = U[0];
                    Qvu = U[1];
                    Qwu = U[2];
                    Quv = V[0];
                    Qvv = V[1];
                    Qwv = V[2];
                    Quw = W[0];
                    Qvw = W[1];
                    Qww = W[2];
                    break;
                case 1:
                    Quu = V[1];
                    Qvu = V[2];
                    Qwu = V[0];
                    Quv = W[1];
                    Qvv = W[2];
                    Qwv = W[0];
                    Quw = U[1];
                    Qvw = U[2];
                    Qww = U[0];
                    break;
                case 2:
                    Quu = W[2];
                    Qvu = W[0];
                    Qwu = W[1];
                    Quv = U[2];
                    Qvv = U[0];
                    Qwv = U[1];
                    Quw = V[2];
                    Qvw = V[0];
                    Qww = V[1];
                    break;
            }
            
            if ((1 + Quu - Qvv - Qww) < 0.01){
                if (Quu < 0.0){
                    if (Qvv > Qww){
                        Quu = -Quu;
                        Qvu = -Qvu;
                        Qwu = -Qwu;
                        Quv = -Quv;
                        Qvv = -Qvv;
                        Qwv = -Qwv;
                    }
                    else{
                        Quu = -Quu;
                        Qvu = -Qvu;
                        Qwu = -Qwu;
                        Quw = -Quw;
                        Qvw = -Qvw;
                        Qww = -Qww;
                    }
                }
                else{
                    Quv = -Quv;
                    Qvv = -Qvv;
                    Qwv = -Qwv;
                    Quw = -Quw;
                    Qvw = -Qvw;
                    Qww = -Qww;
                }
            }
            
            double r = sqrt(1 + Quu - Qvv - Qww);
            if (!(r > (-0.1))){ //Never done
                r = 0.0;
            }
            
            if (r == 0.0) {
                m_r = 1.0;
                m_i[0] = 0.0;
                m_i[1] = 0.0;
                m_i[2] = 0.0;
            }
            else {
                switch (choice)
                {
                    case 0:
                        m_r = (Qwv - Qvw) / (2.0 * r);
                        m_i[0] = r / 2.0;
                        m_i[1] = (Quv + Qvu) / (2.0 * r);
                        m_i[2] = (Qwu + Quw) / (2.0 * r);
                        break;
                    case 1:
                        m_r = (Qwv - Qvw) / (2.0 * r);
                        m_i[1] = r / 2.0;
                        m_i[2] = (Quv + Qvu) / (2.0 * r);
                        m_i[0] = (Qwu + Quw) / (2.0 * r);
                        break;
                    case 2:
                        m_r = (Qwv - Qvw) / (2.0 * r);
                        m_i[2] = r / 2.0;
                        m_i[0] = (Quv + Qvu) / (2.0 * r);
                        m_i[1] = (Qwu + Quw) / (2.0 * r);
                        break;
                }
            }
            this->normalize();
        }
        
        
        /*----------------------------------------------------------------------------*/
        vector<TCoord> Quaternion::getVal()const
        {
            vector<TCoord> v;
            v.reserve(4);
            v.push_back(m_r);
            v.push_back(m_i[0]);
            v.push_back(m_i[1]);
            v.push_back(m_i[2]);
            return v;
        }
        
        /*----------------------------------------------------------------------------*/
        Vector3d Quaternion::axis()const
        {
            double eps = .00001;
            double ratio = sin(acos(m_r));
            if ( fabs(ratio) < eps )
                return Vector3d({0.0, 0.0, 0.0});
            else
                return  m_i/ratio;
        }
        /*----------------------------------------------------------------------------*/
        
        Quaternion Quaternion::SLERP(const Quaternion& AFrom,
                                     const Quaternion& ATo,
                                     const TCoord AParam){
            double eps = .00001;
            Quaternion to;
            
            // calculate cosine
            double c = (AFrom.imaginaryPart().dot(ATo.imaginaryPart()) +
                        AFrom.realPart() + ATo.realPart());
            
            // Adjust signs (if necessary)
            if ( c < 0.0 ) {
                c = -c;
                to = ATo.opposite();
            } else {
                to = ATo;
            }
            
            double s0, s1;
            
            // Calculate coefficients
            if ((1.0 - c) > eps ) {
                double o = acos( c );
                double s = sin( o );
                s0 = sin((1.0 - AParam) * o) / s;
                s1 = sin(AParam * o) / s;
            } else {
                // 'from' and 'to' are very close - just do linear interpolation
                s0 = 1.0 - AParam;
                s1 = AParam;
            }
            return (Quaternion(s0*AFrom.realPart(),
                               s0*AFrom.imaginaryPart()) +
                    Quaternion(s1*to.realPart(),
                               s1*to.imaginaryPart()));
        }
        /*----------------------------------------------------------------------------*/
        
        Quaternion Quaternion::closestImg(const Quaternion & item)const
        {
            const Quaternion candidate = this->conjugate()*item;
            double max = -1;
            unsigned int indexMax;
            for (unsigned int i = 0; i < 24; i++)
            {
                const double dist = fabs(candidate.dot(Quaternion::img[i]));
                if (dist > max)
                {
                    max = dist;
                    indexMax = i;
                }
            }
            if (candidate.dot(Quaternion::img[indexMax]) >= 0.0)
                return (*this)*Quaternion::img[indexMax];
            else
                return ((*this)*Quaternion::img[indexMax]).opposite();
        }
        /*----------------------------------------------------------------------------*/
        
        void Quaternion::setFromEulerAngle(TCoord& AAngleX,
                                           TCoord& AAngleY,
                                           TCoord& AAngleZ)
        {
            
            TCoord cos_z_2 = cos(0.5*AAngleX);//cos(0.5*angleZ);
            TCoord cos_y_2 = cos(0.5*AAngleY);
            TCoord cos_x_2 = cos(0.5*AAngleZ);//cos(0.5*angleX);
            TCoord sin_z_2 = sin(0.5*AAngleX);//sin(0.5*angleZ);
            TCoord sin_y_2 = sin(0.5*AAngleY);
            TCoord sin_x_2 = sin(0.5*AAngleZ);//sin(0.5*angleX);
            
            
            m_r    = cos_z_2*cos_y_2*cos_x_2 - sin_z_2*sin_y_2*sin_x_2;
            m_i[0] = sin_z_2*cos_y_2*cos_x_2 + cos_z_2*sin_y_2*sin_x_2;
            m_i[1] = cos_z_2*sin_y_2*cos_x_2 - sin_z_2*cos_y_2*sin_x_2;
            m_i[2] = cos_z_2*cos_y_2*sin_x_2 + sin_z_2*sin_y_2*cos_x_2;
            
        }
        /*----------------------------------------------------------------------------*/
        
        void Quaternion::toEulerAngle(TCoord& angleX, TCoord& angleY, TCoord& angleZ) const
        {
            //we first select the right Quaternion
            TCoord min = 1000;
            unsigned int indexMin;
            for (unsigned int i = 0; i < 24; i++)
            {
                Quaternion candidate = (*this)*Quaternion::img[i];
                TCoord testTmp = candidate.X() * candidate.J()
                + candidate.I() * candidate.K();
                if (testTmp * testTmp < min){
                    min = testTmp * testTmp;
                    indexMin = i;
                }
            }
            Quaternion candidate = (*this)*Quaternion::img[indexMin];
            
            //TODO VERIFIER ORDRE DES ROTATIONS
            const TCoord PI = 3.141592;
            //	TCoord lw=x, lx= i, ly=j, lz=k;
            //	TCoord w2 = x*x; //q0
            //	TCoord x2 = i*i; //q1
            //	TCoord y2 = j*j; //q2
            //	TCoord z2 = k*k; //q3
            TCoord lw = candidate.X();
            TCoord lx = candidate.I();
            TCoord ly = candidate.J();
            TCoord lz = candidate.K();
            //	TCoord w2 = candidate.X() * candidate.X(); //q0
            TCoord x2 = candidate.I() * candidate.I(); //q1
            TCoord y2 = candidate.J() * candidate.J(); //q2
            TCoord z2 = candidate.K() * candidate.K(); //q3
            
            TCoord phi, theta, psi;
            
            TCoord test = lw * ly + lz * lx;
            
            //les deux cas suivants correspondent aux poles Nord et Sud. C'est le
            //probleme du "Gimbal Lock"
            if (test > 0.49)
            {
                std::cout << "Warning:  gimbal lock !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
                std::cout << "lw " << lw << " lx " << lx << " ly " << ly << " lz " << lz << std::endl;
                phi = -2 * atan2(lw, lx);
                theta = PI / 2;
                psi = PI;
            }
            else if (test < -0.49)
            {
                std::cout << "Warning: gimbal lock !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
                std::cout << "lw " << lw << " lx " << lx << " ly " << ly << " lz " << lz << std::endl;
                phi = -2 * atan2(lw, lx);
                theta = -PI / 2;
                psi = PI;
            }
            else{
                
                phi = atan2(2 * (lw * lx - ly * lz), 1 - 2 * (x2 + y2));
                theta = asin(2 * (lx * lz + ly * lw));
                psi = atan2(2 * (lw * lz - lx * ly), 1 - 2 * (y2 + z2));
            }
            angleZ = psi;//angleX = psi;
            angleY = theta;
            angleX = phi;//angleZ = phi;
            
            
            if (!(angleZ > -5) && !(angleZ<5)){
                std::cout << "PROBLEM ANGLE Z" << std::endl;
            }
            if (!(angleX>-5) && !(angleX<5)){
                std::cout << "PROBLEM ANGLE X" << std::endl;
            }
            if (!(angleY>-5) && !(angleY < 5)){
                std::cout << "PROBLEM ANGLE Y" << std::endl;
            }
            
        }
        //------------------------------------------------------//
        //---------------------OPERATORS------------------------//
        //------------------------------------------------------//
        
        
        bool Quaternion::operator!=(const Quaternion & AQ)const
        {
            return ((this->m_r != AQ.m_r) ||
                    (this->m_i[0] != AQ.m_i[0]) ||
                    (this->m_i[1] != AQ.m_i[1]) ||
                    (this->m_i[2] != AQ.m_i[2]));
        }
        
        ostream & operator << (ostream & AStream, const Quaternion & AQ)
        {
            AStream << "Quat(" << AQ.X() << ", " << AQ.I()
            << ", " << AQ.J() << ", " << AQ.K()<<")";
            return AStream;
        }
        
        //------------------------------------------------------//
        //-------------------TOOL FUNCTIONS---------------------//
        //------------------------------------------------------//
        Quaternion Quaternion::getAleat()
        {
            double q1 = rand() / (double)RAND_MAX;
            double q2 = rand() / (double)RAND_MAX;
            double q3 = rand() / (double)RAND_MAX;
            double q4 = rand() / (double)RAND_MAX;
            return Quaternion(q1, q2, q3, q4);
        }
        /*----------------------------------------------------------------------------*/
        Quaternion Quaternion::simpleMean(const vector<Quaternion> & AQuats,
                                          const vector<TCoord> &     AWeights,
                                          const Quaternion &         ARef)
        {
            TCoord x = 0.0;
            TCoord i = 0.0;
            TCoord j = 0.0;
            TCoord k = 0.0;
            for (unsigned int ci = 0; ci < AQuats.size(); ci++)
            {
                const Quaternion current_quat = AQuats[ci].closestImg(ARef);
                x += current_quat.X()*AWeights[ci];
                i += current_quat.I()*AWeights[ci];
                j += current_quat.J()*AWeights[ci];
                k += current_quat.K()*AWeights[ci];
            }
            return Quaternion(x, i, j, k);
        }
        /*----------------------------------------------------------------------------*/
        Quaternion Quaternion::mean(const vector<Quaternion> & AQuats,
                                    const vector<TCoord> & AWeights)
        {
            const unsigned int nb_steps = 5;
            
            if (AQuats.empty())
            {
                throw GMDSException("The mean of zero Quaternion is undefined.");
            }
            
            if (AWeights.size() != AQuats.size())
            {
                throw GMDSException("There must exactly the same number of Quaternions and coefficients.");
            }
            
            Quaternion ref = Quaternion::simpleMean(AQuats, AWeights, AQuats.front());
            Quaternion prev_ref;
            unsigned int step_index = 0;
            
            while ((step_index < nb_steps) && (ref != prev_ref)) {
                prev_ref = ref;
                ref = Quaternion::simpleMean(AQuats, AWeights, ref);
                step_index++;
            }
            
            return ref;
        }
        /*----------------------------------------------------------------------------*/
        Quaternion Quaternion::mean(const Quaternion& AQ1, const TCoord& AC1,
                                    const Quaternion& AQ2, const TCoord& AC2)
        {
            const Quaternion q1 = AQ1;
            const Quaternion q2 = AQ2.closestImg(q1);
            TCoord x = AC1*q1.X()+AC2*q2.X();
            TCoord i = AC1*q1.I()+AC2*q2.I();
            TCoord j = AC1*q1.J()+AC2*q2.J();
            TCoord k = AC1*q1.K()+AC2*q2.K();
            
            return Quaternion(x, i, j, k);
        }
        /*----------------------------------------------------------------------------*/
        int Quaternion::testSingularity(Quaternion& q1, Quaternion& q2,
                                        Quaternion& q3, Quaternion& q4)
        {
            int sing_type = 0;
            //We value a chart for each quaternion
            
            math::Chart c[4];
            c[0] = math::Chart(q1);
            c[1] = math::Chart(q2);
            c[2] = math::Chart(q3);
            c[3] = math::Chart(q4);

            //FACE 012
            math::Chart::Mapping m01(c[0],c[1]);
            math::Chart::Mapping m12(c[1],c[2]);
            math::Chart::Mapping m20(c[2],c[0]);

            if((m20*m12*m01).isIdentity()==false){
                sing_type++;
            }
            math::Chart::Mapping m13(c[1],c[3]);
            math::Chart::Mapping m30(c[3],c[0]);
            if((m30*m13*m01).isIdentity()==false){
                sing_type++;
            }

            //FACE 023
            math::Chart::Mapping m02(c[0],c[2]);
            math::Chart::Mapping m23(c[2],c[3]);
            if((m30*m23*m02).isIdentity()==false){
                sing_type++;}

            //FACE 123
            math::Chart::Mapping m31(c[3],c[1]);
            if((m31*m23*m12).isIdentity()==false){
                sing_type++;};


            return sing_type;
            
        }
        /*----------------------------------------------------------------------------*/
        int Quaternion::testSingularity(Quaternion& q1, 
                                        Quaternion& q2,
                                        Quaternion& q3)
        {
            int sing_type = 0;
            //We value a chart for each quaternion
            
            math::Chart charts[3];
            charts[0] = math::Chart(q1);
            charts[1] = math::Chart(q2);
            charts[2] = math::Chart(q3);
            
            
            int matchingTab[3][3][3];
            for (unsigned int i = 0; i < 3; i++){
                for (unsigned int j = 0; j < 3; j++){
                    charts[i].matchVectors(charts[j],matchingTab[i][j]);
                }
            }
            
            for (unsigned int i = 0; i < 3; i++){//We take the vectors of charts[0]
                bool goodMatching = true;
                for (unsigned int j = 1; j < 3; j++){
                    for (unsigned int k = 1; k < 3; k++){
                        if (matchingTab[j][k][matchingTab[0][j][i]] != matchingTab[0][k][i]){
                            goodMatching = false;
                        }
                    }
                }
                if (!goodMatching){
                    sing_type++;
                }
            }
            
            return sing_type;
            
        }
        /*----------------------------------------------------------------------------*/
        Quaternion Quaternion::alignWith(const math::Vector3d& AN)
        {
            /* We align *this with AN by making a rotation of c= math::Chart(*this)
             * to be aligned with AN
             *
             * This rotation is minimal; it takes the closest vector v of c with N
             * and computes the rotation that align v with n while keeping
             * v.cross(N) stable (invariant)
             */
            Chart my_chart(*this);
            
            my_chart.align(AN);
            return Quaternion(my_chart);
            
        }
        /*--------------------------------------------------------------------*/
    } // namespace math
    /*------------------------------------------------------------------------*/
    
} //namespace gmds

#undef a
#undef b
/*----------------------------------------------------------------------------*/
