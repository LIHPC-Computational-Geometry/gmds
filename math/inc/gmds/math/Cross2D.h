/*----------------------------------------------------------------------------*/
/*
 * Cross2D.h
 *
 *  Created on: 03/19/2015
 *      Author: ledouxf
 */
/*-----------------------------------------------------------------*/
#ifndef GMDS_MATH_CROSS_2D_H_
#define GMDS_MATH_CROSS_2D_H_
/*-----------------------------------------------------------------*/
#include "gmds/math/Quaternion.h"
#include "gmds/math/Chart.h"
#include <cmath>
/*-----------------------------------------------------------------*/
#include <iostream>
/*-----------------------------------------------------------------*/
using namespace std;
/*-----------------------------------------------------------------*/
namespace gmds {
    /*-------------------------------------------------------------*/
    namespace math {
        /*---------------------------------------------------------*/
        class EXPORT_GMDS Cross2D
        {

        public:

            static int index(const Cross2D& AC1,
                             const Cross2D& AC2,
                             const Cross2D& AC3);
            /*-----------------------------------------------------*/
            /* Constructor
             *  Defines a cross aligned with (1,0,0)
             */
            Cross2D();
            /*-----------------------------------------------------*/
            /* /brief Constructor
             *  Defines a cross from two orthogonal vectors
             */
            Cross2D(Vector3d& AV1, Vector3d& AV2);
            /*-----------------------------------------------------*/
            /* Constructor
             * Defines from a reference vector
             */
            Cross2D(const Vector3d& ARefVec);
            /*-----------------------------------------------------*/
            /* Constructor
             * Defines from a reference angle
             */
            Cross2D(const TCoord& ARefAngle);

            /*-----------------------------------------------------*/
            /* Copy  Constructor
             */
            Cross2D(const Cross2D& AC);

            /*-----------------------------------------------------*/
            /* Return the vector of this, which is the closest of AN
             */
            Vector3d closestComponentVector(const Vector3d& AN) const;
		  
		   /*-----------------------------------------------------*/
            /* Return the vector of this, which is the closest of AN and also its order (is it the 0,1st, 2nd or 3rd compVector?)
             */
            unsigned int closestComponentVectorAsIndex(const Vector3d& AN) const;

            /*-----------------------------------------------------*/
            /* Return the vector representation of this with 4 orthogonal
             * vectors
             */
            std::vector<Vector3d> componentVectors() const;

            /*-----------------------------------------------------*/
            /* Compute and stores the cross vectors
             */
            void computeComponentVectors();

            /*-----------------------------------------------------*/
            /* Indicates if cross vectors are stored or not
             */
            bool hasStoredComponentVectors() const;

            /*-----------------------------------------------------*/
            /* Return the reference vector
             */
            Vector3d referenceVector() const;

            /*-----------------------------------------------------*/
            /* Return the reference angle
             */
            TCoord referenceAngle() const {return m_angle;};
            /*-----------------------------------------------------*/
            /* Return the reference angle
             */
            TCoord orientedReferenceAngle() const
            {
                if(m_angle>Constants::PI)
                    return Constants::PI2-m_angle;
                else
                    return m_angle;
            };

            /*-----------------------------------------------------*/
            /* Return the angle between this and AC
             */
            TCoord angle(const Cross2D& AC) const;

            /*-----------------------------------------------------*/
            /** \brief  Overloaded operator+ to create a new point from 2 points
             */
            friend EXPORT_GMDS Cross2D operator+(const Cross2D&, const Cross2D&);

            /*-----------------------------------------------------*/
            /** \brief A pondered mean
             * This methods implement a k-pass direct mean pondered by the weights
             *
             * \param ACrosses the crosses we want to compute the mean
             * \param Aweights weights associated to each cross in the mean computation
             * \param ANbSteps the max number of step to converge (default = 5)
             */
            static EXPORT_GMDS Cross2D mean(const vector<Cross2D> & ACrosses,
                                            const vector<TCoord> & AWeights,
                                            const TInt ANbSteps=5);
            static EXPORT_GMDS Cross2D mean(const Cross2D & AC1,
                                            const TCoord& AW1,
                                            const Cross2D & AC2,
                                            const TCoord& AW2);
		  static EXPORT_GMDS Cross2D meanNotMedian(const vector<Cross2D> & ACrosses,
                                            const vector<TCoord> & AWeights,
                                            const TInt ANbSteps=5);
        private:
            // A 2D cross is only represented by an angle, which is the
            // angle between its representation vector and the (OX) axis
            TCoord m_angle;

            // the four vector that defines the 2D cross can be kept in
            // mind or computed on the fly. By default they are not
            // stored
            bool m_known_vectors;
            // storage container of cross vectors
            std::vector<Vector3d> m_vectors;
        };


        EXPORT_GMDS ostream & operator << (ostream & AStr,
                                           const Cross2D & AC);
    }//namespace math
    /*-----------------------------------------------------------------*/
} //namespace gmds
/*-----------------------------------------------------------------*/
#endif /* GMDS_MATH_CROSS_2D_H_ */
/*-----------------------------------------------------------------*/
