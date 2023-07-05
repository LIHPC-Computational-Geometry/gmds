/*-----------------------------------------------------------------*/
/*
 * Chart.h
 *
 *  Created on: 03/19/2015
 *      Author: ledouxf
 */
/*-----------------------------------------------------------------*/
#ifndef GMDS_MATH_CHART_H_
#define GMDS_MATH_CHART_H_
/*-----------------------------------------------------------------*/
#include <vector>
/*-----------------------------------------------------------------*/
#include <gmds/math/Vector.h>
#include <gmds/math/Matrix.h>
#include "GMDSMath_export.h"
/*-----------------------------------------------------------------*/
#include <iostream>
/*-----------------------------------------------------------------*/
namespace gmds{
    /*-----------------------------------------------------------------*/
    namespace math{
        class Quaternion;
        /*------------------------------------------------------------*/
        /* \class Chart
         * \brief A char is a right-handed basis in R3
         */
        class GMDSMath_API Chart
        {
        public:
            /*------------------------------------------------------------*/
            /* \brief Defaut constructor
             *  Defines a chart with (1,0,0), (0,1,0), (0,0,1)
             */
            Chart();


            /*------------------------------------------------------------*/
            /* \brief Constructor. Warning there is no check about the parameter
             * values
             */
            Chart(const Vector3d& AX, const Vector3d& AY,const Vector3d& AZ);
            /*------------------------------------------------------------*/
            /* \brief Constructor. Warning there is no check about the parameter
             * values. Last vector is computed on-the-fly as AX x AY
             */
            Chart(const Vector3d& AX, const Vector3d& AY);

            /*------------------------------------------------------------*/
            /*  \brief A constructor of chart from quaternion.
             * This constructor is used for conversion purposes.
             */
            explicit Chart(const Quaternion& AQ);


            Vector3d X() const { return m_v[0];}
            Vector3d Y() const { return m_v[1];}
            Vector3d Z() const { return m_v[2];}
            Vector3d get(int AIndex) const;
            Vector3d operator[](int AIndex) const;
            /*------------------------------------------------------------*/
            /*  \brief Return the six vector corresponding to the chirral
             *  rotation group of this
             */
            std::vector<Vector3d> image()const;

            /*------------------------------------------------------------*/
            /*  \brief Return the 3x3 matrix whose each column is one chart
             *         vector
             *
             * \return a matrix view of *this
             */
            Matrix<3,3,double> toMatrix() const;

            /*------------------------------------------------------------*/
            /*  \brief Return the 3x3 rotation matrix that transforms
             *         *this to chart \p AC
             *
             * \param[in] AC a chart we want to value the rotation matrix to
             *
             * \return a rotation matrix bringing $this to \p AC
             */
            Matrix<3,3,double> computeRotationTo(const Chart& AC) const;


            /*------------------------------------------------------------*/
            /*  \brief align the chart with \p AV following the minimal rotation
             *  \param[in] AV a vector to align (*this) along
             */
            void align(const Vector3d& AV);

            /*------------------------------------------------------------*/
            /*  \brief Give the matching of the vectors between this and
             *         AChart. AMatching[i]=j means the i^th vector of this
             *         matchs the j^th vector of AChart. By matching, we
             *         mean the most aligned with (max absolute dot product)
             *
             * \param[in] AChart a chart to compare with
             * \param[in] AMatching a 3-size tabular storing the computed matching
             */
            void matchVectors(const Chart&AChart, int (&AMatching)[3]) const;

            static int testSingularity(const Chart& AC1, const Chart& AC2, const Chart& AC3,
                                       const Chart& AC4);
            static int testSingularity(const Chart& AC1, const Chart& AC2, const Chart& AC3);
            /*------------------------------------------------------------*/
            /* \brief Nested class storing the mapping from one chart to
             *        another one in a compact manner.
             *
             * \details Considering two charts C1 and C2, a mapping object
             *          indicates which vectors of C2 corresponds to each
             *          vector of C1
             */
            class GMDSMath_API Mapping {
            public:
                /*-------------------------------------------------------*/
                /** \brief Default constructor providing the identity
                 *         mapping
                 */
                Mapping();

                /*-------------------------------------------------------*/
                /** \brief Constructs a new mapping from \p AC1 to \p AC2
                 *
                 * \param[in] AC1 an origin chart
                 * \param[in] AC2 a destination chart
                 */
                Mapping(const Chart& AC1, const Chart& AC2);

                /*--------------------------------------------------------*/
                /** \brief Gets the direction vector
                 *
                 * \return a integer vector
                 */

                const Vector3i& getDirections() const ;
                /*--------------------------------------------------------*/
                /** \brief Gets the index permutation vector
                 *
                 * \return a integer vector
                 */
                const Vector3i& getPermutations() const;

                /*-------------------------------------------------------*/
                /** \brief Provide the inverse mapping
                 *
                 * \return the inverse mapping of *this
                 */
                Mapping inverse() const;
                /*-------------------------------------------------------*/
                /** \brief Check if we have the identity mapping
                 */
                bool isIdentity() const;

                /*-------------------------------------------------------*/
                /** \brief Computes the composition of two mappings, the
                 *         resulting mapping consists in applyng \p AM
                 *         then (*this)
                 *
                 * \param[in] AM another mapping
                 * \return the mapping (*this)*\p AM
                 */
                Mapping operator*(const Mapping& AM) const;
                /*-------------------------------------------------------*/
                /** \brief Transport a solution vector by *this. This
                 *         operation must be understand in the context of
                 *         a chart changement. If a vector \p AV is
                 *         expressed in a chart c1 and (*this) is a mapping
                 *         from c1 to c2 then (*this)*\p AV is the
                 *         expression of \p AV in the basis defined by c2.
                 *
                 * \param[in] AV a vector
                 *
                 * \return the vector \p AV expressed in *this
                 */
                Vector3d operator*(const Vector3d& AV) const;

                /*------------------------------------------------------------*/
                /*  \brief Return the 3x3 rotation matrix corresponding to
                 *         this mapping
                 *
                 * \return a rotation matrix made of 0 and (-)1
                 */
                Matrix<3,3,double> toMatrix() const;
            public:
                Vector3i m_map;
                Vector3i m_dir;
            } ;

        private:
            /** three orthogonal vectors */
            Vector3d m_v[3];
            
        };
        /*-----------------------------------------------------------------*/
        GMDSMath_API  std::ostream & operator<<(std::ostream & op_g,
                                               const Chart & op_d);
        /*-----------------------------------------------------------------*/
    }
}
/*-----------------------------------------------------------------*/
#endif /* GMDS_MATH_CHART_H_ */
/*-----------------------------------------------------------------*/
