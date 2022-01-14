#ifndef METRIC_H
#define METRIC_H
/******************************************************************************/
#include<gmds/utils/CommonTypes.h>
/******************************************************************************/
#include<gmds/hybridMeshAdapt/SimplicesNode.h>
/******************************************************************************/
#include<gmds/math/Point.h>
#include<gmds/math/Vector.h>
#include <gmds/math/Matrix.h>
/******************************************************************************/
#include <Eigen/Eigen>
/******************************************************************************/
#include <math.h>
/******************************************************************************/
namespace gmds
{
    template<typename M>
    class Metric
    {
    public:

      Metric(const M & metric);

      ~Metric();

      bool operator==(const Metric & M1) const ;

      const M& getMetric() const ;

      void setMetric(const M& metric);

      M interpolateMetric(const M& m2, const double t);

      M intersectionMetric(const M& m2);

      double metricDist(const math::Vector3d& coordA, const math::Vector3d& coordB, const Metric& m2) const ;

    private:

      const M& m_metric;
    };



    template<typename M>
    Metric<M>::Metric(const M & metric) : m_metric(metric)
    {
    }

    template<typename M>
    Metric<M>::~Metric()
    {

    }

    template<typename M>
    bool Metric<M>::operator==(const Metric & M1) const
    {
      return (getMetric() == M1.getMetric());
    }

    template<typename M>
    const M& Metric<M>::getMetric() const
    {
      return m_metric;
    }

    template<typename M>
    void Metric<M>::setMetric(const M& metric)
    {
      m_metric = metric;
    }

    template<typename M>
    M Metric<M>::interpolateMetric(const M& m2, const double t)
    {
     //TODO
    }

    template<typename M>
    M Metric<M>::intersectionMetric(const M& m2)
    {
      //TODO
    }

    template<typename M>
    double Metric<M>::metricDist(const math::Vector3d& coordA, const math::Vector3d& coordB, const Metric& m2) const
    {
      double metriclenght = 0.0;
      if(std::is_same<M, Eigen::Matrix3d>::value == true)
      {
        math::Vector3d vecAB = coordB - coordA;
        Eigen::Vector3d evecAB = Eigen::Vector3d(vecAB.X(), vecAB.Y(), vecAB.Z());
        if(*this == m2)
        {
          metriclenght = sqrt(evecAB.transpose() * m_metric * evecAB);
        }
        else
        {
          //unsigned int samplemax = 10;
          //use a loop to have a better approximation of a metric lenght edge
          //https://pages.saclay.inria.fr/frederic.alauzet/cours/cea2010_V3.pdf
          metriclenght = 0.5 * ( sqrt(evecAB.transpose() * m_metric * evecAB) +  sqrt(evecAB.transpose() * m2.getMetric() * evecAB));
        }
      }
      else
      {
        //TODO
      }
      return metriclenght;
    }
};


#endif // METRIC_H
