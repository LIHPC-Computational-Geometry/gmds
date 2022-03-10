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

      const M& getMetric();

      M interpolateMetric(const M& m1, const M& m2, const double t);

      M intersectionMetric(const M& m1, const M& m2);

      double metricDist(const math::Vector3d& coordA, const math::Vector3d& coordB);

      void setMetric(const M& metric);

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
    const M& Metric<M>::getMetric()
    {
      return m_metric;
    }

    template<typename M>
    void Metric<M>::setMetric(const M& metric)
    {
      m_metric = metric;
    }

    template<typename M>
    M Metric<M>::interpolateMetric(const M& m1, const M& m2, const double t)
    {
     //TODO
    }

    template<typename M>
    M Metric<M>::intersectionMetric(const M& m1, const M& m2)
    {
      //TODO
    }

    template<typename M>
    double Metric<M>::metricDist(const math::Vector3d& coordA, const math::Vector3d& coordB)
    {
      double dist = 0.0;
      if(std::is_same<M, Eigen::Matrix3d>::value == true)
      {
        math::Vector3d vecAB = coordB - coordA;
        Eigen::Vector3d evecAB = Eigen::Vector3d(vecAB.X(), vecAB.Y(), vecAB.Z());
        dist = sqrt(evecAB.transpose() * m_metric * evecAB);
      }
      else
      {
        //TODO
      }
      return dist;
    }
};


#endif // METRIC_H
