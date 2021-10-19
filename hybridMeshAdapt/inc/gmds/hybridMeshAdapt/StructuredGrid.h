#ifndef STRUCTURED_GRID_H_
#define STRUCTURED_GRID_H_
/******************************************************************************/
#include <unordered_map>
/******************************************************************************/
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
/******************************************************************************/

namespace gmds
{
  namespace hybrid
  {
    class StructuredGrid
    {
      public :

      // hash functuion for the unordered map...
      struct hash_function
      {
        template<class T>
        std::size_t operator()(const T& a) const
        {
          const std::size_t ha0(std::hash<int>()(a[0]));
          const std::size_t ha1(std::hash<int>()(a[1]));
          const std::size_t ha2(std::hash<int>()(a[2]));
          return  ha0 ^ ha1 ^ ha2;
        }
      };

        StructuredGrid() = delete;

        StructuredGrid(const math::Point& firstPoint, const std::vector<math::Vector3d>& base, const math::Vector3d& sizeGrid, const math::Vector3i& subdivisionXYZ);

        ~StructuredGrid(){}

        const std::unordered_map<math::Vector3i, math::Point, hash_function> & getUnorderedMap() const {return m_unorderedGridCoordinate;}

        math::Vector3i getPreviousNode(const math::Vector3i& currentNodeinGrid) const ;

        math::Vector3i const & getSubdivision() const {return m_subdivisionXYZ;}

      private:

        std::unordered_map<math::Vector3i, math::Point, hash_function> m_unorderedGridCoordinate;

        math::Vector3i m_subdivisionXYZ;
    };
  }
}
#endif //STRUCTURED_GRID_H_
