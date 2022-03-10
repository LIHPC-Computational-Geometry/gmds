/******************************************************************************/
#include <gmds/hybridMeshAdapt/StructuredGrid.h>
/******************************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace math;
/******************************************************************************/
StructuredGrid::StructuredGrid(const Point& firstPoint, const std::vector<Vector3d>& base, const Vector3d& sizeGrid, const Vector3i& subdivisionXYZ)
{
  if(base.size() == 3)
  {
    if(subdivisionXYZ[0] != 0 || subdivisionXYZ[1] != 0 || subdivisionXYZ[2] != 0)
    {
      m_subdivisionXYZ = subdivisionXYZ;
      Vector3i index3D_Point;
      unsigned int numbersPoints =  subdivisionXYZ[0] * subdivisionXYZ[1] * subdivisionXYZ[2];
      for(unsigned int i = 0; i < subdivisionXYZ[0]; i++)
      {
        for(unsigned int j = 0; j < subdivisionXYZ[1]; j++)
        {
          for(unsigned int k = 0; k < subdivisionXYZ[2]; k++)
          {
            index3D_Point[0] = i;
            index3D_Point[1] = j;
            index3D_Point[2] = k;

            Vector3d Ni =  (double)i / (double)subdivisionXYZ[0] * sizeGrid[0] * base[0];
            Vector3d Nj =  (double)j / (double)subdivisionXYZ[1] * sizeGrid[1] * base[1];
            Vector3d Nk =  (double)k / (double)subdivisionXYZ[2] * sizeGrid[2] * base[2];

            Point point = firstPoint + Point(Ni[0] + Nj[0] + Nk[0], Ni[1] + Nj[1] + Nk[1], Ni[2] + Nj[2] + Nk[2]);
            m_unorderedGridCoordinate[index3D_Point] = point;
          }
        }
      }
    }
    else
    {
      /*TODO THROW*/
      std::cout << "can't subdivide the grid because of the subdivisionXYZ == 0 !! " << std::endl;
    }
  }
  else
  {
    /*TODO THROW*/
    std::cout << "base size different of 3 !! " << std::endl;
  }
}



Vector3i StructuredGrid::getPreviousNode(const Vector3i& currentNodeinGrid) const
{
  Vector3i previousNode(0, 0, 0);
  if((currentNodeinGrid[0] < 1 ) || (currentNodeinGrid[0] > m_subdivisionXYZ[0] - 1) ||
      (currentNodeinGrid[1] < 1 ) || (currentNodeinGrid[1] > m_subdivisionXYZ[1] - 1) ||
        (currentNodeinGrid[2] < 1 ) || (currentNodeinGrid[2] > m_subdivisionXYZ[2] - 1))
        {
          previousNode = currentNodeinGrid - Vector3i(-1, -1, -1);
        }
  else
  {
    /*TODO throw...*/
    std::cout << " can not get Previous node of : " << currentNodeinGrid << " one or more of the component is out of border" << std::endl;
  }

  return std::move(previousNode);
}
