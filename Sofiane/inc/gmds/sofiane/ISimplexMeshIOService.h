#ifndef ISIMPLEXMESH_IO_SERVICE_H_
#define ISIMPLEXMESH_IO_SERVICE_H_
/******************************************************************************/
#include <gmds/io/IMeshIOService.h>
/******************************************************************************/
#include <gmds/sofiane/SimplicesNode.h>
/******************************************************************************/
namespace gmds
{
    namespace hybrid
    {
        class SimplexMesh;
    }

    class ISimplexMeshIOService : public IMeshIOService
    {
    public:

      ISimplexMeshIOService(hybrid::SimplexMesh* simplexMesh);

       void createNodes(const std::vector<double>& AX,
                               const std::vector<double>& AY,
                               const std::vector<double>& AZ);

       void createEdge(const TCellID& AID1,
                              const TCellID& AID2){}

       void createTriangle(const TCellID& AID1,
                                  const TCellID& AID2,
                                  const TCellID& AID3);

      /*peut etre a implementer plus tard...*/
       void createQuad(const TCellID& AID1,
                              const TCellID& AID2,
                              const TCellID& AID3,
                              const TCellID& AID4){}

       void createPolygon(const std::vector<TCellID> &ANodes){};

       void createTet(const TCellID& AID1,
                             const TCellID& AID2,
                             const TCellID& AID3,
                             const TCellID& AID4);

       void createHex(const TCellID& AID1,
                             const TCellID& AID2,
                             const TCellID& AID3,
                             const TCellID& AID4,
                             const TCellID& AID5,
                             const TCellID& AID6,
                             const TCellID& AID7,
                             const TCellID& AID8){}

       void createPyramid(const TCellID& AID1,
                                 const TCellID& AID2,
                                 const TCellID& AID3,
                                 const TCellID& AID4,
                                 const TCellID& AID5){}

       void getNodes(std::vector<NodeInfo>& AInfo);
       void getEdges(std::vector<EdgeInfo>& AInfo);
       void getFaces(std::vector<CellInfo>& AInfo);
       void getRegions(std::vector<CellInfo>& AInfo);
       void getDataNodes(DataID& ADataID,
                                std::vector<DataInt>& ADataInt,
                                std::vector<DataReal>& ADataReal,
                                std::vector<DataVector>& ADataVec);

       void getDataEdges(DataID& ADataID,
                                std::vector<DataInt>& ADataInt,
                                std::vector<DataReal>& ADataReal,
                                std::vector<DataVector>& ADataVec);
       void getDataFaces(DataID& ADataID,
                                std::vector<DataInt>& ADataInt,
                                std::vector<DataReal>& ADataReal,
                                std::vector<DataVector>& ADataVec);
       void getDataRegions(DataID& ADataID,
                                  std::vector<DataInt>& ADataInt,
                                  std::vector<DataReal>& ADataReal,
                                  std::vector<DataVector>& ADataVec);

       void addDataIntNodes(DataInt& AData) ;
       void addDataRealNodes(DataReal& AData) ;
       void addDataVectorNodes(DataVector& AData) ;

       void addDataIntEdges(DataInt& AData) ;
       void addDataRealEdges(DataReal& AData) ;
       void addDataVectorEdges(DataVector& AData) ;

       void addDataIntFaces(DataInt& AData) ;
       void addDataRealFaces(DataReal& AData) ;
       void addDataVectorFaces(DataVector& AData) ;

       void addDataIntRegions(DataInt& AData) ;
       void addDataRealRegions(DataReal& AData) ;
       void addDataVectorRegions(DataVector& AData) ;

    private:

      hybrid::SimplexMesh*        m_simplex_mesh;
    };
}

#endif // ISIMPLEXMESH_IO_SERVICE_H_
