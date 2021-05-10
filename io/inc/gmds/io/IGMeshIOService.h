/*----------------------------------------------------------------------------*/
//
// Created by F. Ledoux on 2019-01-12.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_IGMESHIOSERVICE_H
#define GMDS_IGMESHIOSERVICE_H
/*----------------------------------------------------------------------------*/
#include <map>
/*----------------------------------------------------------------------------*/
#include "IMeshIOService.h"
/*----------------------------------------------------------------------------*/
namespace gmds {

    class Mesh;

    class IGMeshIOService : public IMeshIOService {
    public:
        IGMeshIOService(Mesh *AMesh);

        void createNodes(const std::vector<double> &AX,
                         const std::vector<double> &AY,
                         const std::vector<double> &AZ);

        void createEdge(const TCellID &AID1,
                        const TCellID &AID2);

        void createTriangle(const TCellID &AID1,
                            const TCellID &AID2,
                            const TCellID &AID3);

        void createQuad(const TCellID &AID1,
                        const TCellID &AID2,
                        const TCellID &AID3,
                        const TCellID &AID4);


        void createPolygon(const std::vector<TCellID> &ANodes);


        void createTet(const TCellID &AID1,
                       const TCellID &AID2,
                       const TCellID &AID3,
                       const TCellID &AID4);

        void createHex(const TCellID &AID1,
                       const TCellID &AID2,
                       const TCellID &AID3,
                       const TCellID &AID4,
                       const TCellID &AID5,
                       const TCellID &AID6,
                       const TCellID &AID7,
                       const TCellID &AID8);

        void createPyramid(const TCellID &AID1,
                           const TCellID &AID2,
                           const TCellID &AID3,
                           const TCellID &AID4,
                           const TCellID &AID5);

        void getNodes(std::vector<NodeInfo> &AInfo);

        void getEdges(std::vector<EdgeInfo> &AInfo);

        void getFaces(std::vector<CellInfo> &AInfo);

        void getRegions(std::vector<CellInfo> &AInfo);


        void getDataNodes(DataID &ADataID,
                          std::vector<DataInt> &ADataInt,
                          std::vector<DataReal> &ADataReal,
                          std::vector<DataVector>& ADataVec);

        void getDataEdges(DataID &ADataID,
                          std::vector<DataInt> &ADataInt,
                          std::vector<DataReal> &ADataReal,
                          std::vector<DataVector>& ADataVec);

        void getDataFaces(DataID &ADataID,
                          std::vector<DataInt> &ADataInt,
                          std::vector<DataReal> &ADataReal,
                          std::vector<DataVector>& ADataVec);

        void getDataRegions(DataID &ADataID,
                            std::vector<DataInt> &ADataInt,
                            std::vector<DataReal> &ADataReal,
                            std::vector<DataVector>& ADataVec);

        void addDataIntNodes(DataInt &AData);
        void addDataRealNodes(DataReal &AData);
        void addDataVectorNodes(DataVector& AData);

        void addDataIntEdges(DataInt& AData);
        void addDataRealEdges(DataReal& AData);
        void addDataVectorEdges(DataVector& AData);

        void addDataIntFaces(DataInt& AData);
        void addDataRealFaces(DataReal& AData);
        void addDataVectorFaces(DataVector& AData);

        void addDataIntRegions(DataInt& AData);
        void addDataRealRegions(DataReal& AData);
        void addDataVectorRegions(DataVector& AData);


    private:
        /** Mesh instance we work with */
        Mesh *m_mesh;
    };
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_IGMESHIOSERVICE_H
/*----------------------------------------------------------------------------*/
