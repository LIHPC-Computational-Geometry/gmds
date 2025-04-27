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
#include "GMDSIo_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds {

    class  Mesh;

    class GMDSIo_API IGMeshIOService : public IMeshIOService {
    public:
        explicit IGMeshIOService(Mesh *AMesh);

        void createNodes(const std::vector<double> &AX,
                         const std::vector<double> &AY,
                         const std::vector<double> &AZ) override;

        void createEdge(const TCellID &AID1,
                        const TCellID &AID2) override;

        void createTriangle(const TCellID &AID1,
                            const TCellID &AID2,
                            const TCellID &AID3) override;

        void createQuad(const TCellID &AID1,
                        const TCellID &AID2,
                        const TCellID &AID3,
                        const TCellID &AID4) override;


        void createPolygon(const std::vector<TCellID> &ANodes) override;


        void createTet(const TCellID &AID1,
                       const TCellID &AID2,
                       const TCellID &AID3,
                       const TCellID &AID4) override;

        void createHex(const TCellID &AID1,
                       const TCellID &AID2,
                       const TCellID &AID3,
                       const TCellID &AID4,
                       const TCellID &AID5,
                       const TCellID &AID6,
                       const TCellID &AID7,
                       const TCellID &AID8) override;

        void createPyramid(const TCellID &AID1,
                           const TCellID &AID2,
                           const TCellID &AID3,
                           const TCellID &AID4,
                           const TCellID &AID5) override;

        void getNodes(std::vector<NodeInfo> &AInfo) override;

        void getEdges(std::vector<EdgeInfo> &AInfo) override;

        void getFaces(std::vector<CellInfo> &AInfo) override;

        void getRegions(std::vector<CellInfo> &AInfo) override;


        void getDataNodes(DataID &ADataID,
                          std::vector<DataInt> &ADataInt,
                          std::vector<DataReal> &ADataReal,
                          std::vector<DataVector>& ADataVec) override;

        void getDataEdges(DataID &ADataID,
                          std::vector<DataInt> &ADataInt,
                          std::vector<DataReal> &ADataReal,
                          std::vector<DataVector>& ADataVec) override;

        void getDataFaces(DataID &ADataID,
                          std::vector<DataInt> &ADataInt,
                          std::vector<DataReal> &ADataReal,
                          std::vector<DataVector>& ADataVec) override;

        void getDataRegions(DataID &ADataID,
                            std::vector<DataInt> &ADataInt,
                            std::vector<DataReal> &ADataReal,
                            std::vector<DataVector>& ADataVec) override;

        void addDataIntNodes(DataInt &AData) override;
        void addDataRealNodes(DataReal &AData) override;
        void addDataVectorNodes(DataVector& AData) override;

        void addDataIntEdges(DataInt& AData) override;
        void addDataRealEdges(DataReal& AData) override;
        void addDataVectorEdges(DataVector& AData) override;

        void addDataIntFaces(DataInt& AData) override;
        void addDataRealFaces(DataReal& AData) override;
        void addDataVectorFaces(DataVector& AData) override;

        void addDataIntRegions(DataInt& AData) override;
        void addDataRealRegions(DataReal& AData) override;
        void addDataVectorRegions(DataVector& AData) override;


    private:
        /** Mesh instance we work with */
        Mesh *m_mesh;
    };
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_IGMESHIOSERVICE_H
/*----------------------------------------------------------------------------*/
