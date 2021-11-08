/*----------------------------------------------------------------------------*/
//
// Created by F. Ledoux on 2019-01-12.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_IMESHIOSERVICE_H
#define GMDS_IMESHIOSERVICE_H
/*----------------------------------------------------------------------------*/
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include <gmds/utils/Variable.h>
#include "GMDSIo_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{

    class IMeshIOService {
    public:

        enum EVariableType{
            var_int,
            var_double,
            var_double_vec,
            var_cross,
            var_cross_2D,
            var_quaternion,
            var_unknown
        };



        struct NodeInfo{
            TCellID     id;
            math::Point point;
        };

        struct EdgeInfo{
            TCellID     id;
            TCellID     node_ids[2];
        };
        struct CellInfo{
            TCellID              id;
            ECellType            type;
            std::vector<TCellID> node_ids;
        };

        struct DataInt{
            std::string name;
            std::map<TCellID ,int> values;
        };

        struct DataID{
            std::map<TCellID ,TCellID> values;
        };

        struct DataReal{
            std::string name;
            std::map<TCellID ,double> values;
        };

        struct DataVector{
            std::string name;
            std::map<TCellID ,math::Vector3d> values;
        };

        static EVariableType getType(const VariableItf* AVar);

        virtual void createNodes(const std::vector<double>& AX,
                                 const std::vector<double>& AY,
                                 const std::vector<double>& AZ)=0;

        virtual void createEdge(const TCellID& AID1,
                                const TCellID& AID2)=0;

        virtual void createTriangle(const TCellID& AID1,
                                    const TCellID& AID2,
                                    const TCellID& AID3)=0;

        virtual void createQuad(const TCellID& AID1,
                                const TCellID& AID2,
                                const TCellID& AID3,
                                const TCellID& AID4)=0;

        virtual void createPolygon(const std::vector<TCellID> &ANodes)=0;

        virtual void createTet(const TCellID& AID1,
                               const TCellID& AID2,
                               const TCellID& AID3,
                               const TCellID& AID4)=0;

        virtual void createHex(const TCellID& AID1,
                               const TCellID& AID2,
                               const TCellID& AID3,
                               const TCellID& AID4,
                               const TCellID& AID5,
                               const TCellID& AID6,
                               const TCellID& AID7,
                               const TCellID& AID8)=0;

        virtual void createPyramid(const TCellID& AID1,
                                   const TCellID& AID2,
                                   const TCellID& AID3,
                                   const TCellID& AID4,
                                   const TCellID& AID5)=0;

        virtual void getNodes(std::vector<NodeInfo>& AInfo)=0;
        virtual void getEdges(std::vector<EdgeInfo>& AInfo)=0;
        virtual void getFaces(std::vector<CellInfo>& AInfo)=0;
        virtual void getRegions(std::vector<CellInfo>& AInfo)=0;
        virtual void getDataNodes(DataID& ADataID,
                                  std::vector<DataInt>& ADataInt,
                                  std::vector<DataReal>& ADataReal,
                                  std::vector<DataVector>& ADataVec)=0;

        virtual void getDataEdges(DataID& ADataID,
                                  std::vector<DataInt>& ADataInt,
                                  std::vector<DataReal>& ADataReal,
                                  std::vector<DataVector>& ADataVec)=0;
        virtual void getDataFaces(DataID& ADataID,
                                  std::vector<DataInt>& ADataInt,
                                  std::vector<DataReal>& ADataReal,
                                  std::vector<DataVector>& ADataVec)=0;
        virtual void getDataRegions(DataID& ADataID,
                                    std::vector<DataInt>& ADataInt,
                                    std::vector<DataReal>& ADataReal,
                                    std::vector<DataVector>& ADataVec)=0;

        virtual void addDataIntNodes(DataInt& AData) = 0;
        virtual void addDataRealNodes(DataReal& AData) = 0;
        virtual void addDataVectorNodes(DataVector& AData) = 0;

        virtual void addDataIntEdges(DataInt& AData) = 0;
        virtual void addDataRealEdges(DataReal& AData) = 0;
        virtual void addDataVectorEdges(DataVector& AData) = 0;

        virtual void addDataIntFaces(DataInt& AData) = 0;
        virtual void addDataRealFaces(DataReal& AData) = 0;
        virtual void addDataVectorFaces(DataVector& AData) = 0;

        virtual void addDataIntRegions(DataInt& AData) = 0;
        virtual void addDataRealRegions(DataReal& AData) = 0;
        virtual void addDataVectorRegions(DataVector& AData) = 0;
    };
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_IMESHIOSERVICE_H
/*----------------------------------------------------------------------------*/
