/*----------------------------------------------------------------------------*/
/** \file    IReader.h
 *  \author  F. LEDOUX
 *  \date    03/17/2009
 */
/*----------------------------------------------------------------------------*/
#include <fstream>
/*----------------------------------------------------------------------------*/
#include <gmds/io/IReader.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
IReader::IReader(IMeshIOService* AMeshIOService,
				 const MeshModel& ACellModel,
				 const MeshModel& ADataModel)
:m_mesh_service(AMeshIOService),
 m_cell_model(ACellModel),
 m_data_model(ADataModel),
 m_stream(nullptr)
{}
/*----------------------------------------------------------------------------*/
IReader::~IReader()
{
	delete m_stream;
}
/*----------------------------------------------------------------------------*/
void IReader::setCellOptions(const gmds::MeshModel &AModel) {
    m_cell_model = AModel;
}
/*----------------------------------------------------------------------------*/
void IReader::setDataOptions(const gmds::MeshModel &AModel) {
    m_data_model=AModel;
}
/*----------------------------------------------------------------------------*/
void IReader::read(const std::string &AFileName)
{
	//===============================================================
	// First, we check if we can read the file
    m_stream= new std::ifstream(AFileName.c_str(),std::ios::in);
    if(!m_stream){
        std::string mess ="Impossible to read file "+AFileName;
        throw GMDSException(mess);
    }
    //===============================================================
    //then we check if we have some incompatibilities between cells
    //and data to read.
    if(m_data_model.has(N) && !m_cell_model.has(N)){
        std::string mess ="Impossible to get node data without reading nodes";
        throw GMDSException(mess);
    }
    if(m_data_model.has(E) && !m_cell_model.has(E)){
        std::string mess ="Impossible to get edge data without reading edges";
        throw GMDSException(mess);
    }
    if(m_data_model.has(F) && !m_cell_model.has(F)){
        std::string mess ="Impossible to get face data without reading faces";
        throw GMDSException(mess);
    }
    if(m_data_model.has(R) && !m_cell_model.has(R)){
        std::string mess ="Impossible to get region data without reading regions";
        throw GMDSException(mess);
    }
    if(!m_cell_model.has(N)){
        std::string mess ="Impossible to read a mesh file without reading nodes";
        throw GMDSException(mess);
    }
	//===============================================================
	//then we check if the file respects some format usage (delegated)
	//to child classes
	if(!preCheckFormat()){
		std::string mess ="Impossible to read file "+AFileName;
		mess+=" because of format issue";
		throw GMDSException(mess);
	}

	readNodes();

    if(m_cell_model.has(E)) {
        readEdges();
    }

    if(m_cell_model.has(F)) {
        readFaces();
    }

    if(m_cell_model.has(R)) {
        readRegions();
    }


    if(m_data_model.has(N)) {
        readDataNodes();
    }

    if(m_data_model.has(E)) {
        readDataEdges();
    }

    if(m_data_model.has(F)) {
        readDataFaces();
    }

    if(m_data_model.has(R)) {
        readDataRegions();
    }

    m_stream->close();
}
/*----------------------------------------------------------------------------*/
void IReader::readNodes() {
	std::cout<<"Warning - Node reading not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IReader::readEdges() {
	std::cout<<"Warning - Edge reading not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IReader::readFaces() {
	std::cout<<"Warning - Face reading not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IReader::readRegions() {
	std::cout<<"Warning - Region reading not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IReader::readDataNodes() {
    std::cout<<"Warning - Data node reading not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IReader::readDataEdges() {
    std::cout<<"Warning - Data edge reading not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IReader::readDataFaces() {
    std::cout<<"Warning - Data face reading not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IReader::readDataRegions() {
    std::cout<<"Warning - Data region reading not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
bool IReader::moveStreamOntoFirst(const std::string &AString){
	//go to the beginning of the stream
	m_stream->clear();
    m_stream->seekg(0,std::ios::beg);
	std::string str;
	bool found = false;
	while (!found && (*m_stream)>>str)  {
        found = (std::string::npos != str.find(AString));
	}

	return found;
}
//
///*----------------------------------------------------------------------------*/
//void IReader::updateMeshIDContainers()
//{
//	m_mesh.updateIDContainers();
//}
///*----------------------------------------------------------------------------*/
//void IReader::specifyMaxNodeID(const TInt& AID)
//{
//	m_mesh.clearAndResizeNodeIDContainer(AID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::specifyMaxEdgeID(const TInt& AID)
//{
//	m_mesh.clearAndResizeEdgeIDContainer(AID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::specifyMaxFaceID(const TInt& AID)
//{
//	m_mesh.clearAndResizeFaceIDContainer(AID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::specifyMaxRegionID(const TInt& AID)
//{
//	m_mesh.clearAndResizeRegionIDContainer(AID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::newNode(const TCoord& AX, const TCoord& AY,
//					  const TCoord& AZ, const TCellID& AGID)
//{
//	m_mesh.newNodeWithID(AX,AY,AZ,AGID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::newNode(const TCoord& AX, const TCoord& AY,
//					  const TCellID& AGID)
//{
//	m_mesh.newNodeWithID(AX,AY,AGID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::newEdge(const TCellID& AN1, const TCellID& AN2,
//					  const TCellID& AGID)
//{
//	m_mesh.newEdgeWithID(AN1,AN2,AGID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::newTriangle(const TCellID& AN1, const TCellID& AN2,
//						  const TCellID& AN3, const TCellID& AGID)
//{
//	m_mesh.newTriangleWithID(AN1,AN2,AN3,AGID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::newQuad(const TCellID& AN1, const TCellID& AN2,
//					  const TCellID& AN3, const TCellID& AN4,
//					  const TCellID& AGID)
//{
//	m_mesh.newQuadWithID(AN1,AN2,AN3,AN4,AGID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::newPolygon(std::vector<TCellID>& AIDs, const TCellID& AGID)
//{
//	m_mesh.newPolygonWithID(AIDs,AGID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::newFace(std::vector<TCellID>& AIDs, const TCellID& AGID)
//{
//	m_mesh.newFaceWithID(AIDs,AGID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::newTet(const TCellID& AN1, const TCellID& AN2,
//					 const TCellID& AN3, const TCellID& AN4,
//					 const TCellID& AGID)
//{
//	m_mesh.newTetWithID(AN1,AN2,AN3,AN4,AGID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::newPyramid(const TCellID& AN1, const TCellID& AN2,
//						 const TCellID& AN3, const TCellID& AN4,
//						 const TCellID& AN5, const TCellID& AGID)
//{
//	m_mesh.newPyramidWithID(AN1,AN2,AN3,AN4,AN5,AGID);
//}
///*----------------------------------------------------------------------------*/
//void IReader::newHex(const TCellID& AN1, const TCellID& AN2,
//					 const TCellID& AN3, const TCellID& AN4,
//					 const TCellID& AN5, const TCellID& AN6,
//					 const TCellID& AN7, const TCellID& AN8,
//					 const TCellID& AGID)
//{
//	m_mesh.newHexWithID(AN1,AN2,AN3,AN4,AN5,AN6,AN7,AN8,AGID);
//}
///*----------------------------------------------------------------------------*/
