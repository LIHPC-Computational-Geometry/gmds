/*----------------------------------------------------------------------------*/
#include <sstream>
/*----------------------------------------------------------------------------*/
#include "gmds/io/IWriter.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
IWriter::IWriter(IMeshIOService* AMeshIOService,
                 const MeshModel& ACellModel,
                 const MeshModel& ADataModel)
        :m_mesh_service(AMeshIOService),
         m_cell_model(ACellModel),
         m_data_model(ADataModel),
         m_stream(nullptr)
{}
/*----------------------------------------------------------------------------*/
IWriter::~IWriter()
{

	delete m_stream;
}
/*----------------------------------------------------------------------------*/
void IWriter::setCellOptions(const gmds::MeshModel &AModel) {
    m_cell_model = AModel;
}
/*----------------------------------------------------------------------------*/
void IWriter::setDataOptions(const gmds::MeshModel &AModel) {
    m_data_model=AModel;
}
/*----------------------------------------------------------------------------*/
void IWriter::write(const std::string &AFileName)
{

    //===============================================================
    //then we check if we have some incompatibilities between cells
    //and data to write.
    if(m_data_model.has(N) && !m_cell_model.has(N)){
        std::string mess ="Impossible to write node data without writing nodes";
        throw GMDSException(mess);
    }
    if(m_data_model.has(E) && !m_cell_model.has(E)){
        std::string mess ="Impossible to write edge data without writing edges";
        throw GMDSException(mess);
    }
    if(m_data_model.has(F) && !m_cell_model.has(F)){
        std::string mess ="Impossible to write face data without writing faces";
        throw GMDSException(mess);
    }
    if(m_data_model.has(R) && !m_cell_model.has(R)){
        std::string mess ="Impossible to write region data without writing regions";
        throw GMDSException(mess);
    }
    if(!m_cell_model.has(N)){
        std::string mess ="Impossible to write a mesh file without writing nodes";
        throw GMDSException(mess);
    }
    //===============================================================

    //===============================================================
    // First, we check if we can write the file
    std::stringstream file_name;
    file_name << AFileName;

    try{
        initialize(file_name.str());
    } catch (GMDSException& e){
        throw e;
    }
    writeNodes();


    if(m_cell_model.has(E)) {writeEdges();}
    if(m_cell_model.has(F)) {writeFaces();}
    if(m_cell_model.has(R)) {writeRegions();}
    if(m_data_model.has(N)) {writeDataNodes();}
    if(m_data_model.has(E)) {writeDataEdges();}
    if(m_data_model.has(F)) {writeDataFaces();}
    if(m_data_model.has(R)) {writeDataRegions();}

    finalize();
}
/*----------------------------------------------------------------------------*/
void IWriter::writeNodes() {
    std::cout<<"Warning - Node writing not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IWriter::writeEdges() {
    std::cout<<"Warning - Edge writing not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IWriter::writeFaces() {
    std::cout<<"Warning - Face writing not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IWriter::writeRegions() {
    std::cout<<"Warning - Region writing not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IWriter::writeDataNodes() {
    std::cout<<"Warning - Data node writing not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IWriter::writeDataEdges() {
    std::cout<<"Warning - Data edge writing not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IWriter::writeDataFaces() {
    std::cout<<"Warning - Data face writing not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void IWriter::writeDataRegions() {
    std::cout<<"Warning - Data region writing not implemented for this data format"<<std::endl;
}
/*----------------------------------------------------------------------------*/
