/*----------------------------------------------------------------------------*/
/** \file    VTKReader.cpp
 *  \author  F. LEDOUX
 *  \date    09/11/2008
 */
/*----------------------------------------------------------------------------*/
// headers of GMDS files
#include <gmds/io/VTKReader.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
const int VTKReader::m_nb_max_node_per_cell=12;
/*----------------------------------------------------------------------------*/
VTKReader::VTKReader(IMeshIOService *AMeshService)
        :IReader(AMeshService)//,m_nb_imported_faces(0), m_nb_imported_regions(0)
{}
/*----------------------------------------------------------------------------*/
VTKReader::~VTKReader()
{}
/*----------------------------------------------------------------------------*/
bool VTKReader::preCheckFormat() {


    if(!moveStreamOntoFirst("UNSTRUCTURED_GRID")) {
        return false;
    }


    //First we look at cell types. We handle only triangles and quads
    if(!moveStreamOntoFirst("CELL_TYPES")) {
        std::string mess = "VTK read error: no CELL_TYPES keyword found";
        throw GMDSException(mess);
    }

    TInt nb_cell_types;
    *m_stream >> nb_cell_types;

    m_cell_types.resize(nb_cell_types);
    for (auto i = 0; i < nb_cell_types; i++) {

        int cell_type;
        *m_stream >> cell_type;

        switch(cell_type) {
            case 3://LINE
                m_cell_types[i]=GMDS_EDGE;
                break;
            case 5://TRIANGLE
                m_cell_types[i]=GMDS_TRIANGLE;
                break;
            case 9://QUAD
                m_cell_types[i]=GMDS_QUAD;
                break;
            case 7://POLYGON
                m_cell_types[i]=GMDS_POLYGON;
                break;
            case 10://TET
                m_cell_types[i]=GMDS_TETRA;
                break;
		      case 11://VOXEL
			       m_cell_types[i]=GMDS_VOXEL;
			       break;
            case 12://HEX
                m_cell_types[i]=GMDS_HEX;
                break;
            case 14://PYRAMID
                m_cell_types[i]=GMDS_PYRAMID;
                break;
            default:
                m_cell_types[i]=GMDS_UNKNOWN;
                break;
        }
    }


    return true;

}
//    do{
//        input >> current_word;
//    } while ( (!isPolydata && current_word != "CELLS") || (isPolydata && current_word != "POLYGONS"));
//
//    readCells(input,isPolydata);
//
//    readNodesData(AFileName);
//
//    if (m_nb_imported_faces != 0 && m_nb_imported_regions == 0)
//        readFacesData(AFileName);
//    else if (m_nb_imported_faces == 0 && m_nb_imported_regions != 0)
//        readRegionsData(AFileName);
//    else
//        throw GMDSException("Impossible to read a VTK file containing both faces and regions");
//
//    input.close();
//
//
//}

/*----------------------------------------------------------------------------*/
void VTKReader::readNodes()
{
    if(!moveStreamOntoFirst("POINTS")) {
        std::string mess = "VTK read error: no POINTS keyword found";
        throw GMDSException(mess);
    }


    TInt nb_nodes;
    std::string coord_type;
    *m_stream >> nb_nodes >> coord_type;

    std::vector<double> x,y,z;
    x.resize(nb_nodes);
    y.resize(nb_nodes);
    z.resize(nb_nodes);

    for (int i = 0; i < nb_nodes; i++){
        *m_stream >> x[i] >> y[i] >> z[i];
    }

    m_mesh_service->createNodes(x, y, z);

}
/*----------------------------------------------------------------------------*/
void VTKReader::readEdges() {

    //Now we look for cell nodes
     if (!moveStreamOntoFirst("CELLS")) {
        std::string mess = "VTK read error: no CELLS keyword found";
        throw GMDSException(mess);
    }
    TInt nb_cells, nb_values;
    *m_stream >> nb_cells;
    *m_stream >> nb_values;


    //We use a temporary vector with a max size to avoid
    //multiple vector allocations and resizings.
    std::vector<TCellID> nodes;

    nodes.resize(m_nb_max_node_per_cell);
    for (auto i = 0; i < nb_cells; i++) {

        int nb_nodes;
        *m_stream >> nb_nodes;

        for (int j = 0; j < nb_nodes; j++) {
            int val;
            *m_stream >> val;
            nodes[j] = val;
        }

        if (m_cell_types[i] == GMDS_EDGE) {
            m_mesh_service->createEdge(nodes[0],
                                       nodes[1]);
        }
    }

}
/*----------------------------------------------------------------------------*/
void VTKReader::readFaces()
{

    //Now we look for cell nodes
    if (!moveStreamOntoFirst("CELLS")) {
        std::string mess = "VTK read error: no CELLS keyword found";
        throw GMDSException(mess);
    }
    TInt nb_cells, nb_values;
    *m_stream >> nb_cells;
    *m_stream >> nb_values;


    //We use a temporary vector with a max size to avoid
    //multiple vector allocations and resizings.
    std::vector<TCellID> nodes;

    nodes.resize(m_nb_max_node_per_cell);
    for (auto i = 0; i < nb_cells; i++) {

        int nb_nodes;
        *m_stream >> nb_nodes;
        if (nb_nodes > m_nb_max_node_per_cell) {
            std::string mess = "VTKReader is unable to create faces with more than ";
            mess += std::to_string(m_nb_max_node_per_cell);
            mess += ", actual number is ";
            mess += std::to_string(nb_nodes);

            throw GMDSException(mess);
        }
        for (int j = 0; j < nb_nodes; j++) {
            int val;
            *m_stream >> val;
            nodes[j] = val;
        }

        if (m_cell_types[i] == GMDS_TRIANGLE) {
            m_mesh_service->createTriangle(nodes[0],
                                           nodes[1],
                                           nodes[2]);
        } else if (m_cell_types[i] == GMDS_QUAD) {
            m_mesh_service->createQuad(nodes[0],
                                       nodes[1],
                                       nodes[2],
                                       nodes[3]);
        } else if(m_cell_types[i]==GMDS_POLYGON){

            std::vector<TCellID> poly_nodes;
            poly_nodes.resize(nb_nodes);
            for (auto i_nodes = 0; i_nodes < nb_nodes; i_nodes++) {
                poly_nodes[i_nodes] = nodes[i_nodes];
            }
            m_mesh_service->createPolygon(poly_nodes);
        }
    }
}

/*----------------------------------------------------------------------------*/
void VTKReader::readRegions()
{

    //Now we look for cell nodes
    if (!moveStreamOntoFirst("CELLS")) {
        std::string mess = "VTK read error: no CELLS keyword found";
        throw GMDSException(mess);
    }
    TInt nb_cells, nb_values;
    *m_stream >> nb_cells;
    *m_stream >> nb_values;


    //We use a temporary vector with a max size to avoid
    //multiple vector allocations and resizings.
    std::vector<TCellID> nodes;

    nodes.resize(m_nb_max_node_per_cell);
    for (auto i = 0; i < nb_cells; i++) {

        int nb_nodes;
        *m_stream >> nb_nodes;
        if (nb_nodes > m_nb_max_node_per_cell) {
            std::string mess = "VTKReader is unable to create regions with more than ";
            mess += std::to_string(m_nb_max_node_per_cell);
            mess += ", actual number is ";
            mess += std::to_string(nb_nodes);

            throw GMDSException(mess);
        }
        for (int j = 0; j < nb_nodes; j++) {
            int val;
            *m_stream >> val;
            nodes[j] = val;
        }

        if (m_cell_types[i] == GMDS_TETRA) {
            m_mesh_service->createTet(nodes[0],
                                      nodes[1],
                                      nodes[2],
                                      nodes[3]);
		  } else if (m_cell_types[i] == GMDS_VOXEL) {
			   m_mesh_service->createHex(nodes[0],
			                             nodes[1],
			                             nodes[3],
			                             nodes[2],
			                             nodes[4],
			                             nodes[5],
			                             nodes[7],
			                             nodes[6]);
        } else if (m_cell_types[i] == GMDS_HEX) {
            m_mesh_service->createHex(nodes[0],
                                      nodes[1],
                                      nodes[2],
                                      nodes[3],
                                      nodes[4],
                                      nodes[5],
                                      nodes[6],
                                      nodes[7]);
        } else if(m_cell_types[i]==GMDS_PYRAMID){
            m_mesh_service->createPyramid(nodes[0],
                                          nodes[1],
                                          nodes[2],
                                          nodes[3],
                                          nodes[4]);
        }
    }
}
/*----------------------------------------------------------------------------*/
void VTKReader::readDataNodes(){
    //Now we look for cell nodes
    if (!moveStreamOntoFirst("POINT_DATA")) {
        std::string mess = "VTK read error: no POINT_DATA keyword found";
        throw GMDSException(mess);
    }
    TInt nb_values;
    *m_stream >> nb_values;


    //SCALARS VALUE
    std::string current_word ="";

    while (current_word != "SCALARS" &&
           current_word != "FIELD" &&
           current_word != "CELL_DATA" &&
           current_word != "VECTORS" &&
           !m_stream->eof())
    {
        *m_stream >> current_word;
        if (current_word == "SCALARS")
        {
            std::string scalar_name, scalar_type, scalar_nb;
            *m_stream >> scalar_name >> scalar_type >> scalar_nb;
            *m_stream >> current_word; // For LOOKUP_TABLE
            *m_stream >> current_word; // For default
            if (scalar_name != "GMDS_ID"){
                if (scalar_type == "int"){
                    IMeshIOService::DataInt d;
                    d.name = scalar_name;

                    for (int i = 0; i < nb_values; i++){
                        int val;
                        *m_stream >> val;
                        d.values[i] = val;
                    }
                    m_mesh_service->addDataIntNodes(d);
                }
                else if (scalar_type == "float"){
                    IMeshIOService::DataReal d;
                    d.name = scalar_name;

                    for (int i = 0; i < nb_values; i++){
                        double val;
                        *m_stream >> val;
                        d.values[i] = val;
                    }
                    m_mesh_service->addDataRealNodes(d);
                }
            }
        }
        else if (current_word == "VECTORS")
        {
            std::string vector_name, vector_type;
            *m_stream >> vector_name >> vector_type;

            // WARNING: DO NOT REMOVE
            // the next affectation is only used for avoiding loop error
            // in the global while test.
            current_word  = vector_type;

            IMeshIOService::DataVector d;
            d.name = vector_name;
            for (int i = 0; i < nb_values; i++){
                double x, y, z;
                *m_stream >> x >> y >> z;
                d.values[i] = math::Vector3d({x, y, z});
            }
            m_mesh_service->addDataVectorNodes(d);
        }
    }

}
/*----------------------------------------------------------------------------*/
void VTKReader::readDataEdges(){
    //Now we look for cell nodes
    if (!moveStreamOntoFirst("CELL_DATA")) {
        std::string mess = "VTK read error: no CELL_DATA keyword found";
        throw GMDSException(mess);
    }

    //We get the whole number of cells (edges, faces and regions)
    TInt nb_values;
    *m_stream >> nb_values;


    //SCALARS VALUE
    std::string current_word ="";
    while (current_word != "SCALARS" &&
           current_word != "FIELD" &&
           current_word != "POINT_DATA" &&
           current_word != "VECTORS" &&
           !m_stream->eof())
    {
        *m_stream >> current_word;

        if (current_word == "SCALARS") {
            std::string scalar_name, scalar_type, scalar_nb;
            *m_stream >> scalar_name >> scalar_type >> scalar_nb;
            *m_stream >> current_word; // For LOOKUP_TABLE
            *m_stream >> current_word; // For default
            if (scalar_name != "GMDS_ID"){
                if (scalar_type == "int"){
                    std::vector<int> tmp_val;
                    bool only_null = true;
                    tmp_val.reserve(nb_values);
                    int nb_true_cells = 0;
                    for (int i = 0; i < nb_values; i++){
                        int val;
                        *m_stream >> val;
                        if (m_cell_types[i]==GMDS_EDGE) {
                            tmp_val.push_back(val);
                            nb_true_cells++;
                            if(val!=0){
                                only_null=false;
                            }
                        }
                    }
                    //if only_node is true, we expect that this cell data is
                    // on a cell type that is not edge.
                    if(!only_null) {
                        IMeshIOService::DataInt d;
                        d.name = scalar_name;
                        for(int i_cell = 0; i_cell<nb_true_cells;i_cell++){
                            d.values[i_cell]=tmp_val[i_cell];
                        }
                        m_mesh_service->addDataIntEdges(d);
                    }
                }
                else if (scalar_type == "float"){
                    std::vector<double> tmp_val;
                    bool only_null = true;
                    tmp_val.reserve(nb_values);
                    int nb_true_cells = 0;
                    for (int i = 0; i < nb_values; i++){
                        double val;
                        *m_stream >> val;
                        if (m_cell_types[i]==GMDS_EDGE) {
                            tmp_val.push_back(val);
                            nb_true_cells++;
                            if(val!=0.0){
                                only_null=false;
                            }
                        }
                    }
                    //if only_node is true, we expect that this cell data is
                    // on a cell type that is not edge.
                    if(!only_null) {
                        IMeshIOService::DataReal d;
                        d.name = scalar_name;
                        for(int i_cell = 0; i_cell<nb_true_cells;i_cell++){
                            d.values[i_cell]=tmp_val[i_cell];
                        }
                        m_mesh_service->addDataRealEdges(d);
                    }
                }
            }
        }
        else if (current_word == "VECTORS")
        {
            std::string vector_name, vector_type;
            *m_stream >> vector_name >> vector_type;

            std::vector<math::Vector3d> tmp_val;
            bool only_null = true;
            tmp_val.reserve(nb_values);
            int nb_true_cells = 0;
            for (int i = 0; i < nb_values; i++){
                double x=0, y=0, z=0;
                *m_stream >> x;
                *m_stream >> y;
                *m_stream >> z;
                if (m_cell_types[i]==GMDS_EDGE) {

                    tmp_val.push_back({x,y,z});
                    nb_true_cells++;
                    if(!(x==0.0 && y==0 && z==0)){
                        only_null=false;
                    }
                }
            }
            //if only_node is true, we expect that this cell data is
            // on a cell type that is not edge.
            if(!only_null) {
                IMeshIOService::DataVector d;
                d.name = vector_name;
                for(int i_cell = 0; i_cell<nb_true_cells;i_cell++){
                    d.values[i_cell]=tmp_val[i_cell];
                }
                m_mesh_service->addDataVectorEdges(d);
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void VTKReader::readDataFaces(){
    //Now we look for cell nodes
    if (!moveStreamOntoFirst("CELL_DATA")) {
        std::string mess = "VTK read error: no CELL_DATA keyword found";
        throw GMDSException(mess);
    }

    //We get the whole number of cells (edges, faces and regions)
    TInt nb_values;
    *m_stream >> nb_values;

    //SCALARS VALUE
    std::string current_word ="";

    while (current_word != "SCALARS" &&
           current_word != "FIELD" &&
           current_word != "POINT_DATA" &&
           current_word != "VECTORS" &&
           !m_stream->eof())
    {
        *m_stream >> current_word;

        if (current_word == "SCALARS") {
            std::string scalar_name, scalar_type, scalar_nb;
            *m_stream >> scalar_name >> scalar_type >> scalar_nb;
            *m_stream >> current_word; // For LOOKUP_TABLE
            *m_stream >> current_word; // For default
            if (scalar_name != "GMDS_ID"){
                if (scalar_type == "int"){
                    std::vector<int> tmp_val;
                    bool only_null = true;
                    tmp_val.reserve(nb_values);
                    int nb_true_cells = 0;
                    for (int i = 0; i < nb_values; i++){
                        int val;
                        *m_stream >> val;
                        if (m_cell_types[i]==GMDS_QUAD ||
                            m_cell_types[i]==GMDS_TRIANGLE ||
                            m_cell_types[i]==GMDS_POLYGON ) {
                            tmp_val.push_back(val);
                            nb_true_cells++;
                            if(val!=0){
                                only_null=false;
                            }
                        }
                    }
                    //if only_node is true, we expect that this cell data is
                    // on a cell type that is not edge.
                    if(!only_null) {
                        IMeshIOService::DataInt d;
                        d.name = scalar_name;
                        for(int i_cell = 0; i_cell<nb_true_cells;i_cell++){
                            d.values[i_cell]=tmp_val[i_cell];
                        }
                        m_mesh_service->addDataIntFaces(d);
                    }
                }
                else if (scalar_type == "float"){
                    std::vector<double> tmp_val;
                    bool only_null = true;
                    tmp_val.reserve(nb_values);
                    int nb_true_cells = 0;
                    for (int i = 0; i < nb_values; i++){
                        double val;
                        *m_stream >> val;
                        if (m_cell_types[i]==GMDS_QUAD ||
                            m_cell_types[i]==GMDS_TRIANGLE ||
                            m_cell_types[i]==GMDS_POLYGON ) {
                            tmp_val.push_back(val);
                            nb_true_cells++;
                            if(val!=0.0){
                                only_null=false;
                            }
                        }
                    }
                    //if only_node is true, we expect that this cell data is
                    // on a cell type that is not edge.
                    if(!only_null) {
                        IMeshIOService::DataReal d;
                        d.name = scalar_name;
                        for(int i_cell = 0; i_cell<nb_true_cells;i_cell++){
                            d.values[i_cell]=tmp_val[i_cell];
                        }
                        m_mesh_service->addDataRealFaces(d);
                    }
                }
            }
        }
        else if (current_word == "VECTORS")
        {
            std::string vector_name, vector_type;
            *m_stream >> vector_name >> vector_type;

            std::vector<math::Vector3d> tmp_val;
            bool only_null = true;
            tmp_val.reserve(nb_values);
            int nb_true_cells = 0;
            for (int i = 0; i < nb_values; i++){
                double x=0, y=0, z=0;
                *m_stream >> x;
                *m_stream >> y;
                *m_stream >> z;
                if (m_cell_types[i]==GMDS_QUAD ||
                    m_cell_types[i]==GMDS_TRIANGLE ||
                    m_cell_types[i]==GMDS_POLYGON ) {

                    tmp_val.push_back(math::Vector3d({x,y,z}));
                    nb_true_cells++;
                    if(!(x==0.0 && y==0 && z==0)){
                        only_null=false;
                    }
                }
            }
            //if only_node is true, we expect that this cell data is
            // on a cell type that is not edge.
            if(!only_null) {
                IMeshIOService::DataVector d;
                d.name = vector_name;
                for(int i_cell = 0; i_cell<nb_true_cells;i_cell++){
                    d.values[i_cell]=tmp_val[i_cell];
                }
                m_mesh_service->addDataVectorFaces(d);
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void VTKReader::readDataRegions(){
    //Now we look for cell nodes
    if (!moveStreamOntoFirst("CELL_DATA")) {
        std::string mess = "VTK read error: no CELL_DATA keyword found";
        throw GMDSException(mess);
    }

    //We get the whole number of cells (edges, faces and regions)
    TInt nb_values;
    *m_stream >> nb_values;


    //SCALARS VALUE
    std::string current_word ="";

    while (current_word != "SCALARS" &&
           current_word != "FIELD" &&
           current_word != "POINT_DATA" &&
           current_word != "VECTORS" &&
           !m_stream->eof()) {
        *m_stream >> current_word;

        if (current_word == "SCALARS") {
            std::string scalar_name, scalar_type, scalar_nb;
            *m_stream >> scalar_name >> scalar_type >> scalar_nb;
            *m_stream >> current_word; // For LOOKUP_TABLE
            *m_stream >> current_word; // For default
            std::cout<<"scalar : "<<scalar_name<<" "<<scalar_type<<" "<<scalar_nb<<std::endl;
            if (scalar_name != "GMDS_ID"){
                if (scalar_type == "int"){
                    std::vector<int> tmp_val;
                    bool only_null = true;
                    tmp_val.reserve(nb_values);
                    int nb_true_cells = 0;
                    for (int i = 0; i < nb_values; i++){
                        int val;
                        *m_stream >> val;
                        if (m_cell_types[i]==GMDS_TETRA ||
                            m_cell_types[i]==GMDS_HEX ||
                            m_cell_types[i]==GMDS_PYRAMID ||
                            m_cell_types[i]==GMDS_PRISM3 ) {
                            tmp_val.push_back(val);
                            nb_true_cells++;
                            if(val!=0){
                                only_null=false;
                            }
                        }
                    }
                    //if only_node is true, we expect that this cell data is
                    // on a cell type that is not edge.
                    if(!only_null) {
                        IMeshIOService::DataInt d;
                        d.name = scalar_name;
                        for(int i_cell = 0; i_cell<nb_true_cells;i_cell++){
                            d.values[i_cell]=tmp_val[i_cell];
                        }
                        m_mesh_service->addDataIntRegions(d);
                    }
                }
                else if (scalar_type == "float"){
                    std::vector<double> tmp_val;
                    bool only_null = true;
                    tmp_val.reserve(nb_values);
                    int nb_true_cells = 0;
                    for (int i = 0; i < nb_values; i++){
                        double val;
                        *m_stream >> val;
                        if (m_cell_types[i]==GMDS_TETRA ||
                            m_cell_types[i]==GMDS_HEX ||
                            m_cell_types[i]==GMDS_PYRAMID ||
                            m_cell_types[i]==GMDS_PRISM3 ) {
                            tmp_val.push_back(val);
                            nb_true_cells++;
                            if(val!=0.0){
                                only_null=false;
                            }
                        }
                    }
                    //if only_node is true, we expect that this cell data is
                    // on a cell type that is not edge.
                    if(!only_null) {
                        IMeshIOService::DataReal d;
                        d.name = scalar_name;
                        for(int i_cell = 0; i_cell<nb_true_cells;i_cell++){
                            d.values[i_cell]=tmp_val[i_cell];
                        }
                        m_mesh_service->addDataRealRegions(d);
                    }
                }
            }
        }
        else if (current_word == "VECTORS")
        {
            std::string vector_name, vector_type;
            *m_stream >> vector_name >> vector_type;

            std::cout<<"vector : "<<vector_name<<" "<<vector_type<<std::endl;

            std::vector<math::Vector3d> tmp_val;
            bool only_null = true;
            tmp_val.reserve(nb_values);
            int nb_true_cells = 0;
            for (int i = 0; i < nb_values; i++){
                double x=0, y=0, z=0;
                *m_stream >> x;
                *m_stream >> y;
                *m_stream >> z;
                if (m_cell_types[i]==GMDS_TETRA ||
                    m_cell_types[i]==GMDS_HEX ||
                    m_cell_types[i]==GMDS_PYRAMID ||
                    m_cell_types[i]==GMDS_PRISM3 ) {
                    tmp_val.push_back(math::Vector3d({x,y,z}));
                    nb_true_cells++;
                    if(!(x==0.0 && y==0 && z==0)){
                        only_null=false;
                    }
                }
            }
            //if only_node is true, we expect that this cell data is
            // on a cell type that is not edge.
            if(!only_null) {
                IMeshIOService::DataVector d;
                d.name = vector_name;
                for(int i_cell = 0; i_cell<nb_true_cells;i_cell++){
                    d.values[i_cell]=tmp_val[i_cell];
                }
                m_mesh_service->addDataVectorRegions(d);
            }
        }
    }
}
