/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gmds/io/LimaWriter.h>
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
#include <Lima/malipp.h>
#include <Lima/erreur.h>
#include <Lima/polyedre.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
LimaWriter::LimaWriter(Mesh& AMesh)
  : m_mesh(AMesh), m_length_unit(1.), m_writer(nullptr)
{}
/*----------------------------------------------------------------------------*/
LimaWriter::~LimaWriter()= default;
/*----------------------------------------------------------------------------*/
void
LimaWriter::setLengthUnit(double AUnit)
{
	m_length_unit = AUnit;
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::activateZlibCompression()
{
	//WARNING: only available using the hdf145 extension
	//writer_->activer_compression_zlib();
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::write(const std::string& AFileName, MeshModel AModel, int ACompact)
{
	try {
		m_writer = new Lima::MaliPPWriter2(AFileName, 1);

		m_writer->unite_longueur(m_length_unit);
		Lima::dim_t dim;
		if(m_mesh.getDim() == 3) {
			dim = Lima::D3;
		} else if(m_mesh.getDim() == 2) {
			dim = Lima::D2;
		} else {
			dim = Lima::D1;
		}
		m_writer->dimension(dim);
		m_writer->beginWrite();
		
		writeNodes();
		writeEdges();
		writeFaces();
		writeRegions();

		writeClouds();
		writeLines();
		writeSurfaces();
		writeVolumes();

		writeNodesAttributes();
		writeEdgesAttributes();
		writeFacesAttributes();
		writeRegionsAttributes();

		writeCloudsAttributes();
		writeLinesAttributes();
		writeSurfacesAttributes();
		writeVolumesAttributes();

		m_writer->close ( );
	}
	catch(Lima::write_erreur& e) {
		std::cerr<<"LimaWriterAPI::write error: "<<e.what()<<std::endl;
		throw GMDSException(e.what());
	}
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::checkContinuousNodes(){
	bool isContiguous = true;
	TCellID minID = NullID;

	if(m_mesh.getNbNodes() > 0) {
		Mesh::nodes_iterator it_nodes_id = m_mesh.nodes_begin();
		TCellID currentID = *it_nodes_id;
		minID = currentID+1;
		++it_nodes_id;

		for(; it_nodes_id!=m_mesh.nodes_end(); ++it_nodes_id){
			if(*it_nodes_id != currentID+1) {
				isContiguous = false;
				break;
			}
			currentID = *it_nodes_id;
		}
	}
	try {
		m_writer->writeNodesInfo(isContiguous, m_mesh.getNbNodes(),minID);
	}
	catch(Lima::write_erreur& e) {
		std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
		throw GMDSException(e.what());
	}
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::checkContinuousEdges(){
	bool isContiguous = true;
	TCellID minID = NullID;
	if(m_mesh.getNbEdges() > 0) {
		Mesh::edges_iterator it_edges_id = m_mesh.edges_begin();
		TCellID currentID = *it_edges_id;
		minID=currentID+1;
		++it_edges_id;

		for(; it_edges_id !=m_mesh.edges_end(); ++it_edges_id){
			if(*it_edges_id != currentID+1) {
				isContiguous = false;
				break;
			}
			currentID = *it_edges_id;
		}
	}

	try {
		m_writer->writeEdgesInfo(isContiguous, m_mesh.getNbEdges(),minID);
	}
	catch(Lima::write_erreur& e) {
		std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
		throw GMDSException(e.what());
	}
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::checkContinuousFaces(){
	bool isContiguous = true;
	TCellID minID = NullID;
	if(m_mesh.getNbFaces() > 0) {
		Mesh::faces_iterator it_faces_id = m_mesh.faces_begin();
		TCellID currentID = *it_faces_id;
		minID=currentID+1;
		++it_faces_id;

		for(; it_faces_id !=m_mesh.faces_end(); ++it_faces_id){
			if(*it_faces_id != currentID+1) {
				isContiguous = false;
				break;
			}
			currentID = *it_faces_id;
		}
	}

	try {
		m_writer->writeFacesInfo(isContiguous, m_mesh.getNbFaces(),minID);
	}
	catch(Lima::write_erreur& e) {
		std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
		throw GMDSException(e.what());
	}
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::checkContinuousRegions(){
	bool isContiguous = true;
	TCellID minID = NullID;
	if(m_mesh.getNbRegions() > 0) {
		Mesh::regions_iterator it_regions_id = m_mesh.regions_begin();
		TCellID currentID = *it_regions_id;
		minID=currentID+1;
		++it_regions_id;

		for(; it_regions_id !=m_mesh.regions_end(); ++it_regions_id){
			if(*it_regions_id != currentID+1) {
				isContiguous = false;
				break;
			}
			currentID = *it_regions_id;
		}
	}
	try {
		m_writer->writeRegionsInfo(isContiguous, m_mesh.getNbRegions(),minID);
	}
	catch(Lima::write_erreur& e) {
		std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
		throw GMDSException(e.what());
	}

}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeNodes() {
	// check whether the ids are contiguous
	checkContinuousNodes();

	const Lima::id_type LimaWriterAPI_NBNODES_CHUNK = 10000;

	auto* xccords = new double[LimaWriterAPI_NBNODES_CHUNK];
	auto* yccords = new double[LimaWriterAPI_NBNODES_CHUNK];
	auto* zccords = new double[LimaWriterAPI_NBNODES_CHUNK];
	auto* ids = new Lima::id_type[LimaWriterAPI_NBNODES_CHUNK];

	Lima::id_type chunkSize = 0;
	for(auto n_id:m_mesh.nodes()){
		Node n = m_mesh.get<Node>(n_id);
		xccords[chunkSize] = n.X();
		yccords[chunkSize] = n.Y();
		zccords[chunkSize] = n.Z();
		ids[chunkSize] = n.id()+1; // +1 because mli ids begin at 1

		chunkSize++;
		if(chunkSize==LimaWriterAPI_NBNODES_CHUNK) {
			try {
				m_writer->writeNodes(chunkSize,xccords,yccords,zccords,ids);
			}
			catch(Lima::write_erreur& e) {
				std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
				throw GMDSException(e.what());
			}
			chunkSize = 0;
		}
	}

	if(chunkSize>0) {
		try {
			m_writer->writeNodes(chunkSize,xccords,yccords,zccords,ids);
		}
		catch(Lima::write_erreur& e) {
			std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
			throw GMDSException(e.what());
		}
	}

	delete[] xccords;
	delete[] yccords;
	delete[] zccords;
	delete[] ids;
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeEdges() {
	// check whether the ids are contiguous
	checkContinuousEdges();

	const Lima::id_type LimaWriterAPI_NBEDGES_CHUNK = 10000;

	auto* edge2nodeIDs = new Lima::id_type[2*LimaWriterAPI_NBEDGES_CHUNK];
	auto* ids= new Lima::id_type[LimaWriterAPI_NBEDGES_CHUNK];

	Lima::id_type chunkSize = 0;

	for(auto e_id:m_mesh.edges())
	{
		Edge e = m_mesh.get<Edge>(e_id);
		std::vector<TCellID> nodesIDs = e.getAllIDs<Node>();
		edge2nodeIDs[2*chunkSize  ] = nodesIDs[0]+1;
		edge2nodeIDs[2*chunkSize+1] = nodesIDs[1]+1;
		ids[chunkSize] = e.id()+1;

		chunkSize++;
		if(chunkSize==LimaWriterAPI_NBEDGES_CHUNK) {
			try {
				m_writer->writeEdges(chunkSize,edge2nodeIDs,ids);
			}
			catch(Lima::write_erreur& e) {
				std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
				throw GMDSException(e.what());
			}
			chunkSize = 0;
		}
	}

	if(chunkSize>0) {
		try {
			m_writer->writeEdges(chunkSize,edge2nodeIDs,ids);
		}
		catch(Lima::write_erreur& e) {
			std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
			throw GMDSException(e.what());
		}
	}

	delete[] edge2nodeIDs;
	delete[] ids;
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeFaces() {
	// check whether the ids are contiguous
	checkContinuousFaces();

	const Lima::id_type LimaWriterAPI_NBFACES_CHUNK = 10000;
	const Lima::id_type LimaWriterAPI_MAX_NBNODES_PER_FACE = 15; //Lima::MAX_NOEUDS;

	auto* face2nodeIDs = new Lima::id_type[LimaWriterAPI_MAX_NBNODES_PER_FACE*LimaWriterAPI_NBFACES_CHUNK];
	auto* nbNodesPerFace = new Lima::id_type[LimaWriterAPI_NBFACES_CHUNK];
	auto* ids = new Lima::id_type[LimaWriterAPI_NBFACES_CHUNK];

	Lima::id_type chunkSize = 0;
	Lima::id_type currentIndex = 0;

	for(auto f_id:m_mesh.faces()){
		Face f = m_mesh.get<Face>(f_id);
		std::vector<TCellID> nodesIDs = f.getAllIDs<Node>();
		nbNodesPerFace[chunkSize] = nodesIDs.size();

		if(nodesIDs.size() > LimaWriterAPI_MAX_NBNODES_PER_FACE) {
			throw GMDSException("LimaWriterAPI::writeFaces a face has too many nodes (> 15 == Lima::MAX_NOEUDS).");
		}

		for(unsigned int nodesID : nodesIDs) {
			face2nodeIDs[currentIndex] = nodesID+1;
			currentIndex++;
		}

		ids[chunkSize] = f.id()+1;

		chunkSize++;
		if(chunkSize==LimaWriterAPI_NBFACES_CHUNK) {
			try {
				m_writer->writeFaces(chunkSize,face2nodeIDs,nbNodesPerFace,ids);
			}
			catch(Lima::write_erreur& e) {
				std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
				throw GMDSException(e.what());
			}
			chunkSize = 0;
			currentIndex = 0;
		}
	}

	if(chunkSize>0) {
		try {
			m_writer->writeFaces(chunkSize,face2nodeIDs,nbNodesPerFace,ids);
		}
		catch(Lima::write_erreur& e) {
			std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
			throw GMDSException(e.what());
		}
	}

	delete[] face2nodeIDs;
	delete[] nbNodesPerFace;
	delete[] ids;
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeRegions(){
	// check whether the ids are contiguous
	checkContinuousRegions();

	const Lima::id_type LimaWriterAPI_NBREGIONS_CHUNK = 10000;
	const Lima::id_type LimaWriterAPI_MAX_NBNODES_PER_REGION = 15; //Lima::MAX_NOEUDS;

	auto* region2nodeIDs = new Lima::id_type[LimaWriterAPI_MAX_NBNODES_PER_REGION*LimaWriterAPI_NBREGIONS_CHUNK];
	auto* regionTypes = new Lima::Polyedre::PolyedreType[LimaWriterAPI_NBREGIONS_CHUNK];
	auto* ids = new Lima::id_type[LimaWriterAPI_NBREGIONS_CHUNK];

	Lima::id_type chunkSize = 0;
	Lima::id_type currentIndex = 0;

	for(auto r_id:m_mesh.regions()){
		Region r = m_mesh.get<Region>(r_id);
		std::vector<TCellID> nodesIDs = r.getAllIDs<Node>();

		switch(r.type()) {
		case GMDS_TETRA :
			regionTypes[chunkSize] = Lima::Polyedre::TETRAEDRE;
			break;
		case GMDS_PYRAMID :
			regionTypes[chunkSize] = Lima::Polyedre::PYRAMIDE;
			break;
		case GMDS_PRISM3 :
			regionTypes[chunkSize] = Lima::Polyedre::PRISME;
			break;
		case GMDS_HEX :
			regionTypes[chunkSize] = Lima::Polyedre::HEXAEDRE;
			break;
		default:
			throw GMDSException("LimaWriterAPI::writeRegions cell type not handled by Lima.");
		}

		if(nodesIDs.size() > LimaWriterAPI_MAX_NBNODES_PER_REGION) {
			throw GMDSException("LimaWriterAPI::writeRegions a face has too many nodes (> 15 == Lima::MAX_NOEUDS).");
		}

		for(unsigned int nodesID : nodesIDs) {
			region2nodeIDs[currentIndex] = nodesID+1;
			currentIndex++;
		}

		ids[chunkSize] = r.id()+1;

		chunkSize++;
		if(chunkSize==LimaWriterAPI_NBREGIONS_CHUNK) {
			try {
				m_writer->writeRegions(chunkSize,region2nodeIDs,regionTypes,ids);
			}
			catch(Lima::write_erreur& e) {
				std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
				throw GMDSException(e.what());
			}
			chunkSize = 0;
			currentIndex = 0;
		}
	}

	if(chunkSize>0) {
		try {
			m_writer->writeRegions(chunkSize,region2nodeIDs,regionTypes,ids);
		}
		catch(Lima::write_erreur& e) {
			std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
			throw GMDSException(e.what());
		}
	}

	delete[] region2nodeIDs;
	delete[] regionTypes;
	delete[] ids;
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeClouds() {
	const Lima::id_type LimaWriterAPI_NBNODES_CHUNK = 10000;

	std::vector<std::string> names;
	std::vector<Lima::id_type> sizes;

	for(auto i=0; i< m_mesh.getNbGroups<Node>(); i++) {
		CellGroup<Node>* cl = m_mesh.getGroup<Node>(i);
		names.push_back(cl->name());
		sizes.push_back(cl->size());
	}

	Lima::id_type nbClouds = m_mesh.getNbGroups<Node>();

	try {
		m_writer->writeNodeSetInfo(nbClouds,names,sizes);
	}
	catch(Lima::write_erreur& e) {
		std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
		throw GMDSException(e.what());
	}

	Lima::id_type ids[LimaWriterAPI_NBNODES_CHUNK];
	Lima::id_type chunkSize = 0;

	for(auto i=0; i< m_mesh.getNbGroups<Node>(); i++) {
		CellGroup<Node>* cl = m_mesh.getGroup<Node>(i);

		std::vector<TCellID> nodeIDs= cl->cells();

		for(unsigned int nodeID : nodeIDs) {
			ids[chunkSize] = nodeID+1;
			chunkSize++;

			if(chunkSize==LimaWriterAPI_NBNODES_CHUNK) {
				try {
					m_writer->writeNodeSetData(cl->name(),chunkSize,ids);
				}
				catch(Lima::write_erreur& e) {
					std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
					throw GMDSException(e.what());
				}
				chunkSize = 0;
			}
		}

		if(chunkSize > 0) {
			try {
				m_writer->writeNodeSetData(cl->name(),chunkSize,ids);
			}
			catch(Lima::write_erreur& e) {
				std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
				throw GMDSException(e.what());
			}
			chunkSize = 0;
		}

	}
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeLines() {
	const Lima::id_type LimaWriterAPI_NBEDGES_CHUNK = 10000;
	std::vector<std::string> names;
	std::vector<Lima::id_type> sizes;

	for(auto i=0; i< m_mesh.getNbGroups<Edge>(); i++) {
		CellGroup<Edge>* l = m_mesh.getGroup<Edge>(i);
		names.push_back(l->name());
		sizes.push_back(l->size());
	}

	Lima::id_type nbLines = m_mesh.getNbGroups<Edge>();

	try {
		m_writer->writeEdgeSetInfo(nbLines,names,sizes);
	}
	catch(Lima::write_erreur& e) {
		std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
		throw GMDSException(e.what());
	}

	Lima::id_type ids[LimaWriterAPI_NBEDGES_CHUNK];
	Lima::id_type chunkSize = 0;

	for(auto i=0; i< m_mesh.getNbGroups<Edge>(); i++) {
		CellGroup<Edge>* l = m_mesh.getGroup<Edge>(i);

		std::vector<TCellID> edgeIDs= l->cells();

		for(unsigned int edgeID : edgeIDs) {
			ids[chunkSize] = edgeID+1;
			chunkSize++;

			if(chunkSize==LimaWriterAPI_NBEDGES_CHUNK) {
				try {
					m_writer->writeEdgeSetData(l->name(),chunkSize,ids);
				}
				catch(Lima::write_erreur& e) {
					std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
					throw GMDSException(e.what());
				}
				chunkSize = 0;
			}
		}

		if(chunkSize > 0) {
			try {
				m_writer->writeEdgeSetData(l->name(),chunkSize,ids);
			}
			catch(Lima::write_erreur& e) {
				std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
				throw GMDSException(e.what());
			}
			chunkSize = 0;
		}
	}
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeSurfaces() {
	const Lima::id_type LimaWriterAPI_NBFACES_CHUNK = 10000;
	std::vector<std::string> names;
	std::vector<Lima::id_type> sizes;

	for(auto i=0; i< m_mesh.getNbGroups<Face>(); i++) {
		CellGroup<Face>* surf = m_mesh.getGroup<Face>(i);
		names.push_back(surf->name());
		sizes.push_back(surf->size());
	}

	Lima::id_type nbSurfaces = m_mesh.getNbGroups<Face>();

	try {
		m_writer->writeFaceSetInfo(nbSurfaces,names,sizes);
	}
	catch(Lima::write_erreur& e) {
		std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
		throw GMDSException(e.what());
	}

	Lima::id_type ids[LimaWriterAPI_NBFACES_CHUNK];
	Lima::id_type nbNodes[LimaWriterAPI_NBFACES_CHUNK];
	Lima::id_type chunkSize = 0;

	for(auto i=0; i< m_mesh.getNbGroups<Face>(); i++) {
		CellGroup<Face>* surf = m_mesh.getGroup<Face>(i);

		std::vector<TCellID> faceIDs= surf->cells();

		for(unsigned int & faceID : faceIDs) {
			ids[chunkSize] = faceID+1;
			nbNodes[chunkSize] = (m_mesh.get<Face> (faceID)).nbNodes();
			chunkSize++;

			if(chunkSize==LimaWriterAPI_NBFACES_CHUNK) {
				try {
					m_writer->writeFaceSetData(surf->name(),chunkSize,ids,nbNodes);
				}
				catch(Lima::write_erreur& e) {
					std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
					throw GMDSException(e.what());
				}
				chunkSize = 0;
			}
		}

		if(chunkSize > 0) {
			try {
				m_writer->writeFaceSetData(surf->name(),chunkSize,ids,nbNodes);
			}
			catch(Lima::write_erreur& e) {
				std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
				throw GMDSException(e.what());
			}
			chunkSize = 0;
		}

	}
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeVolumes() {
	const Lima::id_type LimaWriterAPI_NBREGIONS_CHUNK = 10000;
	std::vector<std::string> names;
	std::vector<Lima::id_type> sizes;

	for(auto i=0; i< m_mesh.getNbGroups<Region>(); i++) {
		CellGroup<Region>* vol = m_mesh.getGroup<Region>(i);
		names.push_back(vol->name());
		sizes.push_back(vol->size());
	}

	Lima::id_type nbVolumes = m_mesh.getNbGroups<Region>();

	try {
		m_writer->writeRegionSetInfo(nbVolumes,names,sizes);
	}
	catch(Lima::write_erreur& e) {
		std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
		throw GMDSException(e.what());
	}

	Lima::id_type ids[LimaWriterAPI_NBREGIONS_CHUNK];
	Lima::Polyedre::PolyedreType types[LimaWriterAPI_NBREGIONS_CHUNK];
	Lima::id_type chunkSize = 0;

	for(auto i=0; i< m_mesh.getNbGroups<Region>(); i++) {
		CellGroup<Region>* vol = m_mesh.getGroup<Region>(i);

		std::vector<TCellID> regionIDs= vol->cells();

		for(unsigned int & regionID : regionIDs) {
			ids[chunkSize] = regionID+1;
			switch((m_mesh.get<Region>(regionID)).type()) {
			case GMDS_TETRA :
				types[chunkSize] = Lima::Polyedre::TETRAEDRE;
				break;
			case GMDS_PYRAMID :
				types[chunkSize] = Lima::Polyedre::PYRAMIDE;
				break;
			case GMDS_PRISM3 :
				types[chunkSize] = Lima::Polyedre::PRISME;
				break;
			case GMDS_HEX :
				types[chunkSize] = Lima::Polyedre::HEXAEDRE;
				break;
			default:
				throw GMDSException("LimaWriterAPI::writeRegions cell type not handled by Lima.");
			}
			chunkSize++;

			if(chunkSize==LimaWriterAPI_NBREGIONS_CHUNK) {
				try {
					m_writer->writeRegionSetData(vol->name(),chunkSize,ids,types);
				}
				catch(Lima::write_erreur& e) {
					std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
					throw GMDSException(e.what());
				}
				chunkSize = 0;
			}
		}

		if(chunkSize > 0) {
			try {
				m_writer->writeRegionSetData(vol->name(),chunkSize,ids,types);
			}
			catch(Lima::write_erreur& e) {
				std::cerr<<"LimaWriterAPI error: "<<e.what()<<std::endl;
				throw GMDSException(e.what());
			}
			chunkSize = 0;
		}
	}
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeNodesAttributes() {
	m_writer->writeNodeAttributes();
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeEdgesAttributes() {
	m_writer->writeEdgeAttributes();
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeFacesAttributes(){
	m_writer->writeFaceAttributes();
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeRegionsAttributes(){
	m_writer->writeRegionAttributes();
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeCloudsAttributes(){
	m_writer->writeNodeSetsAttributes();
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeLinesAttributes(){
	m_writer->writeEdgeSetsAttributes();
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeSurfacesAttributes(){
	m_writer->writeFaceSetsAttributes();
}
/*----------------------------------------------------------------------------*/
void
LimaWriter::writeVolumesAttributes(){
	m_writer->writeRegionSetsAttributes();
}
/*----------------------------------------------------------------------------*/
