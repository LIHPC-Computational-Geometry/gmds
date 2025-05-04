
#include <fstream>
#include <gmds/io/GMSHWriter.h>
#include <gmds/utils/CommonTypes.h>

namespace gmds {
namespace {

size_t
getElementTypeToGmshTypeMap(ECellType type)
{
	switch (type) {
	case ECellType::GMDS_TRIANGLE:
	case ECellType::GMDS_FACE: return 2;
	case ECellType::GMDS_QUAD: return 3;
	default: throw GMDSException("cell type not supported by the gmsh writer");
	}
}

}     // namespace

GMSHWriter::GMSHWriter(IMeshIOService *AMeshService) : IWriter(AMeshService) {}

void
GMSHWriter::initialize(const std::string &AFileName)
{
	m_stream = new std::ofstream(AFileName, std::ios::out);
	if (!m_stream) {
		std::string s = "Impossible to create a .msh File : " + AFileName;
		throw GMDSException(s);
	}

	*m_stream << "$MeshFormat\n";
	*m_stream << "2.2 0 8\n";
	*m_stream << "$EndMeshFormat\n";
}

void
GMSHWriter::writeNodes()
{
	std::vector<IMeshIOService::NodeInfo> nodesInfo;
	m_mesh_service->getNodes(nodesInfo);

	*m_stream << "$Nodes\n";
	*m_stream << nodesInfo.size() << "\n";
	
	(*m_stream).precision(15);
	for (const auto& info : nodesInfo) {
		const math::Point p = info.point;
		*m_stream << info.id << " " << p.X() << " " << p.Y() << " " << p.Z() << "\n";
	}
	(*m_stream).precision(6);

	*m_stream << "$EndNodes\n";
}

void
GMSHWriter::writeFaces()
{
	std::vector<IMeshIOService::CellInfo> facesInfo;
	m_mesh_service->getFaces(facesInfo);

	*m_stream << "$Elements\n";
	*m_stream << facesInfo.size() << "\n";
	for (const auto& info : facesInfo) {

		*m_stream << info.id << " ";
		*m_stream << getElementTypeToGmshTypeMap(info.type) << " ";
		*m_stream << "2 0 0";     // 0 tag for now (no physical entity, 0 geom id)

		for (const size_t id : info.node_ids) {
			*m_stream << " " << id;
		}
		*m_stream << "\n";
	}
	*m_stream << "$EndElements\n";
}

void
GMSHWriter::finalize()
{
	m_stream->close();
}

}     // namespace gmds
