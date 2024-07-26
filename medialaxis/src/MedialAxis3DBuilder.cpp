//
// Created by chenyt on 25/04/24.
//
/*----------------------------------------------------------------------------*/
#include "gmds/medialaxis/MedialAxis3DBuilder.h"
#include "gmds/ig/MeshDoctor.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include <gmds/math/Point.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace medialaxis {
/*----------------------------------------------------------------------------*/
MedialAxis3DBuilder::MedialAxis3DBuilder(Mesh &AMesh){
	m_mesh = &AMesh;
	// Correspondence between tetras and their groups
	m_mesh->newVariable<int,GMDS_REGION>("Tetra2Group");
	// Mark with 1 nearly flat tetras
	m_mesh->newVariable<int,GMDS_REGION>("is_nearly_flat");
	// Delaunay tetra / medial point correspondence
	m_mesh->newVariable<int,GMDS_REGION>("Tetra2MP");
	// Mark with 1 the problematic tetras (those that have a Steiner point for vertex)
	m_mesh->newVariable<int,GMDS_REGION>("isProblematic");
	// Points classification (1 if the point is considered a detail, 0 if not)
	m_mesh->newVariable<int,GMDS_NODE>("is_detail");
	// Faces classification (1 if the face is considered a detail, 0 if not)
	m_mesh->newVariable<int,GMDS_FACE>("face_is_detail");
	m_ma = new MedialAxis3D();
}
/*----------------------------------------------------------------------------*/
MedialAxis3DBuilder::~MedialAxis3DBuilder()
{
	if(m_ma!= nullptr)
		delete m_ma;
}
/*----------------------------------------------------------------------------*/
MedialAxis3D* MedialAxis3DBuilder::getMedialObject()
{
	return m_ma;
}

/*-------------------------------------------------------------------------*/
Mesh* MedialAxis3DBuilder::getMesh()
{
	return m_mesh;
}

/*-------------------------------------------------------------------------*/
bool MedialAxis3DBuilder::isAnInternalFace(const TCellID &AFaceID)
{
	Face f = m_mesh->get<Face>(AFaceID);
	std::vector<Region> rf = f.get<Region>();
	if (rf.size()==1)
		return false;
	return true;
}

/*-------------------------------------------------------------------------*/
bool MedialAxis3DBuilder::isAnInternalEdge(const TCellID &AEdgeID)
{
	Edge e = m_mesh->get<Edge>(AEdgeID);
	std::vector<Face> ef = e.get<Face>();
	for (auto f:ef)
	{
		if (!isAnInternalFace(f.id()))
			return false;
	}
	return true;
}

/*-------------------------------------------------------------------------*/
bool MedialAxis3DBuilder::shareFace(const gmds::TCellID &ATetraID1, const gmds::TCellID &ATetraID2)
{
	Region Tetra1 = m_mesh->get<Region>(ATetraID1);
	Region Tetra2 = m_mesh->get<Region>(ATetraID2);
	std::vector<Face> faces1 = Tetra1.get<Face>();
	std::vector<Face> faces2 = Tetra2.get<Face>();
	for (auto f1:faces1)
	{
		for (auto f2:faces2)
		{
			if (f1.id() == f2.id())
				return true;
		}
	}
	return false;
}

/*-------------------------------------------------------------------------*/
bool MedialAxis3DBuilder::isDelaunay()
{
	double min_diff = 0;
	TCellID prob_point;
	TCellID prob_tetra;
	for (auto t_id:m_mesh->regions())
	{
		Region Tetra = m_mesh->get<Region>(t_id);
		std::vector<Node> nT = Tetra.get<Node>();
		math::Tetrahedron T(nT[0].point(), nT[1].point(), nT[2].point(), nT[3].point());
		math::Point Center = circumcenterTetra(T);
		double Radius = nT[0].point().distance(Center);
		for (auto n_id:m_mesh->nodes())
		{
			double dist = m_mesh->get<Node>(n_id).point().distance(Center);
			if (dist - Radius < min_diff)
			{
				min_diff = dist - Radius;
				prob_point = n_id;
				prob_tetra = t_id;
			}
		}
	}
	if (min_diff >= 0)
	{
		std::cout<<"The input mesh is of Delaunay"<<std::endl;
		return true;
	}
	std::cout<<"The node "<<prob_point<<" of the input tetrahedralization is too close to the circumcenter of the tetra "<<prob_tetra<<std::endl;
	std::cout<<"The difference between the distance node/circumcenter and the radius should be >= 0 but is "<<min_diff<<std::endl;
	std::cout<<"The input mesh is NOT of Delaunay"<<std::endl;
	return false;
}

/*-------------------------------------------------------------------------*/
std::vector<Node> MedialAxis3DBuilder::minDelProblematicPoints()
{
	std::vector<Node> probNodes;
	for (auto n_id:m_mesh->nodes())
	{
		bool is_problematic = true;
		Node n = m_mesh->get<Node>(n_id);
		std::vector<Face> Faces = n.get<Face>();
		for(auto face:Faces)
		{
			if (face.get<Region>().size() == 1)
			{
				is_problematic = false;
				break;
			}
		}
		if (is_problematic)
			probNodes.push_back(n);
	}
	return probNodes;
}

/*-------------------------------------------------------------------------*/
void MedialAxis3DBuilder::removeProblematicNodes()
{
	std::cout<<"Removing problematic cells...";
	std::vector<Node> probNodes = minDelProblematicPoints();
	for (auto n:probNodes)
	{
		for (auto tetra:n.get<Region>())
		{
			for (auto face:tetra.get<Face>())
				m_mesh->deleteFace(face);
			for (auto edge:tetra.get<Edge>())
				m_mesh->deleteEdge(edge);
			m_mesh->deleteRegion(tetra);
		}
		m_mesh->deleteNode(n);
	}
	// Updating connectivity
	MeshDoctor doc(m_mesh);
	doc.buildEdgesAndX2E();
	doc.buildFacesAndR2F();
	doc.updateUpwardConnectivity();
	std::cout<<" done."<<std::endl;
}

/*-------------------------------------------------------------------------*/
void MedialAxis3DBuilder::markNearlyFlatTetras()
{
	auto IsNearlyFlat = m_mesh->getVariable<int,GMDS_REGION>("is_nearly_flat");
	for (auto t_id:m_mesh->regions())
	{
		Region Tetra = m_mesh->get<Region>(t_id);
		int NbInternalEdges = 0;
		for (auto edge:Tetra.get<Edge>())
		{
			if (isAnInternalEdge(edge.id()))
				NbInternalEdges += 1;
		}
		if (NbInternalEdges == 1)
			IsNearlyFlat->set(Tetra.id(),1);
	}
}

/*-------------------------------------------------------------------------*/
void MedialAxis3DBuilder::dealWithNearlyFlatTetras()
{
	// Requires Tetra2MP, IsNearlyFlat
	auto Tetra2MP = m_mesh->getVariable<int,GMDS_REGION>("Tetra2MP");
	auto IsNearlyFlat = m_mesh->getVariable<int,GMDS_REGION>("is_nearly_flat");
	for (auto tetra_id:m_mesh->regions())
	{
		Region tetra = m_mesh->get<Region>(tetra_id);
		if (IsNearlyFlat->value(tetra.id()) == 1)
		{
			// Find a neighbouring internal tetra
			Region reliableNeighbour;
			bool reliableNeighbourFound = false;
			for (auto face:tetra.get<Face>())
			{
				for (auto tetra2:face.get<Region>())
				{
					if (tetra2.id() != tetra.id())
					{
						reliableNeighbour = tetra2;
						reliableNeighbourFound = true;
						break;
					}
				}
				if (reliableNeighbourFound)
					break;
			}
			TCellID mp1 = Tetra2MP->value(tetra.id());
			TCellID mp2 = Tetra2MP->value(reliableNeighbour.id());
			std::vector<TCellID> neighbour;
			neighbour.push_back(mp2);
			m_ma->averagePoint(mp1,neighbour);
		}
	}
}

/*-------------------------------------------------------------------------*/
double MedialAxis3DBuilder::meshStep()
{
	double h = 0;
	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		if (!isAnInternalEdge(e.id()))
		{
			if (h < e.length())
				h = e.length();
		}
	}
	return h;
}

/*-------------------------------------------------------------------------*/
std::vector<std::vector<TCellID>>
MedialAxis3DBuilder::tetraFilter(const double &ATol)
{
	std::cout<<"> Filtering the tetras"<<std::endl;
	auto Tetra2Group = m_mesh->getVariable<int,GMDS_REGION>("Tetra2Group");
	std::vector<std::vector<TCellID>> tetraGroups;
	int groupID = 0;
	std::vector<int> tetras_to_test(m_mesh->getNbTetrahedra());
	for (int n_id = 0; n_id < m_mesh->getNbTetrahedra(); n_id++)
		tetras_to_test[n_id] = 1;
	for (int n_id = 0; n_id < m_mesh->getNbTetrahedra(); n_id++)
	{
		if (tetras_to_test[n_id] == 1)
		{
			std::vector<TCellID> newGroup;
			newGroup.push_back(n_id);
			Tetra2Group->set(n_id,groupID);
			tetras_to_test[n_id] = 0;
			// Circumcenter
			Region Tetra = m_mesh->get<Region>(n_id);
			std::vector<Node> nT = Tetra.get<Node>();
			math::Tetrahedron T(nT[0].point(), nT[1].point(), nT[2].point(), nT[3].point());
			math::Point Ctr1 = circumcenterTetra(T);
			for (int m_id = n_id + 1; m_id < m_mesh->getNbTetrahedra(); m_id++)
			{
				if (tetras_to_test[m_id] == 1)
				{
					Region Tetra2 = m_mesh->get<Region>(m_id);
					std::vector<Node> nT2 = Tetra2.get<Node>();
					math::Tetrahedron T2(nT2[0].point(), nT2[1].point(), nT2[2].point(), nT2[3].point());
					math::Point Ctr2 = circumcenterTetra(T2);
					double dist = Ctr1.distance(Ctr2);
					if (dist <= ATol)
					{
						newGroup.push_back(m_id);
						Tetra2Group->set(m_id,groupID);
						tetras_to_test[m_id] = 0;
					}
				}
			}
			tetraGroups.push_back(newGroup);
			groupID+=1;
		}
	}
	return tetraGroups;
}

/*----------------------------------------------------------------------------*/
void MedialAxis3DBuilder::buildMedialPoints()
{
	// Delaunay tetra / medial point correspondence
	auto Tetra2MP = m_mesh->getVariable<int,GMDS_REGION>("Tetra2MP");

	// Marking the problematic tetras
	std::vector<Node> probNodes = minDelProblematicPoints();
	auto isProblematic = m_mesh->getVariable<int,GMDS_REGION>("isProblematic");
	for (auto t_id:m_mesh->regions())
	{
		isProblematic->set(t_id,0);
	}
	for (auto n:probNodes)
	{
		for (auto tetra:n.get<Region>())
		{
			isProblematic->set(tetra.id(),1);
		}
	}

	// Adding the medial points
	std::cout<<"> Creating the medial points"<<std::endl;
	// Medial points of the regular tetras
	for (auto n_id:m_mesh->regions()){
		if (isProblematic->value(n_id) == 0)
		{
			Region Tetra = m_mesh->get<Region>(n_id);
			std::vector<Node> nT = Tetra.get<Node>();
			math::Tetrahedron T(nT[0].point(), nT[1].point(), nT[2].point(), nT[3].point());
			math::Point MP = circumcenterTetra(T);
			Node n = m_ma->newMedPoint(MP);
			Tetra2MP->set(Tetra.id(), n.id());
		}
	}
	// Medial points of the problematic tetras
	for (auto n: minDelProblematicPoints())
	{
		for (auto tetra:n.get<Region>())
		{
			Node mp = m_ma->newMedPoint(n.point());
			Tetra2MP->set(tetra.id(),mp.id());
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis3DBuilder::buildMedialEdges()
{
	// Delaunay tetra / medial point correspondence
	auto Tetra2MP = m_mesh->getVariable<int,GMDS_REGION>("Tetra2MP");
	// Create the medial edges
	std::cout<<"> Creating the medial edges"<<std::endl;
	for (auto n_id:m_mesh->regions()) {
		Region Tetra = m_mesh->get<Region>(n_id);
		int mp_id1 = Tetra2MP->value(Tetra.id());
		std::vector<Region> adj_tetras;
		for (auto f : Tetra.get<Face>()) {
			f.get<Region>(adj_tetras);
			for (auto tetra : adj_tetras) {
				if (tetra.id() < n_id) {
					int mp_id2 = Tetra2MP->value(tetra.id());
					Edge me = m_ma->newMedEdge(mp_id1, mp_id2);
				}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis3DBuilder::buildMedialFaces()
{
	// Delaunay tetra / medial point correspondence
	auto Tetra2MP = m_mesh->getVariable<int,GMDS_REGION>("Tetra2MP");
	// Create the medial faces
	std::cout<<"> Creating the medial faces"<<std::endl;
	for (auto n_id:m_mesh->edges())
	{
		if (isAnInternalEdge(n_id)) {
			Edge e = m_mesh->get<Edge>(n_id);
			std::vector<Region> adj_tetras = e.get<Region>();
			// Sorting the tetras
			std::vector<TCellID> adj_tetras_sorted(adj_tetras.size());
			std::vector<int> tetras_to_test(adj_tetras.size());
			for (int i = 1; i < adj_tetras.size(); i++) {
				tetras_to_test[i] = adj_tetras[i].id();
			}
			tetras_to_test[0] = -1;
			adj_tetras_sorted[0] = adj_tetras[0].id();
			for (int i = 0; i < adj_tetras.size() - 1; i++) {
				for (int j = 1; j < adj_tetras.size(); j++) {
					if (tetras_to_test[j] >= 0) {
						if (shareFace(adj_tetras_sorted[i], tetras_to_test[j])) {
							adj_tetras_sorted[i + 1] = tetras_to_test[j];
							tetras_to_test[j] = -1;
							break;
						}
					}
				}
			}
			// Tetras sorted, creating the medial faces
			std::vector<TCellID> mp(adj_tetras.size());
			for (int i = 0; i < adj_tetras.size(); i++) {
				mp[i] = Tetra2MP->value(adj_tetras_sorted[i]);
			}
			Face mf = m_ma->newMedFace(mp);
		}

	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis3DBuilder::repositionProblematicMedialPoints()
{
	// Delaunay tetra / medial point correspondence
	auto Tetra2MP = m_mesh->getVariable<int,GMDS_REGION>("Tetra2MP");

	// Problematic (Steiner) nodes and tetras
	std::vector<Node> probNodes = minDelProblematicPoints();
	auto isProblematic = m_mesh->getVariable<int,GMDS_REGION>("isProblematic");

	// Fixing the position of the problematic medial points
	for (auto n:probNodes)
	{
		std::vector<TCellID> regularMedPoint;
		for (auto tetra:n.get<Region>())
		{
			for (auto face:tetra.get<Face>())
			{
				for (auto tetra2:face.get<Region>())
				{
					if (isProblematic->value(tetra2.id()) == 0)
					{
						bool alreadySeen = false;
						for (auto id:regularMedPoint)
						{
							if (Tetra2MP->value(tetra2.id()) == id)
								alreadySeen = true;
						}
						if (!alreadySeen)
							regularMedPoint.push_back(Tetra2MP->value(tetra2.id()));
					}
				}
			}
		}
		for (auto tetra:n.get<Region>())
		{
			m_ma->averagePoint(Tetra2MP->value(tetra.id()), regularMedPoint);
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis3DBuilder::setMedialRadius()
{
	// Requires the medial axis already constructed
	auto Tetra2MP = m_mesh->getVariable<int,GMDS_REGION>("Tetra2MP");
	for (auto t_id:m_mesh->regions())
	{
		Region Tetra = m_mesh->get<Region>(t_id);
		std::vector<Node> adj_nodes = Tetra.get<Node>();
		math::Point MP = m_ma->getMedPoint(Tetra2MP->value(Tetra.id())).point();
		double r = MP.distance(adj_nodes[0].point());
		m_ma->setMedialRadius(Tetra2MP->value(Tetra.id()),r);
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis3DBuilder::markGeometryDetailsOnPoints()
{
	auto tetra2MP = m_mesh->getVariable<int,GMDS_REGION>("Tetra2MP");
	auto isDetail = m_mesh->getVariable<int,GMDS_NODE>("is_detail");
	for (auto n_id:m_mesh->nodes())
		isDetail->set(n_id,0);
	for (auto t_id:m_mesh->regions())
	{
		TCellID medPoint = tetra2MP->value(t_id);
		if (m_ma->correspondsToDetail(medPoint) == 1)
		{
			Region tetra = m_mesh->get<Region>(t_id);
			for (auto n:tetra.get<Node>())
				isDetail->set(n.id(),1);
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis3DBuilder::markGeometryDetailsOnFaces()
{
	auto isDetail = m_mesh->getVariable<int,GMDS_NODE>("is_detail");
	auto faceIsDetail = m_mesh->getVariable<int,GMDS_FACE>("face_is_detail");
	for (auto f_id:m_mesh->faces())
	{
		faceIsDetail->set(f_id,0);
		Face f = m_mesh->get<Face>(f_id);
		bool detail = true;
		for (auto n:f.get<Node>())
		{
			if (isDetail->value(n.id()) == 0)
			{
				detail = false;
				break;
			}
		}
		if (detail && (f.get<Region>().size() == 1))
			faceIsDetail->set(f_id,1);
	}

}

/*----------------------------------------------------------------------------*/
double MedialAxis3DBuilder::maxDelEdgeLength()
{
	double maxLength = 0.;
	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		if (e.length() > maxLength)
			maxLength = e.length();
	}
	return maxLength;
}

/*----------------------------------------------------------------------------*/
void MedialAxis3DBuilder::identifyDetails()
{
	std::cout<<" "<<std::endl;
	std::cout<<"== Identifying details on the geometry using the medial axis =="<<std::endl;

	// Identify medial surfaces corresponding to details
	m_ma->identifyDetailsSurfaces(maxDelEdgeLength()/20.);
	//std::cout<<"Max Del edge length : "<<maxDelEdgeLength()<<std::endl;

	// Mark geometry details on points
	markGeometryDetailsOnPoints();

	// Mark geometry details on faces
	markGeometryDetailsOnFaces();

	std::cout<<"==============================================================="<<std::endl;
	std::cout<<" "<<std::endl;
}

/*----------------------------------------------------------------------------*/
MedialAxis3DBuilder::STATUS

MedialAxis3DBuilder::execute()
{
//	MeshDoctor doc(m_mesh);
//	doc.buildEdgesAndX2E();
//	doc.buildFacesAndR2F();
//	doc.updateUpwardConnectivity();
	// Check if the tetrahedralization is of Delaunay
	//isDelaunay();
	// Look for problematic (internal) nodes in m_mesh
	//std::cout<<"There is "<<minDelProblematicPoints().size()<<" problematic (ie internal) nodes in the input tetrahedralization."<<std::endl;
	// Removing the problematic cells
	//removeProblematicNodes();

	// Build the medial points
	buildMedialPoints();

	// Build medial edges
	buildMedialEdges();

	// Build the medial faces
	buildMedialFaces();

	// Reposition medial points coming from Steiner tetras
	repositionProblematicMedialPoints();


	// Setting connectivity of the medial axis
	m_ma->updateConnectivity();

	// Set medial edges type
	m_ma->setMedialEdgeType();

	// Set medial points type
	m_ma->setMedialPointType();

	// Set medial faces type
	m_ma->setMedialFaceType();

	// Set medial surface IDs on faces
	m_ma->setMedialSurfaceIdOnFaces();

	// Set medial surface IDs on points
	m_ma->setMedialSurfaceIdOnPoints();

	// Set medial surface type on faces
	m_ma->setMedialSurfaceTypeOnFaces();

	// Set medial surface size on faces
	m_ma->setMedialSurfaceSizeOnFaces();

	// Marking the nearly flat tetras
	markNearlyFlatTetras();

	// Deal with the nearly flat tetras
	dealWithNearlyFlatTetras();

	// Flatten unwanted faces
	//m_voronoi_medax->flattenUnwantedFaces();

	// Set medial radii
	setMedialRadius();

	// Identify details on geometry using the (Voronoï) medial axis
	identifyDetails();



	return MedialAxis3DBuilder::SUCCESS;
}

/*----------------------------------------------------------------------------*/
MedialAxis3DBuilder::STATUS

MedialAxis3DBuilder::executeFilter()
{
	// Correspondence between tetras and medial points
	auto Tetra2MP = m_mesh->getVariable<int,GMDS_REGION>("Tetra2MP");
	// Create the medial points
	double h = meshStep();
	std::cout<<"Mesh step: "<<h<<std::endl;
	std::vector<std::vector<TCellID>> tetraGroups = tetraFilter(-1);
	std::cout<<"> Creating the medial points"<<std::endl;
	for (int g_id = 0; g_id < tetraGroups.size(); g_id++)
	{
		std::vector<TCellID> tetraGroup = tetraGroups[g_id];
		int n_id = tetraGroup[0];
		Region Tetra = m_mesh->get<Region>(n_id);
		std::vector<Node> nT = Tetra.get<Node>();
		math::Tetrahedron T(nT[0].point(), nT[1].point(), nT[2].point(), nT[3].point());
		math::Point MP = circumcenterTetra(T);
		Node n = m_ma->newMedPoint(MP);
		for (auto id:tetraGroup)
		{
			Tetra2MP->set(id,n.id());
		}
	}
	// Create the medial edges
	std::cout<<"> Creating the medial edges"<<std::endl;
	auto Tetra2Group = m_mesh->getVariable<int,GMDS_REGION>("Tetra2Group");
	for (int g_id = 0; g_id < tetraGroups.size(); g_id++)
	{
		std::vector<TCellID> tetraGroup = tetraGroups[g_id];
		int n_id;
		int groupID = Tetra2Group->value(tetraGroup[0]);
		int mp_id1 = Tetra2MP->value(tetraGroup[0]);
		Region Tetra;
		std::vector<Region> adj_tetras;
		std::vector<int> visited_groups;
		bool already_visited;
		std::vector<TCellID> adj_mp;
		visited_groups.push_back(groupID);
		for (int i = 0; i < tetraGroup.size(); i++)
		{
			n_id = tetraGroup[i];
			Tetra = m_mesh->get<Region>(n_id);
			for (auto f : Tetra.get<Face>())
			{
				f.get<Region>(adj_tetras);
				for (auto tetra : adj_tetras)
				{
					int groupID2 = Tetra2Group->value(tetra.id());
					already_visited = false;
					for (auto groupID3:visited_groups)
					{
						if (groupID2 == groupID3)
							already_visited = true;
					}
					if (!already_visited)
					{
						if (groupID2 < groupID)
						{
							int mp_id2 = Tetra2MP->value(tetra.id());
							Edge me = m_ma->newMedEdge(mp_id1, mp_id2);
							visited_groups.push_back(groupID2);
						}
					}
				}
			}
		}

	}
	// Create the medial faces
	std::cout<<"> Creating the medial faces"<<std::endl;
	for (auto n_id:m_mesh->edges())
	{
		if (isAnInternalEdge(n_id)) {
			Edge e = m_mesh->get<Edge>(n_id);
			std::vector<Region> adj_tetras = e.get<Region>();
			// Counting the number of adjacent groups
			int NbrGroups = 0;
			std::vector<int> adj_groups;
			int groupID;
			bool already_seen;
			for (auto tetra:adj_tetras)
			{
				bool already_seen = false;
				groupID = Tetra2Group->value(tetra.id());
				for (auto groupID2:adj_groups)
				{
					if (groupID2 == groupID)
						already_seen = true;
				}
				if (!already_seen)
				{
					NbrGroups+=1;
					adj_groups.push_back(groupID);
				}
			}
			if (NbrGroups >= 3)
			{
				// Sorting the tetras
				std::vector<TCellID> adj_tetras_sorted(adj_tetras.size());
				std::vector<int> tetras_to_test(adj_tetras.size());
				for (int i = 1; i < adj_tetras.size(); i++) {
					tetras_to_test[i] = adj_tetras[i].id();
				}
				tetras_to_test[0] = -1;
				adj_tetras_sorted[0] = adj_tetras[0].id();
				for (int i = 0; i < adj_tetras.size() - 1; i++) {
					for (int j = 1; j < adj_tetras.size(); j++) {
						if (tetras_to_test[j] >= 0) {
							if (shareFace(adj_tetras_sorted[i], tetras_to_test[j])) {
								adj_tetras_sorted[i + 1] = tetras_to_test[j];
								tetras_to_test[j] = -1;
								break;
							}
						}
					}
				}
				// Tetra sorted,filtering the tetras
				std::vector<TCellID> adj_tetras_filtered;
				std::vector<int> groups;
				for (auto id:adj_tetras_sorted)
				{
					groupID = Tetra2Group->value(id);
					already_seen = false;
					for (auto group:groups)
					{
						if (groupID == group)
							already_seen = true;
					}
					if (!already_seen)
					{
						adj_tetras_filtered.push_back(id);
						groups.push_back(groupID);
					}

				}
				// Tetras filtered, creating the medial faces
				std::vector<TCellID> mp(adj_tetras_filtered.size());
				for (int i = 0; i < adj_tetras_filtered.size(); i++) {
					mp[i] = Tetra2MP->value(adj_tetras_filtered[i]);
				}
				Face mf = m_ma->newMedFace(mp);
			}

		}
	}

	// Set medial edges type
	m_ma->setMedialEdgeType();
	// Set medial points type
	m_ma->setMedialPointType();
	// Set medial faces type
	m_ma->setMedialFaceType();
	// Set medial surface IDs on faces
	m_ma->setMedialSurfaceIdOnFaces();
	// Set medial surface type on faces
	m_ma->setMedialSurfaceTypeOnFaces();

	return MedialAxis3DBuilder::SUCCESS;
}
/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
