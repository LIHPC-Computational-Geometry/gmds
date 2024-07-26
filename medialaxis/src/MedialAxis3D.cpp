//
// Created by chenyt on 25/04/24.
//
/*----------------------------------------------------------------------------*/
#include "gmds/medialaxis/MedialAxis3D.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/medialaxis/MedialAxisMath.h"
#include <stack>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace medialaxis {
/*----------------------------------------------------------------------------*/
MedialAxis3D::MedialAxis3D(){
	m_mesh_representation = new Mesh(MeshModel(DIM3 | F | E | N | //R |
	                                            F2E | E2F | //R2F | F2R | R2E | E2R | R2N | N2R |
	                                           F2N | E2N | N2E | N2F));
	// Medial points type: tangency of their associated triangle (1: end point, 2: regular point, 3: intersection point)
	m_mesh_representation->newVariable<int,GMDS_NODE>("medial_point_type");
	// Medial edges type: (1: end edge, 2: regular edge, 3 or more: intersection edge)
	m_mesh_representation->newVariable<int,GMDS_EDGE>("medial_edge_type");
	// Medial faces type: (1: end face, 2: regular face, 3 or more: intersection face)
	m_mesh_representation->newVariable<int,GMDS_FACE>("medial_face_type");
	// Medial surface IDs on faces
	m_mesh_representation->newVariable<int,GMDS_FACE>("medial_surface_id");
	// Medial surface IDs on points
	m_mesh_representation->newVariable<int,GMDS_NODE>("medial_surface_id_on_points");
	// Medial surface type on faces (0 if interior, 1 if adjacent to the boundary)
	m_mesh_representation->newVariable<int,GMDS_FACE>("medial_surface_type");
	// Medial surface size on faces (number of faces that compose its surface)
	m_mesh_representation->newVariable<int,GMDS_FACE>("medial_surface_size");
	// Medial radius
	m_mesh_representation->newVariable<double,GMDS_NODE>("medial_radius");
	// Mark with 1 the faces belonging to surfaces corresponding to details
	m_mesh_representation->newVariable<int,GMDS_FACE>("face_corresponding_to_details");
	// Mark with 1 the points belonging to surfaces corresponding to details
	m_mesh_representation->newVariable<int,GMDS_NODE>("point_corresponding_to_details");
}
MedialAxis3D::~MedialAxis3D()
{
	if(m_mesh_representation!= nullptr)
		delete m_mesh_representation;
}
/*----------------------------------------------------------------------------*/
Node MedialAxis3D::newMedPoint(const math::Point &APnt)
{
	Node n = m_mesh_representation->newNode(APnt);
	//std::vector<Node> nvF = {n, n, n,};
	//m_mesh_representation->newFace(nvF);
	//m_mesh_representation->newTet(n, n, n, n);
	return n;
}

/*----------------------------------------------------------------------------*/
Edge MedialAxis3D::newMedEdge(const TCellID &AN1, const TCellID &AN2)
{
	Edge e = m_mesh_representation->newEdge(AN1, AN2);
	return e;
}

/*----------------------------------------------------------------------------*/
Face MedialAxis3D::newMedFace(const std::vector<TCellID> &ANodes)
{
	Face f = m_mesh_representation->newFace(ANodes);
	return f;
}

/*----------------------------------------------------------------------------*/
Node MedialAxis3D::getMedPoint(const gmds::TCellID &AN)
{
	Node n = m_mesh_representation->get<Node>(AN);
	return n;
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::deleteMedPoint(const gmds::TCellID &AN)
{
	Node medPoint = m_mesh_representation->get<Node>(AN);
	for (auto face:medPoint.get<Face>())
	{
		m_mesh_representation->deleteFace(face);
	}
	for (auto edge:medPoint.get<Edge>())
	{
		m_mesh_representation->deleteEdge(edge);
	}
	m_mesh_representation->deleteNode(medPoint);
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::deleteMedEdge(const gmds::TCellID &AE)
{
	m_mesh_representation->deleteEdge(AE);
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::deleteMedFace(const gmds::TCellID &AF)
{
	m_mesh_representation->deleteFace(AF);
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::averagePoint(const gmds::TCellID &AN, const std::vector<TCellID> ANodes)
{
	Node n = m_mesh_representation->get<Node>(AN);
	if (ANodes.size()>0)
	{
		double x = 0.;
		double y = 0.;
		double z = 0.;
		for (auto nodeID:ANodes)
		{
			Node node = m_mesh_representation->get<Node>(nodeID);
			x = x + node.X();
			y = y + node.Y();
			z = z + node.Z();
		}
		n.X() = x/double(ANodes.size());
		n.Y() = y/double(ANodes.size());
		n.Z() = z/double(ANodes.size());
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::updateConnectivity()
{
	MeshDoctor doc(m_mesh_representation);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::setMedialPointType()
{
	// Requires medial edges type set
	std::cout<<"> Setting medial point type"<<std::endl;
	auto medEdgeType = m_mesh_representation->getVariable<int,GMDS_EDGE>("medial_edge_type");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	for (auto n_id:m_mesh_representation->nodes())
	{
		Node n = m_mesh_representation->get<Node>(n_id);
		std::vector<Edge> adj_edges = n.get<Edge>();
		medPointType->set(n_id,2);
		for (auto e:adj_edges)
		{
			if (medEdgeType->value(e.id())==1)
			{
				medPointType->set(n_id,1);
				break;
			}
			if (medEdgeType->value(e.id()) > medPointType->value(n_id))
			{
				medPointType->set(n_id,medEdgeType->value(e.id()));
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::setMedialEdgeType()
{
//	std::cout<<"> Setting E2F connectivity"<<std::endl;
//	MeshDoctor doc(m_mesh_representation);
//	doc.buildEdgesAndX2E();
//	doc.updateUpwardConnectivity();
	std::cout<<"> Setting medial edge type"<<std::endl;
	//auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	auto medEdgeType = m_mesh_representation->getVariable<int,GMDS_EDGE>("medial_edge_type");
	for (auto n_id:m_mesh_representation->edges())
	{
		Edge e = m_mesh_representation->get<Edge>(n_id);
		medEdgeType->set(e.id(),e.get<Face>().size());
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::setMedialFaceType()
{
	// Requires medial edges type set
	std::cout<<"> Setting medial face type"<<std::endl;
	auto medFaceType = m_mesh_representation->getVariable<int,GMDS_FACE>("medial_face_type");
	auto medEdgeType = m_mesh_representation->getVariable<int,GMDS_EDGE>("medial_edge_type");
	for (auto n_id:m_mesh_representation->faces())
	{
		Face f = m_mesh_representation->get<Face>(n_id);
		int faceType = 2;
		int edgeType;
		for (auto e:f.get<Edge>())
		{
			edgeType = medEdgeType->value(e.id());
			if (edgeType == 1)
			{
				faceType = 1;
				break;
			}
			if (edgeType > faceType)
				faceType = edgeType;
		}
		medFaceType->set(f.id(),faceType);
	}
}

/*----------------------------------------------------------------------------*/
int MedialAxis3D::setMedialSurfaceIdOnFaces()
{
	// Requires medial edges type set
	std::cout<<"> Setting medial surface IDs"<<std::endl;
	auto medSurfId = m_mesh_representation->getVariable<int,GMDS_FACE>("medial_surface_id");
	auto medEdgeType = m_mesh_representation->getVariable<int,GMDS_EDGE>("medial_edge_type");
	std::vector<TCellID> alreadyVisited(m_mesh_representation->getNbFaces());
	for (int i = 0; i<m_mesh_representation->getNbFaces(); i++)
		alreadyVisited[i] = 0;
	int Id = 0;
	for (auto f_id:m_mesh_representation->faces())
	{
		if (alreadyVisited[f_id] == 0)
		{
			std::stack<TCellID> surface;
			surface.push(f_id);
			alreadyVisited[f_id] = 1;
			while (!surface.empty())
			{
				TCellID i = surface.top();
				surface.pop();
				medSurfId->set(i,Id);
				Face f = m_mesh_representation->get<Face>(i);
				for (auto e:f.get<Edge>())
				{
					if (medEdgeType->value(e.id()) == 2)
					{
						for (auto face:e.get<Face>())
						{
							if (alreadyVisited[face.id()] == 0)
							{
								surface.push(face.id());
								alreadyVisited[face.id()] = 1;
							}
						}
					}
				}
			}
			Id += 1;
		}
	}
	std::cout<<"NB medial surfaces: "<<Id<<std::endl;
	return Id;
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::setMedialSurfaceIdOnPoints()
{
	// Requires medial surface IDs on faces
	auto medSurfIdOnFaces = m_mesh_representation->getVariable<int,GMDS_FACE>("medial_surface_id");
	auto medSurfIdOnPoints = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_surface_id_on_points");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	for (auto n_id:m_mesh_representation->nodes())
	{
		Node n = m_mesh_representation->get<Node>(n_id);
		std::vector<Face> adj_faces = n.get<Face>();
		if (medPointType->value(n.id()) == 1 || medPointType->value(n.id()) == 2)
			medSurfIdOnPoints->set(n.id(),medSurfIdOnFaces->value(adj_faces[0].id()));
		else
			medSurfIdOnPoints->set(n.id(),-1);
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::identifyDetailsSurfaces(const double &ATol)
{
	// Requires medial surfaces IDs, medial points type, medial radii
	auto surfIdOnPoints = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_surface_id_on_points");
	auto surfIdOnFaces = m_mesh_representation->getVariable<int,GMDS_FACE>("medial_surface_id");
	auto medPointType = m_mesh_representation->getVariable<int,GMDS_NODE>("medial_point_type");
	auto medRadius = m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius");
	auto faceCorrespondingToDetails = m_mesh_representation->getVariable<int,GMDS_FACE>("face_corresponding_to_details");
	auto pointCorrespondingToDetails = m_mesh_representation->getVariable<int,GMDS_NODE>("point_corresponding_to_details");
	// Number of surfaces
	int NbSurfaces = setMedialSurfaceIdOnFaces();
	// Vector stocking the medial points of type >= 3 surrounding the surface for each surface
	std::vector<std::vector<Node>> intersectionPoints(NbSurfaces);
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (medPointType->value(n_id) >= 3)
		{
			Node n = m_mesh_representation->get<Node>(n_id);
			std::vector<Face> adj_faces = n.get<Face>();
			std::vector<int> medialSurfaces;
			for (auto face:adj_faces)
			{
				int surfID = surfIdOnFaces->value(face.id());
				bool alreadySeen = false;
				for (auto ID:medialSurfaces)
				{
					if (ID == surfID)
						alreadySeen = true;
				}
				if (!alreadySeen)
					medialSurfaces.push_back(surfID);
			}
			for (auto ID:medialSurfaces)
				intersectionPoints[ID].push_back(n);
		}
	}
	// Vector stocking the max ranges for each surface
	std::vector<double> maxRanges(NbSurfaces);
	// Initializing max ranges to 0
	for (int i = 0; i<NbSurfaces; i++)
		maxRanges[i] = 0.;
	// Computing the max ranges
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (medPointType->value(n_id) == 2 || medPointType->value(n_id) == 1)
		{
			Node n = m_mesh_representation->get<Node>(n_id);
			int surfID = surfIdOnPoints->value(n_id);
			std::vector<Node> intPoints = intersectionPoints[surfID];
			if (intPoints.size() > 0)
			{
				double range = maxRange(intPoints[0].point(),medRadius->value(intPoints[0].id()),n.point(),medRadius->value(n.id()));
				for (auto intPoint:intPoints)
				{
					double range1 = maxRange(intPoint.point(),medRadius->value(intPoint.id()),n.point(),medRadius->value(n.id()));
					if (range1 < range)
						range = range1;
				}
				if (range > maxRanges[surfID])
					maxRanges[surfID] = range;
			}
			else
				// This means that this surface doesn't intersect any other surface, so it is not considered as corresponding to details
				maxRanges[surfID] = ATol + 1.;
		}
	}
	// Marking the faces whose surface corresponds to details
	for (auto f_id:m_mesh_representation->faces())
	{
		if (maxRanges[surfIdOnFaces->value(f_id)] < ATol)
			faceCorrespondingToDetails->set(f_id,1);
		else
			faceCorrespondingToDetails->set(f_id,0);
	}
	// Marking the points whose surface corresponds to details
	for (auto n_id:m_mesh_representation->nodes())
	{
		Node n = m_mesh_representation->get<Node>(n_id);
		if (medPointType->value(n.id()) != 1 && medPointType->value(n.id()) != 2)
			pointCorrespondingToDetails->set(n.id(),0);
		else
			pointCorrespondingToDetails->set(n.id(),faceCorrespondingToDetails->value(n.get<Face>()[0].id()));
	}
}

/*----------------------------------------------------------------------------*/
int MedialAxis3D::correspondsToDetail(const gmds::TCellID &APointID)
{
	auto pointCorrespondingToDetails = m_mesh_representation->getVariable<int,GMDS_NODE>("point_corresponding_to_details");
	return (pointCorrespondingToDetails->value(APointID));
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::setMedialRadius(const TCellID &APointID, double AValue)
{
	auto medRadius = m_mesh_representation->getVariable<double,GMDS_NODE>("medial_radius");
	Node n = m_mesh_representation->get<Node>(APointID);
	medRadius->set(n.id(),AValue);
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::setMedialSurfaceTypeOnFaces()
{
	// Requires medial face type and medial surfaces IDs set
	auto medFaceType = m_mesh_representation->getVariable<int,GMDS_FACE>("medial_face_type");
	auto medSurfId = m_mesh_representation->getVariable<int,GMDS_FACE>("medial_surface_id");
	auto medSurfType = m_mesh_representation->getVariable<int,GMDS_FACE>("medial_surface_type");
	for (auto f_id:m_mesh_representation->faces())
	{
		medSurfType->set(f_id,0);
	}
	for (auto f_id:m_mesh_representation->faces())
	{
		if (medFaceType->value(f_id) == 1)
		{
			for (auto f2_id:m_mesh_representation->faces())
			{
				if (medSurfId->value(f_id) == medSurfId->value(f2_id))
					medSurfType->set(f2_id,1);
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::setMedialSurfaceSizeOnFaces()
{
	// Requires medial surface ID
	auto medSurfId = m_mesh_representation->getVariable<int,GMDS_FACE>("medial_surface_id");
	auto medSurfSize = m_mesh_representation->getVariable<int,GMDS_FACE>("medial_surface_size");
	int NbSurfaces = setMedialSurfaceIdOnFaces();
	for (int Id = 0; Id<NbSurfaces; Id++)
	{
		int size = 0;
		for (auto face_id:m_mesh_representation->faces())
		{
			if (medSurfId->value(face_id) == Id)
				size += 1;
		}
		for (auto face_id:m_mesh_representation->faces())
		{
			if (medSurfId->value(face_id) == Id)
				medSurfSize->set(face_id,size);
		}
	}
}

/*----------------------------------------------------------------------------*/
void MedialAxis3D::flattenUnwantedFaces()
{
	int NbFacesFlattened = 0;
	for (auto face_id:m_mesh_representation->faces())
	{
		Face face = m_mesh_representation->get<Face>(face_id);
		std::vector<Node> adj_nodes = face.get<Node>();
		std::vector<double> criterium;
		for (int i = 0; i<adj_nodes.size(); i++)
		{
			double c = 0.;
			double max_dist = 0.;
			for (auto n:adj_nodes)
			{
				double dist = n.point().distance(adj_nodes[i].point());
				c += dist;
				if (dist > max_dist)
					max_dist = dist;
			}
			c -= max_dist;
			c = fabs(c);
			criterium.push_back(c);
		}
		// test
		if (face_id == 12994)
		{
			std::cout<<"Criterium : ";
			for (auto c:criterium)
				std::cout<<c<<" ";
			std::cout<<std::endl;
			std::cout<<adj_nodes.size()<<std::endl;
			std::cout<<"Points: ";
			for (auto n:adj_nodes)
				std::cout<<n.id()<<" ";
			std::cout<<std::endl;
			std::cout<<adj_nodes.size()<<std::endl;
			std::cout<<"Distances : "<<m_mesh_representation->get<Node>(73877).point().distance(m_mesh_representation->get<Node>(72652).point())<<std::endl;
		}
		// Find the two nodes with the highest values of c
		int max1 = 0;
		int max2 = 0;
		double max_c = 0;
		for (int i = 0; i<adj_nodes.size(); i++)
		{
			if (criterium[i] > max_c)
			{
				max1 = i;
				max_c = criterium[i];
			}
		}
		max_c = 0;
		for (int i = 0; i<adj_nodes.size(); i++)
		{
			if (criterium[i] > max_c && i != max1)
			{
				max2 = i;
				max_c = criterium[i];
			}
		}
		// test
		if (face_id == 12994)
		{
			std::cout<<"Max: "<< criterium[max1] << " "<<criterium[max2]<<std::endl;
		}
		if (criterium[max1] > 0 && criterium[max2]/criterium[max1] < 1./10.)
		{
			TCellID badPoint = adj_nodes[max1].id();
			std::vector<TCellID> goodPoints;
			for (auto n:adj_nodes)
			{
				if (n.id() != badPoint)
					goodPoints.push_back(n.id());
			}
			averagePoint(badPoint,goodPoints);
			NbFacesFlattened += 1;
		}
	}
	std::cout<<"NB faces flattened : "<<NbFacesFlattened<<std::endl;
}


/*----------------------------------------------------------------------------*/
void MedialAxis3D::write(std::string AFileName)
{
	IGMeshIOService ioService(m_mesh_representation);
	VTKWriter vtkWriter(&ioService);
	//vtkWriter.setCellOptions(N| E| F| R);
	//vtkWriter.setDataOptions(N| E| F| R);
	std::cout<<"> Writing the medial axis"<<std::endl;
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
	std::cout<<"NB medial points : "<<m_mesh_representation->getNbNodes()<<std::endl;
	std::cout<<"NB medial edges : "<<m_mesh_representation->getNbEdges()<<std::endl;
	std::cout<<"NB medial faces : "<<m_mesh_representation->getNbFaces()<<std::endl;
}

/*----------------------------------------------------------------------------*/
}  // end namespace medialaxis
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
      /*----------------------------------------------------------------------------*/