//
// Created by rochec on 13/03/2023.
//

/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/claire/MeshAlignment_2D.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

MeshAlignment_2D::MeshAlignment_2D(Mesh *AMeshTri, Variable<math::Vector3d>* A_VectorField, Mesh *AMeshQuad) :
  m_meshTri(AMeshTri),
  m_fl(m_meshTri),
  m_meshQuad(AMeshQuad),
  m_VectorField(A_VectorField)
{
	m_var_deviation = m_meshQuad->newVariable<double, GMDS_NODE>("MeshAlignment");
}
/*------------------------------------------------------------------------*/
MeshAlignment_2D::~MeshAlignment_2D()
{
	m_meshQuad->deleteVariable(GMDS_NODE, "MeshAlignment");
}
/*------------------------------------------------------------------------*/
MeshAlignment_2D::STATUS
MeshAlignment_2D::execute()
{

	// Compute the flow deviation
	for (auto const &n_id:m_meshQuad->nodes())
	{
		Node n = m_meshQuad->get<Node>(n_id);
		std::vector<Face> n_faces_in_quad = n.get<Face>() ;

		if (n_faces_in_quad.size() > 2) {

			gmds::Cell::Data data = m_fl.find(n.point());

			Node n_closer_tri = m_meshTri->get<Node>(data.id);
			std::vector<Face> faces = n_closer_tri.get<Face>();

			bool triangle_found(false);
			TCellID f_id_in_Tri_mesh(NullID);
			for (auto const &f_tri : faces) {
				std::vector<Node> f_nodes = f_tri.get<Node>();
				triangle_found = math::Utils::isInTriangle(f_nodes[0].point(), f_nodes[1].point(), f_nodes[2].point(), n.point());
				if (triangle_found) {
					f_id_in_Tri_mesh = f_tri.id();
				}
			}

			if (f_id_in_Tri_mesh != NullID) {
				math::Vector3d v_ref = ComputeInterpolatedVector(f_id_in_Tri_mesh, n.point());
				m_var_deviation->set(n_id, maxAngleDeviationAtNode(n_id, v_ref) );
			}

		}

	}

	return MeshAlignment_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/
Variable<double>*
MeshAlignment_2D::getVariableDeviation()
{
	return m_var_deviation;
}
/*------------------------------------------------------------------------*/
math::Vector3d
MeshAlignment_2D::ComputeInterpolatedVector(TCellID tri_id, const math::Point& p)
{
	math::Vector3d v_ref;
	Face f = m_meshTri->get<Face>(tri_id);
	std::vector<Node> f_nodes = f.get<Node>();
	v_ref.setX(math::Utils::linearInterpolation2D3Pt(f_nodes[0].point(), f_nodes[1].point(), f_nodes[2].point(), p,
	                                                 m_VectorField->value(f_nodes[0].id()).X(), m_VectorField->value(f_nodes[1].id()).X(),
	                                                 m_VectorField->value(f_nodes[2].id()).X()));
	v_ref.setY(math::Utils::linearInterpolation2D3Pt(f_nodes[0].point(), f_nodes[1].point(), f_nodes[2].point(), p,
	                                                 m_VectorField->value(f_nodes[0].id()).Y(), m_VectorField->value(f_nodes[1].id()).Y(),
	                                                 m_VectorField->value(f_nodes[2].id()).Y()));

	v_ref.setZ(0.0);

	v_ref.normalize();

	return v_ref;
}
/*------------------------------------------------------------------------*/
double
MeshAlignment_2D::minAngleBtwVectorandCross(const math::Vector3d& v, const math::Vector3d& v_cross)
{
	double min_angle(90);

	math::Vector3d v_cross_ortho({-v_cross.Y(), v_cross.X(), 0.0});
	v_cross_ortho.normalize();

	if (acos(v.dot(v_cross)) * 180 / M_PI < min_angle) {
		min_angle = acos(v.dot(v_cross)) * 180 / M_PI;
	}
	if (acos(v.dot(-v_cross)) * 180 / M_PI < min_angle) {
		min_angle = acos(v.dot(-v_cross)) * 180 / M_PI;
	}
	if (acos(v.dot(v_cross_ortho)) * 180 / M_PI < min_angle) {
		min_angle = acos(v.dot(v_cross_ortho)) * 180 / M_PI;
	}
	if (acos(v.dot(-v_cross_ortho)) * 180 / M_PI < min_angle) {
		min_angle = acos(v.dot(-v_cross_ortho)) * 180 / M_PI;
	}

	return min_angle;
}
/*------------------------------------------------------------------------*/
double
MeshAlignment_2D::maxAngleDeviationAtNode(TCellID n_id, const math::Vector3d& v_ref)
{
	double max_deviation(0);

	Node n = m_meshQuad->get<Node>(n_id);
	std::vector<Edge> edges = n.get<Edge>();
	std::vector<double> min_angles;

	for (auto const &e : edges) {
		std::vector<Node> e_nodes = e.get<Node>();
		math::Vector3d vec_edge = (e_nodes[0].point() - e_nodes[1].point()).normalize();
		double min_angle = minAngleBtwVectorandCross(vec_edge, v_ref) ;
		min_angles.push_back(min_angle);
	}

	for (auto angle:min_angles)
	{
		if (angle > max_deviation)
		{
			max_deviation = angle;
		}
	}

	return max_deviation;
}
/*------------------------------------------------------------------------*/