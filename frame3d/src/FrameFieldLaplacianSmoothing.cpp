/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// STL Header files
#include <iostream>
#include <map>
#include <set>
/*---------------------------------------------------------------------------*/
#include <gmds/frame3d/FrameFieldLaplacianSmoothing.h>
/*---------------------------------------------------------------------------*/
// GMDS Header files
#include "gmds/math/Matrix.h"
#include "gmds/math/Cross.h"
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
FrameFieldLaplacianSmoothing::FrameFieldLaplacianSmoothing(
	gmds::Mesh* AMesh,
	gmds::Variable<gmds::math::Quaternion>*& AField,
	gmds::Variable<gmds::math::Vector3d>*& ANormal)
	: m_mesh(AMesh), m_frame_field(AField), m_surf_normal(ANormal),
	m_mark_candidates(0), m_candidate_mark_initialized(false),
	m_boundary_marks_initialized(false)
{
	m_mesh = AMesh;
}
/*---------------------------------------------------------------------------*/
void FrameFieldLaplacianSmoothing::selectNodes(const int AMark)
{
	m_mark_candidates = AMark;
	m_candidate_mark_initialized = true;
}
/*---------------------------------------------------------------------------*/
void FrameFieldLaplacianSmoothing::initBoundaryMarks(
	const int AMarkPnt, const int AMarkCurve,
	const int AMarkSurf)
{
	m_mark_point = AMarkPnt;
	m_mark_curve = AMarkCurve;
	m_mark_surface = AMarkSurf;
	m_boundary_marks_initialized = true;
}
/*---------------------------------------------------------------------------*/
void FrameFieldLaplacianSmoothing::execute()
{
	if (!m_candidate_mark_initialized)
		throw GMDSException("Candidate mark is not initialized");
	if (!m_boundary_marks_initialized)
		throw GMDSException("Boundary marks are not initialized");

	smooth();
}

/*---------------------------------------------------------------------------*/
void FrameFieldLaplacianSmoothing::initCandidates()
{
	int index = 0;
	m_candidates.clear();
	for (auto n_id: m_mesh->nodes())
	{
		if (m_mesh->isMarked<Node>(n_id, m_mark_candidates))
		{
			m_candidates.push_back(m_mesh->get<Node>(n_id));
			m_candidates_index[n_id] = index++;
		}
	}
	std::cout << "Nb candidates (" << m_candidates.size()
		<< " / " << m_mesh->getNbNodes() << ")" << std::endl;

	//Definition of the neighborhood
	m_candidates_neighborhood.clear();
	m_candidates_neighborhood.resize(m_candidates.size());

	for (unsigned int i = 0; i < m_candidates.size(); i++)
	{
		Node n = m_candidates[i];
		std::vector<Edge> adj_edges = n.get<Edge>();

		for (unsigned int j = 0; j < adj_edges.size(); j++)
		{
			Edge e_j = adj_edges[j];

			std::vector<Node> e_j_nodes = e_j.get<Node>();
			Node other_node =
				(e_j_nodes[0].id() == n.id()) ? e_j_nodes[1] : e_j_nodes[0];

			m_candidates_neighborhood[i].push_back(other_node);
		}
	}

}
/*---------------------------------------------------------------------------*/
void FrameFieldLaplacianSmoothing::smooth(gmds::Node& ANode)
{
	int candidate_index = m_candidates_index[ANode.id()];

	std::vector<Node> adj_nodes = m_candidates_neighborhood[candidate_index];

	int nb_quaternions = adj_nodes.size();

	bool onSurface = false;
	math::Vector3d normal;
	math::Quaternion qref = (*m_frame_field)[ANode.id()];
	if (m_mesh->isMarked(ANode, m_mark_surface))
	{
		onSurface = true;
		normal = (*m_surf_normal)[ANode.id()];
	}
	if (onSurface)
	{
		//current cross
		math::Cross cref(qref, normal);

		//we get the adjacent quaternions associated to nodes that are not
		//in the volume
		std::vector<math::Quaternion> quats;
		std::vector<TCoord> coefs;
		for (int i = 0; i < adj_nodes.size(); i++)
		{
			Node n = adj_nodes[i];
			if ((m_mesh->isMarked(n, m_mark_surface) ||
				m_mesh->isMarked(n, m_mark_curve) ) &&
				!m_mesh->isMarked(n, m_mark_point) )
			{
				math::Quaternion q = (*m_frame_field)[n.id()];

				q = q.alignWith(normal);
				quats.push_back(q);
				coefs.push_back(1.0);
			}
		}

		math::Quaternion res = math::Quaternion::mean(quats, coefs);

		
			res = res.alignWith(normal);

		(*m_frame_field)[ANode.id()] = res;

	}
}
/*---------------------------------------------------------------------------*/
void FrameFieldLaplacianSmoothing::smooth()
{
	//We collect the nodes we want to smooth the quaternion on
	initCandidates();

	int nb_steps = 3;
	for (unsigned int i = 0; i < m_candidates.size(); i++){
		smooth(m_candidates[i]);
	}
}

