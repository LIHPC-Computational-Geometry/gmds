/*----------------------------------------------------------------------------*/
#include <gmds/smoothy/EllipticSmoother2D.h>
#include <gmds/smoothy/EllipticSmoothing.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::smoothy;
/*----------------------------------------------------------------------------*/
EllipticSmoother2D::EllipticSmoother2D(Mesh *AMesh)
: m_mesh(AMesh), m_max_iterations(1000), m_theta(1e-3)
{
	m_lock= m_mesh->newMark<Node>();
}
/*----------------------------------------------------------------------------*/
EllipticSmoother2D::~EllipticSmoother2D(){
	m_mesh->unmarkAll<Node>(m_lock);
	m_mesh->freeMark<Node>(m_lock);
}
/*----------------------------------------------------------------------------*/
void EllipticSmoother2D::initLock()
{
	m_mesh->unmarkAll<Node>(m_lock);
}
/*----------------------------------------------------------------------------*/
void
EllipticSmoother2D::lock(const TInt AMark)
{
	initLock();
	for(auto i:m_mesh->nodes())
		if(m_mesh->isMarked<Node>(i,AMark))
			m_mesh->mark<Node>(i,m_lock);
}
/*----------------------------------------------------------------------------*/
void
EllipticSmoother2D::lock(const std::vector<TCellID> &AV)
{
	initLock();
	for(auto i:AV)
		m_mesh->mark<Node>(i,m_lock);
}
/*----------------------------------------------------------------------------*/
void
EllipticSmoother2D::lockBoundary()
{
	initLock();
	BoundaryOperator2D op(m_mesh);
	auto mark_node_NAN = m_mesh->newMark<Node>();
	auto mark_node_on_pnt =  m_mesh->newMark<Node>();
	auto mark_node_on_crv =  m_mesh->newMark<Node>();
	auto mark_edge_on_crv =  m_mesh->newMark<Edge>();

	op.markCellOnGeometry(mark_edge_on_crv, mark_node_on_crv, mark_node_on_pnt, mark_node_NAN);
	for (auto n_id :  m_mesh->nodes()) {
		if ( m_mesh->isMarked<Node>(n_id, mark_node_on_crv) ||  m_mesh->isMarked<Node>(n_id, mark_node_on_pnt)) {
			m_mesh->mark<Node>(n_id, m_lock);
		}
	}
}
/*----------------------------------------------------------------------------*/
bool EllipticSmoother2D::isValid() const
{
	MeshModel model = m_mesh->getModel();
	if(model.has(R)) {
		return false;
	}
	if (!model.has(F) ||
	    !model.has(N) ||
	    !model.has(F2N) ||
	    !model.has(N2F)  ){
		return false;
	}
	return true;
}
/*----------------------------------------------------------------------------*/
void EllipticSmoother2D::execute()
{
	//------------------------------------------------------------------------
	//	Step 1: prepare data for the EllipticSmoothing requirement
	//------------------------------------------------------------------------
	std::vector<double> node_coords(2*m_mesh->getNbNodes());
	std::vector<bool> lock_nodes(2*m_mesh->getNbNodes());
	for(auto n_id:m_mesh->nodes()){
		math::Point pi = m_mesh->get<Node>(n_id).point();
		node_coords[2*n_id]=pi.X();
		node_coords[2*n_id+1]=pi.Y();
		lock_nodes[2*n_id]= (m_mesh->isMarked<Node>(n_id, m_lock));
		lock_nodes[2*n_id+1]= (m_mesh->isMarked<Node>(n_id, m_lock));
	}

	std::vector<std::array<int, 3>> triangles(4*m_mesh->getNbFaces());
	constexpr int quadrature[4][3]= {{0,1,3},{1,2,0},{2,3,1},{3,0,2}};
	for(auto f_id:m_mesh->faces()){
		std::vector<TCellID> nids = m_mesh->get<Face>(f_id).getIDs<Node>();
		for(auto i=0;i<4;i++) {
			triangles[4 * f_id + i] = {static_cast<int>(nids[quadrature[i][0]]),
			                           static_cast<int>(nids[quadrature[i][1]]),
			                           static_cast<int>(nids[quadrature[i][2]])};
		}
	}
	EllipticSmoothingMeshGlue2D tri(node_coords,triangles,lock_nodes);

	EllipticSmoothingOptions options = ESO_default2D;
	options.maxiter=m_max_iterations;
	options.eps_from_theorem=true;
	options.theta = m_theta;
	options.bfgs_threshold = 1e-9;
	options.debug=0;
	EllipticSmoothing2D smooth(tri, triangles.size(),options);
	//	smoother2D.start_eps = 1E-4;
	smooth.execute();
	tri.get_verts(node_coords);
	for(auto n_id:m_mesh->nodes()){
		m_mesh->get<Node>(n_id).setXYZ(node_coords[2*n_id],node_coords[2*n_id+1],.0);
	}
}
/*----------------------------------------------------------------------------*/
