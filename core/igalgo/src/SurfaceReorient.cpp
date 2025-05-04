/*----------------------------------------------------------------------------*/
#include <gmds/igalgo/SurfaceReorient.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
SurfaceReorient::SurfaceReorient(Mesh* AMesh, const TInt ADim)
        :m_mesh(AMesh), m_dim(ADim)
{}
/*----------------------------------------------------------------------------*/
SurfaceReorient::~SurfaceReorient()
{}
/*----------------------------------------------------------------------------*/
bool SurfaceReorient::isValid() const
{
    if(m_dim==3)
        return (m_mesh->getModel()==(DIM3|F|N|F2N|N2F));
    else if(m_dim==2)
        return (m_mesh->getModel()==(DIM3|F|N|F2N) ||
                m_mesh->getModel()==(DIM2|F|N|F2N));

    //dimension error
    return false;
}
/*----------------------------------------------------------------------------*/
int SurfaceReorient::execute() {
	 return (m_dim==2)?orient2d():orient3d();
}
/*----------------------------------------------------------------------------*/
int SurfaceReorient::orient2d()
{
	 auto nb_reorientation =0;
	 for(auto f_id:m_mesh->faces()){
		  Face f = m_mesh->get<Face>(f_id);
		  if (orient2d(f))
			  nb_reorientation++;
	 }
	 return nb_reorientation;
}
/*------------------------------------------------------------------------*/
bool SurfaceReorient::orient2d(Face& AF)
{
	 bool isReoriented = false;
	 std::vector<Node> nodes = AF.get<Node>();
	 TCoord orientation=0;
	 if(AF.type()==GMDS_TRIANGLE) {
		  orientation = isLeft(nodes[0],nodes[1],nodes[2]);
	 }
	 else {
		  //find the rightmost lowest vertex of the polygon
		  unsigned int index_min=0;
		  TCoord x_min = nodes[0].X();
		  TCoord y_min = nodes[0].Y();
		  for(unsigned int
		          i=0;i<nodes.size();i++)
		  {
			  if(nodes[i].Y()>y_min)
				  continue;
			  if(nodes[i].Y()==y_min) {	// just as low
				  if(nodes[i].X()<x_min)  // and to left
					  continue;
			  }

			  index_min =i;
			  x_min = nodes[i].X();
			  y_min = nodes[i].Y();
		  }

		  if(index_min==0)
			  orientation = isLeft(nodes[nodes.size()-1],nodes[0],nodes[1]);
		  else if (index_min==nodes.size()-1)
			  orientation = isLeft(nodes[index_min-1],nodes[index_min],nodes[0]);
		  else
			  orientation = isLeft(nodes[index_min-1],nodes[index_min],nodes[index_min+1]);
	 }
	 if(orientation>0.0) // clockwise or degenerated (=0)
	 {
		  isReoriented= true;
		  std::vector<Node> nodes_inv;
		  nodes_inv.resize(nodes.size());
		  auto node_size = nodes.size();
		  for(unsigned int i=0;i<node_size;i++)
			  nodes_inv[i] = nodes[node_size-1-i];
		  AF.set<Node>(nodes_inv);
	 }
	 return isReoriented;
}
/*----------------------------------------------------------------------------*/
TCoord SurfaceReorient::isLeft(Node& AN1, Node& AN2, Node& AN3)
{
	 return ( (AN2.X()-AN1.X()) * (AN3.Y()-AN1.Y()) -
	         (AN3.X()-AN1.X()) * (AN2.Y()-AN1.Y()) );
}
/*----------------------------------------------------------------------------*/
int SurfaceReorient::orient3d()
{
	 auto nb_reorientation =0;
	 auto mark_done = m_mesh->newMark<Face>();

	 bool keep_working =true;
	 while (keep_working) {
		 //We look for the first unmarked face

		 auto it_face = m_mesh->faces_begin();
		 while (m_mesh->isMarked<Face>(*it_face,mark_done) &&
		    it_face!=m_mesh->faces_end()){
			 ++it_face;
		 }
		 if(it_face==m_mesh->faces_end()){
			 keep_working=false;
		 }
		 if(keep_working) {
			 // We take the first face as the orientation reference seed
			 auto seed_id = *it_face;
			 Face seed = m_mesh->get<Face>(seed_id);
			    // And now, we go through all the faces in an advancing-front manner to
			    // orient faces according to the seed
			    std::vector<TCellID> front;
			 front.push_back(seed.id());
			 m_mesh->mark(seed, mark_done);
			 while (!front.empty()) {
				 auto current_id = front.back();
				 front.pop_back();
				 auto current_face = m_mesh->get<Face>(current_id);
				 // the current face is considered as well oriented and we orient the face sharing
				 //  an edge with it in a valid manner
				 std::vector<Node> node_faces = current_face.get<Node>();
				 for (auto i = 0; i < node_faces.size(); i++) {
					 Node ni = node_faces[i];
					 Node nj = node_faces[(i + 1) % node_faces.size()];
					 auto fij_ids = m_mesh->getCommonFaces(ni, nj);
					 if (fij_ids.size() == 2) {
						 // means there is another face
						 auto other_face_id = (fij_ids[0] == current_id) ? fij_ids[1] : fij_ids[0];
						 auto other_face = m_mesh->get<Face>(other_face_id);
						 if (!m_mesh->isMarked(other_face, mark_done)) {
							 std::vector<TCellID> other_face_node_ids = other_face.getIDs<Node>();
							 // we check if the node ni and nj are traversed in the same way.
							 bool same_direction = false;
							 for (auto k = 0; k < other_face_node_ids.size(); k++) {

								 auto id_k = other_face_node_ids[k];
								 auto id_l = other_face_node_ids[(k + 1) % other_face_node_ids.size()];
								 if (id_k == ni.id() && id_l == nj.id()) same_direction = true;
							 }
							 if (same_direction) {
								 std::reverse(other_face_node_ids.begin(), other_face_node_ids.end());
								 other_face.set<Node>(other_face_node_ids);
								 nb_reorientation++;
							 }
							 m_mesh->mark(other_face, mark_done);
							 front.push_back(other_face_id);
						 }
					 }
				 }
			 }
		 }
	 }
	 m_mesh->negateMaskMark<Node>(mark_done);
	 m_mesh->freeMark<Node>(mark_done);
	 //warning, we must ensure that the mesh is reoriented even if it has several
	 //connex part
	 return nb_reorientation;
}