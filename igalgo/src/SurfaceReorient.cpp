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
        return (m_mesh->getModel()==(DIM3|R|N|R2N));
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
	 //We take the first face as the orientation reference seed
	 auto seed = m_mesh->get<Face>(*m_mesh->faces_begin());
	 // And now, we go through all the faces in an advancing-front manner to
	 // orient faces according to the seed
	 auto mark_done = m_mesh->newMark<Face>();


	 //warning, we must ensure that the mesh is reoriented even if it has several
	 //connex part
	 return nb_reorientation;
}