/*----------------------------------------------------------------------------*/
#include <gmds/igalgo/BoundaryOperator.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
BoundaryOperator::BoundaryOperator(Mesh* AMesh, double AAngle)
        :m_mesh(AMesh), m_surface_angle_dot(AAngle)
{}
/*----------------------------------------------------------------------------*/
BoundaryOperator::~BoundaryOperator()
= default;
/*----------------------------------------------------------------------------*/
void BoundaryOperator::setSurfaceAngleDot(double AD) {
    m_surface_angle_dot=AD;
}
/*----------------------------------------------------------------------------*/
double BoundaryOperator::getSurfaceAngleDot() const {
    return m_surface_angle_dot;
}
/*----------------------------------------------------------------------------*/
bool BoundaryOperator::isValid() const
{
    MeshModel model = m_mesh->getModel();
    if(model.has(R)) {
        if (!model.has(F2R))
            return false;
        if (!model.has(F2E))
            return false;
        if (!model.has(E2F))
            return false;

        return true;
    }
    else if(model.has(F)) {
        if (!model.has(F2N))
            return false;
        if (!model.has(N2F))
            return false;

        return true;
    }

    return false;
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::getBoundaryNodes(std::vector<TCellID>& ANodeIDs)
{
    std::set<TCellID> elts;
    MeshModel model = m_mesh->getModel();
    if(model.has(R)) {
        throw GMDSException("Not yet implemented in 3D");
    }
    else {
        for (auto e_id : m_mesh->edges())  {
            Edge e = m_mesh->get<Edge>(e_id);
            std::vector<TCellID> faces_e = e.getIDs<Face>();
            if(faces_e.size()==1){
                std::vector<TCellID> nodes_e = e.getIDs<Node>();
                elts.insert(nodes_e.begin(),nodes_e.end());
            }
        }

    }
    ANodeIDs.clear();
    ANodeIDs.insert(ANodeIDs.end(),elts.begin(),elts.end());
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
markCellOnGeometry(int AMarkFOnSurf,
                   int AMarkEOnSurf,
                   int AMarkNOnSurf,
                   int AMarkEOnCurve,
                   int AMarkNOnCurve,
                   int AMarkNOnPnt,
                   int AMarkAN)
{
    MeshModel model = m_mesh->getModel();
    if(model.has(R)) {
        markCellsOnSurfaces(AMarkFOnSurf, AMarkEOnSurf, AMarkNOnSurf);
        markCellsOnCurves(AMarkFOnSurf, AMarkEOnSurf, AMarkEOnCurve, AMarkNOnCurve);
        markNodesOnPoint(AMarkEOnCurve, AMarkNOnCurve,AMarkNOnPnt);
        markAloneNodes(AMarkAN);
        colorFaces(AMarkFOnSurf,AMarkEOnCurve);
    }
    else {

        //all the faces, edges and nodes are on the surface, so we invert the mesh
        //mask
        m_mesh->negateMaskMark<Node>(AMarkNOnSurf);
        m_mesh->negateMaskMark<Edge>(AMarkEOnSurf);
        m_mesh->negateMaskMark<Face>(AMarkFOnSurf);
        if(model.has(DIM2)) {
            markCellsOnCurves(AMarkEOnCurve, AMarkNOnCurve);

        }
        else{
            markCellsOnCurves(AMarkFOnSurf, AMarkEOnSurf, AMarkEOnCurve, AMarkNOnCurve);

        }
        markNodesOnPoint(AMarkEOnCurve, AMarkNOnCurve,AMarkNOnPnt);

        markAloneNodes(AMarkAN);

        colorFaces(AMarkFOnSurf,AMarkEOnCurve);
    }
}

/*----------------------------------------------------------------------------*/
void BoundaryOperator::markAloneNodes(int AMarkAlone)
{
    int cpt = 0;
    for (auto n_id:m_mesh->nodes())
    {
        Node n = m_mesh->get<Node>(n_id);
        if (n.get<Edge>().empty()) {
            m_mesh->mark(n, AMarkAlone);
            cpt++;
        }

    }
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
markCellsOnSurfaces(int AMarkBF, int AMarkBE, int AMarkBN)
{
    int cpt1 = 0, cpt2 = 0;

    for (auto f_id:m_mesh->faces())
    {
        Face f = m_mesh->get<Face>(f_id);
        cpt1++;
        if (f.get<Region>().size() == 1)
        {
            m_mesh->mark(f, AMarkBF);
            for (auto e_id:f.getIDs<Edge>())
                m_mesh->mark<Edge>(e_id, AMarkBE);

            for (auto n_id:f.getIDs<Node>())
                m_mesh->mark<Node>(n_id, AMarkBN);

            cpt2++;
        }
    }
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
markCellsOnCurves(int AMarkBF,  //mark for faces on surfaces //IN
                  int AMarkBE,  //mark for edges on surfaces //IN
                  int AMarkCE,  //mark for edges on curves //OUT
                  int AMarkCN)  //mark for nodes on curves //OUT
{

    int cpt1 = 0, cpt2 = 0;

    if(!m_mesh->hasVariable(GMDS_NODE, "BND_CURVE_COLOR") || !m_mesh->hasVariable(GMDS_NODE, "BND_VERTEX_COLOR"))
    {
      for (auto e_id: m_mesh->edges())  {

          if (m_mesh->isMarked<Edge>(e_id, AMarkBE)) {
              Edge e = m_mesh->get<Edge>(e_id);
              cpt1++;
              // We compute normals to the adjacent boundary faces and their
              // scalar products
              std::vector<Face> adj_faces = e.get<Face>();
              if(adj_faces.size()==1) {
                  //2D boundary edge
                  m_mesh->mark(e, AMarkCE);
                  std::vector<Node> e_nodes = e.get<Node>();
                  m_mesh->mark(e_nodes[0], AMarkCN);
                  m_mesh->mark(e_nodes[1], AMarkCN);
              } // if(adj_faces.size()==1) {
              else{
                  std::vector<Face> boundary_adj_faces;

                  for (auto current_face : adj_faces) {
                      if (m_mesh->isMarked(current_face, AMarkBF))
                          boundary_adj_faces.push_back(current_face);
                  }

                  if (boundary_adj_faces.size() != 2){
                      throw GMDSException("a boundary edge should be adjacent to 2 boundary faces!!!");
                  }

                  Face f0 = boundary_adj_faces[0];
                  Face f1 = boundary_adj_faces[1];

                  //LA SUITE DES VECTEURS A CONSTRUIRE A BASE DE POINTS
                  math::Vector3d n0 = /*f0.normal();*/getOutputNormalOfABoundaryFace(f0);
                  math::Vector3d n1 = /*f1.normal();*/getOutputNormalOfABoundaryFace(f1);


                  double dotProduct = n0.dot(n1);

                  if (dotProduct < m_surface_angle_dot)//(sqrt(2.0) / 2.0))
                  {
                      cpt2++;
                      m_mesh->mark(e, AMarkCE);
                      std::vector<TCellID > e_nodes = e.getIDs<Node>();
                      m_mesh->mark<Node>(e_nodes[0], AMarkCN);
                      m_mesh->mark<Node>(e_nodes[1], AMarkCN);
                  }
              }
          }//else{
      } //for (; !it.isDone(); it.next())  {
    }
    else
    {
      Variable<int>* VERTEX_NODE = m_mesh->getVariable<int, GMDS_NODE>("BND_VERTEX_COLOR"  );
      Variable<int>* CURVE_NODE  = m_mesh->getVariable<int, GMDS_NODE>("BND_CURVE_COLOR"  );

      for (auto e_id: m_mesh->edges())  {
        if (m_mesh->isMarked<Edge>(e_id, AMarkBE)) {
            Edge e = m_mesh->get<Edge>(e_id);
            std::vector<Node> e_nodes = e.get<Node>();

            if((*CURVE_NODE)[e_nodes.front().id()] != 0 || (*CURVE_NODE)[e_nodes.back().id()] != 0)
            {
              if(((*CURVE_NODE)[e_nodes.front().id()] == (*CURVE_NODE)[e_nodes.back().id()]) ||
                (((*VERTEX_NODE)[e_nodes.front().id()]) != 0 || (*VERTEX_NODE)[e_nodes.back().id()] != 0))
                {
                    m_mesh->mark(e, AMarkCE);
                    m_mesh->mark(e_nodes.front(), AMarkCN);
                    m_mesh->mark(e_nodes.back(), AMarkCN);
                }
            }
            else if((*VERTEX_NODE)[e_nodes.front().id()] != 0 && (*VERTEX_NODE)[e_nodes.back().id()] != 0)
            {
              m_mesh->mark(e, AMarkCE);
              m_mesh->mark(e_nodes.front(), AMarkCN);
              m_mesh->mark(e_nodes.back(), AMarkCN);
            }
          }
      }
    }
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
markCellsOnCurves(int AMarkCE, //mark for edges on curves //OUT
                  int AMarkCN) //mark for nodes on curves //OUT
{
    int cpt1 = 0, cpt2 = 0;
    for (auto e_id: m_mesh->edges())  {
        Edge e = m_mesh->get<Edge>(e_id);
        cpt1++;
        std::vector<Node> ne = e.get<Node>();
        std::vector<TCellID> adj_faces = m_mesh->getCommonFaces(ne[0],ne[1]);
        if(adj_faces.size()==1) {
            //2D boundary edge
            m_mesh->mark(e, AMarkCE);
            std::vector<TCellID > e_nodes = e.getIDs<Node>();
            m_mesh->mark<Node>(e_nodes[0], AMarkCN);
            m_mesh->mark<Node>(e_nodes[1], AMarkCN);
            cpt2++;
        }
    } //for (; !it.isDone(); it.next())  {
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::markNodesOnPoint(int AMarkCE,// edge on curve IN
                                        int AMarkCN,// node on curve IN
                                        int AMarkPN)// node on vertex OUT
{
    int cpt1=0, cpt2=0;

    if(!m_mesh->hasVariable(GMDS_NODE, "BND_VERTEX_COLOR"))
    {
      for (auto n_id:m_mesh->nodes()) {
          cpt1++;
          if (m_mesh->isMarked<Node>(n_id, AMarkCN)){
              Node n = m_mesh->get<Node>(n_id);
              //We have a node on curve
              std::vector<Edge> adj_edges = n.get<Edge>();
              int cpt_tmp = 0;
              for (const auto& ei:adj_edges){
                  if (m_mesh->isMarked(ei, AMarkCE))
                      cpt_tmp++;
              }
              if (cpt_tmp > 2) {
                  m_mesh->mark(n, AMarkPN);
                  cpt2++;
              }
              else if (cpt_tmp == 2){
                  //check if we have a brutal normal change

                  //First we get all the nodes connected to n by a
                  // boundary edge
                  std::vector<Node> connected_nodes;

                  for (const auto& ei:adj_edges){
                      if (m_mesh->isMarked(ei, AMarkCE)){
                          std::vector<Node> edge_nodes = ei.get<Node>();
                          for (const auto& nj:edge_nodes){
                              if (nj != n)
                                  connected_nodes.push_back(nj);
                          }
                      }
                  }

                  Node n0 = connected_nodes[0];
                  Node n1 = connected_nodes[1];

                  math::Vector3d v0= n0.point()- n.point();
                  math::Vector3d v1= n1.point()- n.point();
                  v0.normalize();
                  v1.normalize();
                  double dotProduct = v0.dot(v1);
                  // if we have a brutal normal change
                  if (dotProduct > -m_surface_angle_dot){
                      m_mesh->mark(n, AMarkPN);
                      cpt2++;
                  }

              }
          }
      }
    }
    else
    {
      Variable<int>* VERTEX_NODE = m_mesh->getVariable<int, GMDS_NODE>("BND_VERTEX_COLOR"  );

      for (auto n_id:m_mesh->nodes()) {
        if((*VERTEX_NODE)[n_id] != 0)
        {
          Node n = m_mesh->get<Node>(n_id);
          m_mesh->mark(n, AMarkPN);
        }
      }
    }
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
colorFaces(int AMarkFOnSurf, int AMarkEOnCurv,
           Variable<int>* AColor)
{
    Variable<int>* var_color = AColor;
    if(var_color==nullptr) {
        try {
            var_color = m_mesh->newVariable<int, GMDS_FACE>("BND_SURFACE_COLOR");
        }
        catch (GMDSException &e) {
            var_color = m_mesh->getVariable<int, GMDS_FACE>("BND_SURFACE_COLOR");
        }
    }

    int color = 0; //Default value is 0
    TInt markDone = m_mesh->newMark<Face>();
    for (auto f_id: m_mesh->faces())
    {
        //on ne considere que les faces au bord
        //qui n'ont pas encore ete traitees
        if ( m_mesh->isMarked<Face>(f_id, AMarkFOnSurf) &&
             !m_mesh->isMarked<Face>(f_id, markDone))
        {
            Face f = m_mesh->get<Face>(f_id);
            //new surface
            color++; // so new color
            m_mesh->mark(f, markDone);
            (*var_color)[f_id] = color;

            //on se propage et on marque
            std::vector<Face> next;
            next.push_back(f);

            while (!next.empty()){
                Face current = next.back();
                next.pop_back();
                //recuperation des faces voisines non traitees et appartenant a la surface
                std::vector<Edge> current_edges = current.get<Edge>();

                for (const auto& ei:current_edges) {
                    if (!m_mesh->isMarked(ei, AMarkEOnCurv))//si ce n'est pas une arete au bord
                    {
                        std::vector<Face> f_edges = ei.get<Face>();
                        for (const auto& fj:f_edges){
                            if (m_mesh->isMarked(fj, AMarkFOnSurf) &&
                                !m_mesh->isMarked(fj, markDone)){
                                m_mesh->mark(fj, markDone);
                                (*var_color)[fj.id()] = color;
                                next.push_back(fj);

                            }
                        }
                    }
                }
            }
        }
    }
    m_mesh->unmarkAll<Face>(markDone);
    m_mesh->freeMark<Face>(markDone);
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
colorEdges(int AMarkEOnCurv, int AMarkNOnPnt,
           Variable<int>* AColor)
{
    Variable<int>* var_color = AColor;
    if(var_color==nullptr) {
        try {
            var_color = m_mesh->newVariable<int, GMDS_EDGE>("BND_CURVE_COLOR");
        }
        catch (GMDSException &e) {
            var_color = m_mesh->getVariable<int, GMDS_EDGE>("BND_CURVE_COLOR");

        }
    }
    int color = 0; //Default value is 0
    TInt markDone = m_mesh->newMark<Edge>();

    std::vector<Edge> done_edges;
    for (auto e_id:m_mesh->edges())
    {
        // We only go throug edges classigied on curves and that have not been
        // yet handled
        if ( m_mesh->isMarked<Edge>(e_id, AMarkEOnCurv) &&
             !m_mesh->isMarked<Edge>(e_id, markDone)) {
            Edge e = m_mesh->get<Edge>(e_id);
            //new curve
            color++; // so new color
            m_mesh->mark(e, markDone);
            (*var_color)[e_id] = color;

            //propagation to curve edges sharing a point with e
            std::vector<Edge> next;
            next.push_back(e);

            while (!next.empty()){
                Edge current = next.back();
                next.pop_back();
                //We get the ajacent edges that are on a curve but not yet done
                std::vector<Node> current_nodes = current.get<Node>();

                for (const auto& ni: current_nodes) {

                    if (!m_mesh->isMarked(ni, AMarkNOnPnt)){
                        //If it is not a node classified on a point, we can found
                        // a next edge on this curve

                        std::vector<Edge> n_edges = ni.get<Edge>();
                        for (const auto& ej:n_edges){
                            if (m_mesh->isMarked(ej, AMarkEOnCurv) &&
                                !m_mesh->isMarked(ej, markDone)){
                                m_mesh->mark(ej, markDone);
                                (*var_color)[ej.id()] = color;
                                next.push_back(ej);

                            }
                        }
                    }
                }
            }
        }
    }
    m_mesh->unmarkAll<Edge>(markDone);
    m_mesh->freeMark<Edge>(markDone);
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
colorNodes(int AMarkNOnPnt, Variable<int>* AColor)
{
    Variable<int>* v_color=AColor;
    if(v_color==nullptr)
    {
        try {
            v_color = m_mesh->newVariable<int, GMDS_NODE>("BND_VERTEX_COLOR");
        }
        catch (GMDSException &e) {
            v_color = m_mesh->getVariable<int, GMDS_NODE>("BND_VERTEX_COLOR");

        }
    }

    int color = 0; //Default value is 0
    for (auto n_id : m_mesh->nodes() )
    {
        if (m_mesh->isMarked<Node>(n_id, AMarkNOnPnt)) {
            //new point
            color++; // so new color
            (*v_color)[n_id] = color;
        }

    }
}
/*----------------------------------------------------------------------------*/
math::Vector3d BoundaryOperator::getOutputNormal(Face& AFace, Region& ARegion)
{
    std::vector<Node> region_nodes = ARegion.get<Node>();
    std::vector<Node> face_nodes = AFace.get<Node>();

    math::Point region_center = ARegion.center();
    math::Point face_center = AFace.center();
    math::Vector3d face_normal = AFace.normal();

    math::Vector3d inward= region_center-face_center;
    if (inward.dot(face_normal)>0.0)
    {
        return math::Vector3d({
                -face_normal.X(),
                -face_normal.Y(),
                -face_normal.Z()});
    }
    return face_normal;
}
/*----------------------------------------------------------------------------*/
math::Vector3d BoundaryOperator::
getOutputNormalOfABoundaryFace(Face& AFace)
{
    std::vector<Region> adj_regions = AFace.get<Region>();
    if (adj_regions.size() != 1)
        throw GMDSException("A boundary face must be adjacent to 1 region!!!");

    return getOutputNormal(AFace, adj_regions[0]);
}
/*----------------------------------------------------------------------------*/
math::Vector3d BoundaryOperator::
getOutputNormalOfABoundaryNode(const Node& ANode)
{
    std::vector<Face> adj_faces = ANode.get<Face>();
    std::vector<math::Vector3d> weighted_normals;
    for(auto fi:adj_faces){
        std::vector<Region> fi_regions = fi.get<Region>();
        if (fi_regions.size() == 1){
            // we have a boundary face
            math::Vector3d vi = getOutputNormal(fi, fi_regions[0]);
            TCoord ai = fi.area();
            weighted_normals.push_back(ai*vi);

        }
    }
    // now we compute the normal vector at ANode
    math::Vector3d n = weighted_normals[0];
    for(auto i=1; i <weighted_normals.size();i++){
        n = n+weighted_normals[i];
    }
    n.normalize();
    return n;
}
/*----------------------------------------------------------------------------*/
