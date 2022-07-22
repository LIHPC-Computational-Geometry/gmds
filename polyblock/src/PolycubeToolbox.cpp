/*--------------------------------------------------------------------------*/
#include <gmds/polyblock/PolycubeToolbox.h>
/*--------------------------------------------------------------------------*/
#include <gmds/math/Vector.h>
#include <gmds/ig/Mesh.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/io/VTKWriter.h>
/*--------------------------------------------------------------------------*/
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <ctime>
#include <fstream>
#include <gmds/io/IGMeshIOService.h>
/*--------------------------------------------------------------------------*/
#include "GCoptimization.h"
/*--------------------------------------------------------------------------*/
using namespace gmds;
using namespace Eigen;
/*--------------------------------------------------------------------------*/
typedef SparseMatrix<double> SpMat;
/*--------------------------------------------------------------------------*/
// MATHMECA
//  efficient GMRES for Ax = b sparse-solving
/**This code is cortesy of V. Delmas, implementing the algorithm suggested by
Ayachour, E. (2003). A fast implementation for gmres method. */
void gmresAya(SparseMatrix<double>& A, SparseVector<double>& x,
              SparseVector<double>& b, double tol=1E-13,
              int maxIter=1000, int m_max=80, int precond=1){
	assert(A.cols()==A.rows());
		const int n = A.cols();
		assert(m_max>2);


		std::ofstream residus("aya.dat");

		SparseMatrix<double> V(n,m_max+1);
		SparseMatrix<double> H(m_max+1,m_max);

		SparseVector<double> y(m_max+1);
		SparseVector<double> r(n);
		SparseVector<double> u(m_max);

		double beta;
		SparseMatrix<double> rot(m_max+1,2);

	    assert(n>3);
	    int d(10), m_min(20), m(m_max);
	    double c(1.), max_c(cos(8*3.14159265/180.)), min_c(cos(80*3.14159265/180.));

		int nbIter(0), mloc(m);
		do{
	        if(c>max_c){m=m_max;}
	        else if(c<min_c){m=m_min;}
	        else{
	            if((m-d)>m_min){m=m-d;}
	            else{m=m_max/2;}
	        }
	        //cout << m << endl;
			mloc = m;
			r=b-A*x;

			if(precond==1){for(int k=0;k<n;k++){r.coeffRef(k)/=A.coeffRef(k,k);}}
			residus << nbIter << " " << beta << " " << r.norm() << endl;
			beta=r.norm();

			V.col(0) = r/beta;

			for(int j=0;j<m;j++){
				SparseVector<double> Av = A*V.col(j);
				if(precond==1){for(int k=0;k<n;k++){Av.coeffRef(k)/=A.coeffRef(k,k);}}
				V.col(j+1) = Av;
				for(int i=0;i<=j;i++){
					H.coeffRef(i,j) = Av.dot(V.col(i));
					V.col(j+1) -= H.coeffRef(i,j)*V.col(i);
				}
				H.coeffRef(j+1,j) = V.col(j+1).norm();
				if(H.coeffRef(j+1,j)<1e-13){mloc=j;break;}
				V.col(j+1) /= H.coeffRef(j+1,j);

				u.coeffRef(j) = H.coeffRef(0,j);
				for(int p=0;p<j;p++){u.coeffRef(j) -= H.coeffRef(p+1,j)*u.coeffRef(p);}
				u.coeffRef(j) /= H.coeffRef(j+1,j);
				c = 1./sqrt(1.+u.head(m).dot(u.head(m)));
				//cout << u << endl << c << endl;
				nbIter++;

				if(beta*c<tol){
					mloc = j+1;
					break;
				}

			}
			beta *=c;
			//cout << c << endl;
			//Backsolve de Hk*y = beta*c*u
			for(int k=mloc-1;k>=0;k--){
				y.coeffRef(k) = beta*c*u.coeffRef(k);
				for(int p=k+1;p<mloc;p++){y.coeffRef(k) -= H.coeffRef(k+1,p)*y.coeffRef(p);}
				y.coeffRef(k) /= H.coeffRef(k+1,k);
				x+= V.col(k)*y.coeffRef(k);
			}


	}while((beta>tol)and(nbIter<maxIter));

}
// END MATHMECA

Eigen::Vector3d math2eigen(math::Vector v)
{
		return Eigen::Vector3d(v.X(),v.Y(),v.Z());
}


/*----------------------------------------------------------------------------*/
PolycubeToolbox::PolycubeToolbox(Mesh *mesh)
{
	this->mesh = mesh;

	mesh->deleteVariable(GMDS_FACE, "normal");
	mesh->deleteVariable(GMDS_FACE, "ChosenNormal_ID");
	mesh->deleteVariable(GMDS_FACE, "ChosenNormal");
	mesh->deleteVariable(GMDS_FACE, "chart");
	mesh->deleteVariable(GMDS_FACE, "constrained");

	face_normals        = mesh->newVariable<math::Vector,GMDS_FACE>("normal");
	face_chosenNormal_ID= mesh->newVariable<int         ,GMDS_FACE>("ChosenNormal_ID");
	face_chosenNormal   = mesh->newVariable<math::Vector,GMDS_FACE>("ChosenNormal");
	face_chart          = mesh->newVariable<int         ,GMDS_FACE>("chart");
	face_constrained    = mesh->newVariable<int         ,GMDS_FACE>("constrained");



	mesh->deleteVariable(GMDS_NODE, "normal");
	mesh->deleteVariable(GMDS_NODE, "ChosenNormal");
	mesh->deleteVariable(GMDS_NODE, "rotationMatrix");
	mesh->deleteVariable(GMDS_NODE, "Quaternion");

	point_normals       = mesh->newVariable<math::Vector      ,GMDS_NODE>("normal");
	point_chosenNormal  = mesh->newVariable<math::Vector      ,GMDS_NODE>("ChosenNormal");
	point_rotationMatrix= mesh->newVariable<Eigen::Matrix3d          ,GMDS_NODE>("rotationMatrix");
	point_quaternion    = mesh->newVariable<Quaternion<double>,GMDS_NODE>("Quaternion");

	mesh->deleteVariable(GMDS_NODE, "updatedX");
	mesh->deleteVariable(GMDS_NODE, "updatedY");
	mesh->deleteVariable(GMDS_NODE, "updatedZ");

	node_updatedX = mesh->newVariable<int,GMDS_NODE>("updatedX");
	node_updatedY = mesh->newVariable<int,GMDS_NODE>("updatedY");
	node_updatedZ = mesh->newVariable<int,GMDS_NODE>("updatedZ");

	mesh->deleteVariable(GMDS_NODE, "chartCorner");
	mesh->deleteVariable(GMDS_NODE, "turningPoint");

	node_chartCorner  = mesh->newVariable<int,GMDS_NODE>("chartCorner");
	node_turningPoint = mesh->newVariable<int,GMDS_NODE>("turningPoint");


	mesh->deleteVariable(GMDS_EDGE, "normal");

	edge_normal = mesh->newVariable<math::Vector,GMDS_EDGE>("normal");


	mesh->deleteVariable(GMDS_NODE, "browsed");

	node_browsed = mesh->newVariable<int,GMDS_NODE>("browsed");


	double minX,minY,minZ;
	double maxX,maxY,maxZ;
	minX = minY = minZ = 1E80;
	maxX = maxY = maxX = -1E80;

	for (auto n_id: mesh->nodes()) {
		Node node = mesh->get<Node>(n_id);
		if (node.X()>maxX){
			maxX= node.X();
		}
		if (node.Y()>maxY){
			maxY= node.Y();
		}
		if (node.Z()>maxZ){
			maxZ= node.Z();
		}
		if (node.X()<minX){
			minX= node.X();
		}
		if (node.Y()<minY){
			minY= node.Y();
		}
		if (node.Z()<minZ){
			minZ= node.Z();
		}
	}

	BOUNDING_BOX_LENTH = sqrt((maxX-minX)*(maxX-minX) + (maxY-minY)*(maxY-minY) + (maxZ-minZ)*(maxZ-minZ));
	std::cout << "The mesh has a bounding box with a diagonal of " << BOUNDING_BOX_LENTH << '\n';

	ANGLE_DETERMINENT = 0.05;
	GREGSON2011_IS_INIT = false;

	markOnBoundary();
}
/*----------------------------------------------------------------------------*/
void PolycubeToolbox::markOnBoundary()
{
	int mark_edge_on_curv = mesh->newMark<Edge>();

	int mark_node_on_surf = mesh->newMark<Node>();
	int mark_node_on_curv = mesh->newMark<Node>();
	int mark_node_on_pnt  = mesh->newMark<Node>();
	int mark_node_alone   = mesh->newMark<Node>();

	mark_smooth_bnd_nodes = mesh->newMark<Node>();
	mark_bnd_nodes = mesh->newMark<Node>();
	mark_bnd_faces = mesh->newMark<Face>();
	mark_sharp_edges = mesh->newMark<Edge>();
    mark_bnd_edges = mesh->newMark<Edge>();

    BoundaryOperator bnd_op(mesh);
    mesh->deleteVariable(GMDS_FACE, "BND_SURFACE_COLOR");

    bnd_op.markCellOnGeometry(mark_bnd_faces,
                              mark_bnd_edges,
                              mark_node_on_surf,
                              mark_sharp_edges,
                              mark_node_on_curv,
                              mark_node_on_pnt,
                              mark_node_alone);

	for (auto n_id:mesh->nodes()) {
		Node n = mesh->get<Node>(n_id);
		if(mesh->isMarked(n,mark_node_on_surf)){
			mesh->mark(n,mark_bnd_nodes);
			if(!mesh->isMarked(n,mark_node_on_curv) &&
				 !mesh->isMarked(n,mark_node_on_pnt) ){
				mesh->mark(n,mark_smooth_bnd_nodes);
			}
		}
	}

	mesh->unmarkAll<Node>(mark_node_on_surf);
	mesh->unmarkAll<Node>(mark_node_on_curv);
	mesh->unmarkAll<Node>(mark_node_on_pnt);
	mesh->unmarkAll<Node>(mark_node_alone);
	mesh->freeMark<Node>(mark_node_on_surf);
	mesh->freeMark<Node>(mark_node_on_curv);
	mesh->freeMark<Node>(mark_node_on_pnt);
	mesh->freeMark<Node>(mark_node_alone);

	mesh->unmarkAll<Edge>(mark_edge_on_curv);
	mesh->freeMark<Edge>(mark_edge_on_curv);
}
/*----------------------------------------------------------------------------*/
// fully updating the boundary, might lead to non convex geometry
void PolycubeToolbox::updateBoundary(){
	int mark_edge_on_surf = mesh->newMark<Edge>();
	int mark_edge_on_curv = mesh->newMark<Edge>();

	int mark_node_on_surf = mesh->newMark<Node>();
	int mark_node_on_curv = mesh->newMark<Node>();
	int mark_node_on_pnt  = mesh->newMark<Node>();
	int mark_node_alone   = mesh->newMark<Node>();


	BoundaryOperator bnd_op(mesh);
    mesh->deleteVariable(GMDS_FACE, "BND_SURFACE_COLOR");

    bnd_op.markCellOnGeometry(mark_bnd_faces,
                              mark_edge_on_surf,
                              mark_node_on_surf,
                              mark_sharp_edges,
                              mark_node_on_curv,
                              mark_node_on_pnt,
                              mark_node_alone);

    for (auto n_id : mesh->nodes()) {
        Node n = mesh->get<Node>(n_id);
		if(mesh->isMarked(n,mark_node_on_surf)){
			mesh->mark(n,mark_bnd_nodes);
			if(!mesh->isMarked(n,mark_node_on_curv) &&
				 !mesh->isMarked(n,mark_node_on_pnt) ){
				mesh->mark(n,mark_smooth_bnd_nodes);
			}
		}
	}

	mesh->unmarkAll<Node>(mark_node_on_surf);
	mesh->unmarkAll<Node>(mark_node_on_curv);
	mesh->unmarkAll<Node>(mark_node_on_pnt);
	mesh->unmarkAll<Node>(mark_node_alone);
	mesh->freeMark<Node>(mark_node_on_surf);
	mesh->freeMark<Node>(mark_node_on_curv);
	mesh->freeMark<Node>(mark_node_on_pnt);
	mesh->freeMark<Node>(mark_node_alone);

	mesh->unmarkAll<Edge>(mark_edge_on_surf);
	mesh->unmarkAll<Edge>(mark_edge_on_curv);
	mesh->freeMark<Edge>(mark_edge_on_surf);
	mesh->freeMark<Edge>(mark_edge_on_curv);
}
/*----------------------------------------------------------------------------*/
// compute poisson with n on diag and -1 where link to other node
SpMat PolycubeToolbox::computePoisson(){
	SpMat poisson_matrix(mesh->getNbNodes(),mesh->getNbNodes());

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(2*mesh->getNbEdges()+mesh->getNbNodes());

	std::vector<TCellID> nodes;

    for (auto e_id : mesh->edges()) {
        Edge e = mesh->get<Edge>(e_id);	
        nodes = e.getIDs<Node>();
		if(nodes[0] != 0){
			tripletList.push_back(T(nodes[0],nodes[1],-1));
		}
		if(nodes[1] != 0){
			tripletList.push_back(T(nodes[1],nodes[0],-1));
		}
	}


    for (auto n_id : mesh->nodes()) {
        Node n = mesh->get<Node>(n_id);
        if(n.id()!= 0){
            tripletList.push_back(T(n.id(),n.id(),n.nbEdges()));
        }
        else{
            tripletList.push_back(T(n.id(),n.id(),1));
        }
    }

	poisson_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
	poisson_matrix.makeCompressed();

	return poisson_matrix;
}
/*----------------------------------------------------------------------------*/
// compute adjacency with -n on diagonal where angle = 0, and diag 1 otherwise
SpMat PolycubeToolbox::computeAdjacency(){
	SpMat adjacency_matrix(mesh->getNbNodes(),mesh->getNbNodes());

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(2*mesh->getNbEdges()+mesh->getNbNodes());

	std::vector<TCellID> nodes;

    for (auto e_id : mesh->edges()) {
        Edge e = mesh->get<Edge>(e_id);
        nodes = e.getIDs<Node>();
		if(!mesh->isMarked(mesh->get<Node>(nodes[0]),mark_smooth_bnd_nodes)){
			tripletList.push_back(T(nodes[0],nodes[1],1));
		}
		if(!mesh->isMarked(mesh->get<Node>(nodes[1]),mark_smooth_bnd_nodes)){
			tripletList.push_back(T(nodes[1],nodes[0],1));
		}
	}


    for (auto n_id : mesh->nodes()) {
        Node n = mesh->get<Node>(n_id);
        if(!mesh->isMarked(n,mark_smooth_bnd_nodes)){
			tripletList.push_back(T(n.id(),n.id(),-n.nbEdges()));
		}
		else{
			tripletList.push_back(T(n.id(),n.id(),1));
		}
	}
	adjacency_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
	adjacency_matrix.makeCompressed();
	return adjacency_matrix;
}
/*----------------------------------------------------------------------------*/
// compute adjacency with -n on diagonal
SpMat PolycubeToolbox::computeLaplace(){
	SpMat adjacency_matrix(mesh->getNbNodes(),mesh->getNbNodes());

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(2*mesh->getNbEdges()+mesh->getNbNodes());

	std::vector<TCellID> nodes;

    for (auto e_id : mesh->edges()) {
        Edge e = mesh->get<Edge>(e_id);
        nodes = e.getIDs<Node>();
		if(!mesh->isMarked(mesh->get<Node>(nodes[0]),mark_bnd_nodes)){
			tripletList.push_back(T(nodes[0],nodes[1],-1));
		}
		if(!mesh->isMarked(mesh->get<Node>(nodes[1]),mark_bnd_nodes)){
			tripletList.push_back(T(nodes[1],nodes[0],-1));
		}
	}


    for (auto n_id : mesh->nodes()) {
        Node n = mesh->get<Node>(n_id);
		if(!mesh->isMarked(n,mark_bnd_nodes)){
			tripletList.push_back(T(n.id(),n.id(),n.nbEdges()));
		}
		else{
			tripletList.push_back(T(n.id(),n.id(),1));
		}
	}
	adjacency_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
	adjacency_matrix.makeCompressed();
	return adjacency_matrix;
}
/*----------------------------------------------------------------------------*/
void PolycubeToolbox::graphCutBasedFaceSegmentation(){
  //First we take all the faces we are going to work on.
	std::cout << "==giving a first labelling via graph cut==" << '\n';

	int COMPTEUR=0;
	linker_ID2bnd_faces.clear();
	bnd_faces.clear();

	std::cout << "generating graph ..." << '\n';
	int num=0;
	for(auto f_id : mesh->faces()){
		Face f = mesh->get<Face>(f_id);
		if (mesh->isMarked(f,mark_bnd_faces)) {
			bnd_faces.push_back(f.id());
			linker_ID2bnd_faces.insert( std::pair<TCellID,int>(f.id(),num) ); num++;
		}
	}

	int nb_faces = bnd_faces.size();
	gc_face_is_constrained.reserve(nb_faces*6);
	for (size_t i = 0; i < nb_faces*6; i++) {
		gc_face_is_constrained[i]=false;
	}

	std::map<std::pair<TCellID,int>,double> bnd_faces_gc_face_is_constrained; // useful to know to which id was assigned a face.

	bool OverloadedCorner = true, NonPlannarCharts = true;
	while ((OverloadedCorner || NonPlannarCharts) && COMPTEUR<8 )
	{
		GCoptimizationGeneralGraph graphCut(nb_faces,6);

		for (size_t i = 0; i < nb_faces; i++) {
			std::vector<TCellID> edges = ( mesh->get<Face>(bnd_faces[i]) ).getIDs<Edge>();
			for (auto e : edges) {
				std::vector<TCellID> faces = ( mesh->get<Edge>(e) ).getIDs<Face>();
				for (auto f : faces) {
					if (mesh->isMarked<Face>(f,mark_bnd_faces) && bnd_faces[i] != f) {
						graphCut.setNeighbors(i, linker_ID2bnd_faces[f]);
					}
				}
			}
		}

		// I use an array  to avoid static fonction pb
		std::vector<double> dataCost;
		for (auto f : bnd_faces){
			for (size_t i = 0; i < 6; i++) {
				double nts=  (( mesh->get<Face>(f) ).normal().dot( vect[i] ) - 1) /0.2;
				dataCost.push_back(1-exp(-0.5*nts*nts));
			}
		}
		// take constraint into account;
		for (auto constraint : gc_face_constraint) {
			std::pair<TCellID,int> pos = constraint.first;
			dataCost[linker_ID2bnd_faces[pos.first]*6 + pos.second] = constraint.second;
		}
		graphCut.setDataCost(dataCost.data());

		std::vector<math::Vector> extraData;
		for (size_t i = 0; i < nb_faces; i++) {
			extraData.push_back( (mesh->get<Face>(bnd_faces[i])).normal() );
		}

		graphCut.setSmoothCost(&PolycubeToolbox::graphCutCpq, (void *)&extraData);

		std::cout << "Running graphCut algorithm ...";
		std::cout.flush();
		const clock_t map_time = clock();
        try {
            graphCut.expansion(300);
            std::cout << "energy :" << graphCut.swap(1000) <<'\n';
        } catch (const GCException& e) {std::cout << e.message << '\n'; }

        std::cout << "Computing time : " << float( clock () - map_time ) /  CLOCKS_PER_SEC << "sec" << endl;


        for(auto f_id : mesh->faces()){

            Face f = mesh->get<Face>(f_id);
            if (mesh->isMarked(f,mark_bnd_faces)) {
                (*face_normals)[f.id()] = f.normal();
                (*face_chosenNormal_ID)[f.id()] = graphCut.whatLabel(linker_ID2bnd_faces[f.id()]);
                (*face_chosenNormal)[f.id()] = vect[(*face_chosenNormal_ID)[f.id()]];
            } else {
                (*face_normals)[f.id()] = f.normal();
                (*face_chosenNormal_ID)[f.id()] = -1;
                (*face_chosenNormal)[f.id()] = vect[(*face_chosenNormal_ID)[f.id()]];
            }
        }
    }
    this->chartsGeneration();
    NonPlannarCharts = this->checkingChartDivision();
    OverloadedCorner = this->checkingChartCornerOverloading();



    this->compute_turning_point();
    std::cout << "Writing results ..." << '\n';
    IGMeshIOService ioService(mesh);
    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(N|F);
    vtkWriter.setDataOptions(N| F);
    vtkWriter.write("result_check"  + std::to_string(COMPTEUR++));

    if (OverloadedCorner&& NonPlannarCharts) {
        std::cout << "Non plannar chart division and overloaded chart corner detected. Restarting labelling with more constraint :" << '\n';
    } else if (OverloadedCorner){
        std::cout << "Overloaded chart corner detected. Restarting labelling with more constraint :" << '\n';
    } else if (NonPlannarCharts){
        std::cout << "Non plannar chart division detected. Restarting labelling with more constraint :" << '\n';
    } else {
        std::cout << "Valid first segmentation." << '\n';
    }
}


/*----------------------------------------------------------------------------*/
bool PolycubeToolbox::checkingChartCornerOverloading(){
	std::cout << "Looking for overloaded chart's corner ..." << '\n';
	/*we look for all the corner with a valence of more than 3, and apply a constraint
	on them */
	// we compute the valence of each corner
	bool OverloadedCorner = false;
	map<TCellID,int> corner_valence;
	for (auto e : charts_edges){
		if ( corner_valence.find(e.corner1) == corner_valence.end() ){
			corner_valence.insert( std::pair<TCellID,int>(e.corner1,1) );
		}
		else {
			corner_valence[e.corner1]++;
		}
		if ( corner_valence.find(e.corner2) == corner_valence.end() ){
			corner_valence.insert(std::pair<TCellID,int>(e.corner2,1));
		}
		else {
			corner_valence[e.corner2]++;
		}
	}

	// now for each corner with a valence >3 we will apply constaint
	for (auto& corner : corner_valence){
		// corner.second/=2;

		if(corner.second>3) {
			OverloadedCorner = true;

			// separate the faces between detected mesh edge
			int marked_face = mesh->newMark<Face>();

			std::vector<TCellID>  faces = ( mesh->get<Node>(corner.first) ).getIDs<Face>();
			for (auto f : faces) {
				if (mesh->isMarked<Face>(f,mark_bnd_faces)){
					mesh->mark<Face>(f, marked_face);
				}
			}

			int marked_edge = mesh->newMark<Edge>();
			std::vector<TCellID> edges = (mesh->get<Node>(corner.first)).getIDs<Edge>();
			TCellID current_edge;
			for (auto e : edges){
				if (mesh->isMarked<Edge>(e,mark_bnd_edges)){
					mesh->mark<Edge>(e, marked_edge);
					current_edge = e;
				}
			}
			TCellID first_edge = current_edge;

			std::vector<std::vector<TCellID>> different_faces_set;
			different_faces_set.push_back(std::vector<TCellID>());
			bool continu = true;
			while ( continu ){
				faces.clear(); faces = (mesh->get<Edge>(current_edge)).getIDs<Face>();
				mesh->unmark<Edge>(current_edge,marked_edge);

				TCellID current_face;
				for (auto f : faces){
					if (mesh->isMarked<Face>(f,marked_face)){
						current_face = f;
					}
				}
				mesh->unmark<Face>(current_face,marked_face);
				different_faces_set.back().push_back(current_face);
				continu=false;
				edges.clear(); edges = (mesh->get<Face>(current_face)).getIDs<Edge>();
				for (auto e : edges){
					if (mesh->isMarked<Edge>(e,marked_edge)){
						current_edge = e; continu= true;
					}
				}
				if (mesh->isMarked<Edge>(current_edge, mark_sharp_edges) && continu){
					different_faces_set.push_back(std::vector<TCellID>());
				}
			}
			//we check if we have to merge the first and last part
			if (!mesh->isMarked<Edge>(first_edge, mark_sharp_edges) &&  different_faces_set.size() > 1){
				for (auto f : different_faces_set.back()){
					different_faces_set.front().push_back(f);
				}
				different_faces_set.pop_back();
			}

			// if we have a node with just a number of sharp edge to high :
			if (different_faces_set.size() >= corner.second){
				OverloadedCorner =false;
			}
			mesh->freeMark<Face>(marked_face);
			mesh->freeMark<Edge>(marked_edge);
			for (auto on_bnd_faces : different_faces_set) {
				double area[6] = {0,0,0,0,0,0};
				// we check if there is some constraint faces, we want to keep it the same
				for (auto f : on_bnd_faces){
					for (size_t i = 0; i < 6; i++) {
						if (gc_face_is_constrained[linker_ID2bnd_faces[f]*6+i] ){
							// we set node constraint at 10 and chart division at 50, so we prioritise chart division, to avoid looping
							area[i] -= gc_face_constraint[std::pair<TCellID,int>(f,i)];
						}
					}
				}

				//if not, we take the label with the maximum area
				for (auto f : on_bnd_faces){
					area[(*face_chosenNormal_ID)[f]] += (mesh->get<Face>(f)).area();
				}
				double max = area[0]; int forced_labbelling = 0;
				for (size_t i = 1; i < 6; i++) {
					if (max < area[i]){
						max = area[i];
						forced_labbelling = i;
					}
				}
				// we now make the constraint
				for (auto f : on_bnd_faces){
					(*face_constrained)[f]=1;
					for (size_t i = 0; i < 6; i++) {
						if (i == forced_labbelling){
							//new constraint
							if (gc_face_is_constrained[linker_ID2bnd_faces[f]*6+i]) {
								gc_face_constraint[std::pair<TCellID,int>(f,i)] = 0;
							} else {
								gc_face_constraint.insert( std::pair<std::pair<TCellID,int>,double>(std::pair<TCellID,int>(f,i) ,0.) );
							}
						} else {
							//new constraint
							if (gc_face_is_constrained[linker_ID2bnd_faces[f]*6+i]) {
								gc_face_constraint[std::pair<TCellID,int>(f,i)] = 10;
							} else {
								gc_face_constraint.insert( std::pair<std::pair<TCellID,int>,double>(std::pair<TCellID,int>(f,i) ,10.) );
							}
						}
						gc_face_is_constrained[linker_ID2bnd_faces[f]*6+i] = true;

					}
				}
			}


		}
	}
	return OverloadedCorner;
}

/*----------------------------------------------------------------------------*/
bool PolycubeToolbox::checkingChartDivision(){
	bool NonPlannarCharts = false;
	std::cout << " checking chart division " << '\n';
	/* this part is equivalent to checking if two adjacent chart have the same
	label, if we were to compute a polycube, it would give bad result. we will
	therefore constrain on side of the edge to another value*/
	for (auto chart_edge : charts_edges){
		if ( (*face_chosenNormal_ID)[ (charts[chart_edge.chart1])[0] ]
		 ==  (*face_chosenNormal_ID)[ (charts[chart_edge.chart2])[0] ]) {
			 int chart_to_modify;
			 std::vector<TCellID> chart1_close_faces,chart2_close_faces;
			 // we take the smallest chart to enforce value on its edge
			 for (auto edge : chart_edge.edge){
				 std::vector<TCellID> nodes = (mesh->get<Edge>(edge)).getIDs<Node>();
				 for (auto n : nodes){
					 std::vector<TCellID> wanted_faces,faces = (mesh->get<Node>(n)).getIDs<Face>();
					 for (auto f: faces){
						 if (mesh->isMarked<Face>(f,mark_bnd_faces)){
							 wanted_faces.push_back(f);
						 }
					 }
					 for (auto f : wanted_faces){
						 if ((*face_chart)[f] == chart_edge.chart1) {
							 chart1_close_faces.push_back(f);
						 } else if ((*face_chart)[f] == chart_edge.chart2){
							 chart2_close_faces.push_back(f);
						 }
					 }
				 }
			 }

			 math::Vector averageVector1({0., 0., 0.}),averageVector2({0., 0., 0.});
			 for (auto f : chart1_close_faces) {
				 averageVector1 = averageVector1 + (*face_normals)[f];
			 }
			 averageVector1 = averageVector1 / chart1_close_faces.size();
			 for (auto f : chart2_close_faces) {
				 averageVector2 = averageVector2 + (*face_normals)[f];
			 }
			 averageVector2 = averageVector2 / chart2_close_faces.size();
			 // we set constraint to the set of face which are the furthest from their normal,
			 // on average
			 if (averageVector2.dot((*face_chosenNormal)[ (charts[chart_edge.chart1])[0] ])
				> averageVector1.dot((*face_chosenNormal)[ (charts[chart_edge.chart1])[0] ])){
					for (auto f : chart1_close_faces){
						(*face_constrained)[f] = 1;
						//new constraint
						if (gc_face_is_constrained[linker_ID2bnd_faces[f]*6+(*face_chosenNormal_ID)[f]]) {
							gc_face_constraint[std::pair<TCellID,int>(f,(*face_chosenNormal_ID)[f])] = 50;
						} else {
							gc_face_constraint.insert( std::pair<std::pair<TCellID,int>,double>(std::pair<TCellID,int>(f,(*face_chosenNormal_ID)[f]) ,50.) );
						}
						gc_face_is_constrained[linker_ID2bnd_faces[f]*6+(*face_chosenNormal_ID)[f]] = true;
						NonPlannarCharts =true;

					}
				}else {
					for (auto f : chart2_close_faces){
						(*face_constrained)[f] = 1;

						//new constraint
						if (gc_face_is_constrained[linker_ID2bnd_faces[f]*6+(*face_chosenNormal_ID)[f]]) {
							gc_face_constraint[std::pair<TCellID,int>(f,(*face_chosenNormal_ID)[f])] = 50;
						} else {
							gc_face_constraint.insert( std::pair<std::pair<TCellID,int>,double>(std::pair<TCellID,int>(f,(*face_chosenNormal_ID)[f]) ,50.) );
						}
						gc_face_is_constrained[linker_ID2bnd_faces[f]*6+(*face_chosenNormal_ID)[f]] = true;
						NonPlannarCharts =true;
					}
				}
		}
	}
	return NonPlannarCharts;
}
/*----------------------------------------------------------------------------*/
/* we take advantage of the fact that we only need the vect and label  to compute
	the cost, not needing any normal */
double PolycubeToolbox::graphCutCpq(int p, int q, int label_p, int label_q, void *extraData)
{
	double sigma = 0.25;
	double cons =1;
	if (label_q == label_p){
		return 0;
	} else if ((vect[label_p]).dot( vect[label_q]) < -0.5) { // if the two labels are opposites
		return 50;
	} else {
		double npnq = ( ((*(std::vector<math::Vector> *)(extraData))[p]).dot( (*(std::vector<math::Vector> *)(extraData))[q] ) -1)/sigma;
		return cons*exp(-0.5*npnq*npnq);
	}

}
/*----------------------------------------------------------------------------*/
double cpqTurningPoint(int p, int q, int label_p, int label_q, void *extraData){
	double sigma = 1.2;
	double cons =3;
	if (label_q == label_p){
		return 0;
	} else{
		double npnq = ( ((*(std::vector<math::Vector> *)(extraData))[p]).dot( (*(std::vector<math::Vector> *)(extraData))[q] ) -1)/sigma;
		return cons*exp(-0.5*npnq*npnq);
	}
}
/*----------------------------------------------------------------------------*/
void PolycubeToolbox::compute_turning_point(bool speaking){
	if (charts_edges.empty()){
        this->chartsGeneration();
	}
	if (speaking)
	std::cout << "computing turning points ..." << '\n';
	for (auto n_id : mesh->nodes()) {
		(*node_turningPoint)[n_id] = 0;
	}

	int nb_turning_point = 0;
	int nb_unary_edge = 0;


	for (auto& TEDGE : charts_edges){
		if (TEDGE.edge.size()>1) {

			int supposed_axis_of_edge;
			int axis1 = (*face_chosenNormal_ID)[charts[TEDGE.chart1].back()];
			int axis2 = (*face_chosenNormal_ID)[charts[TEDGE.chart2].back()];

			if ((axis1 == 0 || axis1 == 1) && (axis2 ==2 || axis2 ==3)) {supposed_axis_of_edge=4; }
			if ((axis1 == 0 || axis1 == 1) && (axis2 ==4 || axis2 ==5)) {supposed_axis_of_edge=2; }
			if ((axis1 == 2 || axis1 == 3) && (axis2 ==4 || axis2 ==5)) {supposed_axis_of_edge=0; }

			if ((axis2 == 0 || axis2 == 1) && (axis1 ==2 || axis1 ==3)) {supposed_axis_of_edge=4; }
			if ((axis2 == 0 || axis2 == 1) && (axis1 ==4 || axis1 ==5)) {supposed_axis_of_edge=2; }
			if ((axis2 == 2 || axis2 == 3) && (axis1 ==4 || axis1 ==5)) {supposed_axis_of_edge=0; }


			int nb_edge = TEDGE.edge.size();


			GCoptimizationGeneralGraph graphCut(nb_edge,2);
			for (size_t i = 0; i < nb_edge-1; i++) {
				graphCut.setNeighbors(i,i+1);
			}
			/* if we have a chart isolated into another chart
			we have looping edge*/
			if (TEDGE.corner1 == TEDGE.corner2){
				graphCut.setNeighbors(0,nb_edge-1);
			}
			std::vector<double> dataCost;
			for (auto& e : TEDGE.edge){
				math::Vector normal = (*edge_normal)[e];

				double nts=  (normal.dot( vect[supposed_axis_of_edge] ) - 1) /0.9;
				// std::cout << normal << '\n';
				// std::cout << nts << '\n';
				dataCost.push_back(1-exp(-0.5*nts*nts));
				nts=  (normal.dot( vect[supposed_axis_of_edge+1] ) - 1) /0.9;
				// std::cout << nts << '\n';

				dataCost.push_back(1-exp(-0.5*nts*nts));

			}

			graphCut.setDataCost(dataCost.data());


			std::vector<math::Vector> extraData;
			for (auto& e : TEDGE.edge){
				extraData.push_back( (*edge_normal)[e] );
			}

			graphCut.setSmoothCost(&cpqTurningPoint, (void *)&extraData);


			try { //graphCut.expansion(300);
				graphCut.swap(300);
			} catch (const GCException& e) {std::cout << e.message << '\n'; }
			int precedent_label = graphCut.whatLabel(0);
			for (size_t i = 1; i < nb_edge; i++) {
				if (graphCut.whatLabel(i) != precedent_label){
					std::vector<TCellID> nodes1 = (mesh->get<Edge>( TEDGE.edge[i])).getIDs<Node>();
					std::vector<TCellID> nodes2 = (mesh->get<Edge>( TEDGE.edge[i-1])).getIDs<Node>();
					if (nodes1[0] == nodes2[0]) { (*node_turningPoint)[nodes1[0]]=1;
						TEDGE.turningPoints.push_back(nodes1[0]); nb_turning_point++;} else
					if (nodes1[0] == nodes2[1]) { (*node_turningPoint)[nodes1[0]]=1;
						TEDGE.turningPoints.push_back(nodes1[0]); nb_turning_point++;} else
					if (nodes1[1] == nodes2[0]) { (*node_turningPoint)[nodes1[1]]=1;
						TEDGE.turningPoints.push_back(nodes1[1]); nb_turning_point++;} else
					if (nodes1[1] == nodes2[1]) { (*node_turningPoint)[nodes1[1]]=1;
						TEDGE.turningPoints.push_back(nodes1[1]); nb_turning_point++;}

					// (*node_browsed)[nodes1[0]] = 1;
					// (*node_browsed)[nodes1[1]] = 1;
					// (*node_browsed)[nodes2[0]] = 1;
					// (*node_browsed)[nodes2[1]] = 1;

				}
				precedent_label = graphCut.whatLabel(i);
			}
		} else {nb_unary_edge++;}

    }
    if (speaking)
        std::cout << "We have "<< nb_turning_point << " turning points and "
                  << nb_unary_edge/2 << " unary egdes."<< '\n';

}
/*----------------------------------------------------------------------------*/
void PolycubeToolbox::gregson2011VariablesGeneration(){

	std::cout << "==updating variables :==" << '\n';
	std::cout << "computing faces' normal ..." << '\n';
	// computing faces variable, straight up from the mesh class
	for(auto f_id : mesh->faces()){
		Face f = mesh->get<Face>(f_id);
		(*face_normals)         [f_id] = f.normal();
		(*face_chosenNormal_ID) [f_id] = closestNormal((*face_normals)[f_id]);
		(*face_chosenNormal)    [f_id] = vect[(*face_chosenNormal_ID)[f_id]];
	}

	//smoothing up the result
	this->smoothing();
    this->computeNodeInformation();

	// marking particularly not on boundary faces
    for(auto f_id : mesh->faces()){
        Face f = mesh->get<Face>(f_id);
		if (!mesh->isMarked(f, mark_bnd_faces)){
			(*face_chosenNormal_ID)[f.id()]=-1;
		}
	}
	std::cout << "==Finished updating==" << '\n';

}
/*----------------------------------------------------------------------------*/
inline double PolycubeToolbox::angleFace2Node(TCellID face, TCellID node){
	Node n1 = mesh->get<Node>(node);
	std::vector<TCellID> nodes = (mesh->get<Face>(face)).getIDs<Node>();
	std::vector<math::Vector> v;
	for (auto n : nodes){
		if (n != node){
			Node n2 = mesh->get<Node>(n);
			v.push_back(math::Vector({n2.X() - n1.X(), n2.Y() - n1.Y(), n2.Z() - n1.Z()}));
		}
	}
	return v[0].angle(v[1]);
}


/*----------------------------------------------------------------------------*/
void PolycubeToolbox::computeNodeInformation(){
	std::cout << "computing node information ..." << '\n';
	// compute normal and quaternion of regular boundary point
	for (auto n_id : mesh->nodes()) {
		math::Vector normal({0, 0, 0});
		math::Vector ChosenNormal({0, 0, 0});
        Node n = mesh->get<Node>(n_id);
        if(mesh->isMarked(n,mark_smooth_bnd_nodes)){
            std::vector<TCellID> faces =n.getIDs<Face>();
            for (auto i : faces){
                if (mesh->isMarked<Face>(i,mark_bnd_faces)){
                    double angle = angleFace2Node(i, n.id());
                    normal=normal + angle * (*face_normals)[i];
                    ChosenNormal = ChosenNormal + angle * (*face_chosenNormal)[i];
                }
            }
            normal.normalize();
            ChosenNormal.normalize();
            (*point_normals)[n.id()] = normal;
            (*point_chosenNormal)[n.id()] = ChosenNormal;
            (*point_quaternion)[n.id()].setFromTwoVectors(math2eigen((*point_normals)[n.id()]),
                                                          math2eigen((*point_chosenNormal)[n.id()]));

        }
    }

}
/*----------------------------------------------------------------------------*/
void PolycubeToolbox::gregson2011Iteration(){

	if (GREGSON2011_IS_INIT){
		std::cout << "------------" << "Iteration " << ++GREGSON2011_ITERATION_NUMBER<< "------------" << '\n';
		if (GREGSON2011_RECOMPUTING_VARIABLES_IN_ITERATIONS){
            this->gregson2011VariablesGeneration();
		}
        this->computeNodeInformation();

		Matrix<double,Dynamic,4> smooth_quaternion(mesh->getNbNodes(),4);
		Matrix<double,Dynamic,4> new_quaternion(mesh->getNbNodes(),4);

		Matrix<double,Dynamic,3> scnd_member(mesh->getNbNodes(),3);
		Matrix<double,Dynamic,3> node_new_centers(mesh->getNbNodes(),3);

		smooth_quaternion << Matrix<double,Dynamic,4>::Zero(mesh->getNbNodes(),4); // initialisation to 0
		// making a eigen matrix to resolve from quaternion
		for (auto n_id : mesh->nodes()) {
			math::Vector normal({0, 0, 0});
			Node n = mesh->get<Node>(n_id);
			if(mesh->isMarked(n,mark_smooth_bnd_nodes)){
				smooth_quaternion(n.id(),0) = (*point_quaternion)[n_id].w();
				smooth_quaternion(n.id(),1) = (*point_quaternion)[n_id].x();
				smooth_quaternion(n.id(),2) = (*point_quaternion)[n_id].y();
				smooth_quaternion(n.id(),3) = (*point_quaternion)[n_id].z();
			} // end if
		} // end for

		std::cout << "Solving laplace parse system ...";
		std::cout.flush();
		const clock_t adjacency_system_time = clock();
		new_quaternion = adjacency_solver.solve(smooth_quaternion);
		std::cout << " ok. Computing time : " << float( clock () - adjacency_system_time ) /  CLOCKS_PER_SEC << "sec" << endl;


		//convert back to quaternion
		for (auto n_id : mesh->nodes()) {
			Node n = mesh->get<Node>(n_id);
            (*point_quaternion)[n.id()] = Eigen::Quaternion<double>(new_quaternion(n_id,0),
                                                                    new_quaternion(n_id,1),
                                                                    new_quaternion(n_id,2),
                                                                    new_quaternion(n_id,3));
            (*point_quaternion)[n.id()].normalize();
			(*point_rotationMatrix)[n.id()] = (*point_quaternion)[n.id()].toRotationMatrix();
		} // end for


		//compute the second member of the poisson equation
		scnd_member <<   Matrix<double,Dynamic,3>::Zero(mesh->getNbNodes(),3);
		for (auto n_id : mesh->nodes()) {
			Node n = mesh->get<Node>(n_id);
			std::vector<TCellID> edges = n.getIDs<Edge>();
			for (auto e:edges){
				std::vector<TCellID> nodes = (mesh->get<Edge>(e)).getIDs<Node>();
				for (auto e_n : nodes) {
					if (!(e_n == n.id())){
						Eigen::Matrix3d m = ((*point_rotationMatrix)[n.id()]+(*point_rotationMatrix)[e_n])/2;
						scnd_member.row(n.id()) += m *Vector3d(n.X() - (mesh->get<Node>(e_n)).X(),
						n.Y() - (mesh->get<Node>(e_n)).Y(),
						n.Z() - (mesh->get<Node>(e_n)).Z());
					} // end if
				} // end for (auto e_n : nodes)
			} // end for (auto e:edges)
		} // end for (IGMesh::node_iterator itn = mesh->nodes_begin(); !itn.isDone(); itn.next())
		scnd_member.row(0) = Vector3d(0, 0,0);

		std::cout << "Solving Poisson parse system ...";
		std::cout.flush();
		const clock_t poisson_system_time = clock();
		node_new_centers = poisson_solver.solve(scnd_member);
		std::cout << " ok. Computing time : " << float( clock () - poisson_system_time ) /  CLOCKS_PER_SEC << "sec" << endl;


		for (auto n_id : mesh->nodes()) {
			Node n = mesh->get<Node>(n_id);
			n.setX(node_new_centers(n_id,0));
			n.setY(node_new_centers(n_id,1));
			n.setZ(node_new_centers(n_id,2));
		}

		if (GREGSON2011_WRITING_RESULT){
			std::cout << "Writing results ..." << '\n';

            gmds::IGMeshIOService ioService(mesh);
            gmds::VTKWriter w(&ioService);
            w.setCellOptions(gmds::N|gmds::F);
            w.setDataOptions(gmds::N|gmds::F);
            w.write(string("result") + std::to_string(GREGSON2011_ITERATION_NUMBER));
        }
	} // IF IS INIT
	else
	{
		std::cout << "YOU TRIED RUNNING GREGSON2011 ITERATION WITHOUT INITIALISATION" << '\n';
	}
}
/*----------------------------------------------------------------------------*/
void PolycubeToolbox::gregson2011SetWritingResult(bool b){
	std::cout << "GREGSON2011 writing result set on : " << b << '\n';
	this->GREGSON2011_WRITING_RESULT = b;
}
/*----------------------------------------------------------------------------*/
void PolycubeToolbox::gregson2011SetVariableComputationInIteration(bool b){
	std::cout << "GREGSON2011 variable computation in iteration set on : " << b << '\n';
	if (!b){
		std::cout << "Be careful, no further check will be done, and without required "
		<< "computation, errors will occur." << '\n';
	}
	this->GREGSON2011_RECOMPUTING_VARIABLES_IN_ITERATIONS = b;
}
/*----------------------------------------------------------------------------*/
void PolycubeToolbox::gregson2011Init(){
	std::cout << "====INITIALISATION OF GREGSON2011====" << '\n';
	// computing adjacency and poisson matrix
	std::cout << "Creating point's adjacency and poisson's matrix ...";
	std::cout.flush();
	const clock_t matrice_computation_time = clock();
	adjacency_matrix = computeAdjacency();

	poisson_matrix = computePoisson();

	std::cout << " ok. Computing time : " << float( clock () - matrice_computation_time ) /  CLOCKS_PER_SEC << "sec" << endl;
	std::cout << "We are working on " << mesh->getNbNodes() << "x"
						<< mesh->getNbNodes() << " matrices with " << adjacency_matrix.nonZeros()
						<<" and "<< poisson_matrix.nonZeros() << " non zeros values."<< '\n';


	std::cout << "Computing LU decompositions : " << '\n';

	std::cout << "Poisson ...";
	std::cout.flush();
	const clock_t poisson_LU_time = clock();
	poisson_solver.compute(poisson_matrix);
	std::cout << " ok. Computing time : " << float( clock () - poisson_LU_time ) /  CLOCKS_PER_SEC << "sec" << endl;

	std::cout << "Laplace ...";
	std::cout.flush();
	const clock_t laplace_LU_time = clock();
	adjacency_solver.compute(adjacency_matrix);
	std::cout << " ok. Computing time : " << float( clock () - laplace_LU_time ) /  CLOCKS_PER_SEC << "sec" << endl;

	GREGSON2011_ITERATION_NUMBER = 0;
	GREGSON2011_WRITING_RESULT = true;
	GREGSON2011_IS_INIT = true;
	GREGSON2011_RECOMPUTING_VARIABLES_IN_ITERATIONS = true;
	std::cout << "Writing result : True" << '\n';
	std::cout << "Recomputing variables in iterations : True" << '\n';

	std::cout << "====INITIALISATION FINISHED====" << '\n';
}
/*----------------------------------------------------------------------------*/
void  PolycubeToolbox::gregson2011Run(int nb_iteration){
	std::cout << "==========RUNNING GREGSON2011 ALGORITHM ==========" << '\n';
    this->gregson2011Init();
	std::cout << "====     ITERATING    ====" << '\n';
	std::cout << "We will now have " << nb_iteration << " iterations to convert "
						<<"the mesh to a closer approximation of a Polycube :" << '\n';
	for (int iteration = 0; iteration < nb_iteration; iteration++)
	{
        this->gregson2011Iteration();

	} // end (int iteration = 0; iteration < n; iteration++)
	std::cout << "====FINISHED ITERATING====" << '\n';
    this->gregson2011VariablesGeneration();
	std::cout << "====PROCESSING THE POLYCUBE====" << '\n';
    this->projectionPreprocess();
    this->projectOnPolycube();
	std::cout << "====END POLYCUBE PROCESSING====" << '\n';

	std::cout << "==========END GREGSON2011 ALGORITHM ==========" << '\n';


}
/*----------------------------------------------------------------------------*/
typedef std::vector<TCellID> chart;
/*----------------------------------------------------------------------------*/
void PolycubeToolbox::projectOnPolycube(){
	std::cout << "Now starting projection onto a polycube." << '\n';

	if (charts.empty()){
        this->chartsGeneration();
	}

	// we now have charts of faces.
	// Each chart will be project the farthest into a certain direction
	// we will therefore found a maximum in this direction for each chart
	// then project every node of the chart to this maximum

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList_X,tripletList_Y,tripletList_Z;


	SparseVector<double> node_scnd_member_X(mesh->getNbNodes(),1);
	SparseVector<double>node_scnd_member_Y(mesh->getNbNodes(),1);
	SparseVector<double> node_scnd_member_Z(mesh->getNbNodes(),1);

	Matrix<double,Dynamic,1> new_center_X(mesh->getNbNodes(),1);
	Matrix<double,Dynamic,1> new_center_Y(mesh->getNbNodes(),1);
	Matrix<double,Dynamic,1> new_center_Z(mesh->getNbNodes(),1);


	// node_scnd_member_X = Matrix<double,Dynamic,1>::Zero(mesh->getNbNodes(),1);
	// node_scnd_member_Y = Matrix<double,Dynamic,1>::Zero(mesh->getNbNodes(),1);
	// node_scnd_member_Z = Matrix<double,Dynamic,1>::Zero(mesh->getNbNodes(),1);

	Vector3d working_vector;
	std::vector<TCellID> chart_nodes;
	double max_value, set_value;

	std::cout << "Computing faces of polycube ..." << '\n';
	// Now we want to push each chart to a border of a polycube
	int chart_id=0;
	for (auto c : charts)
	{
		// we take an arbitrary small value, to make sure we don't encounter smaller value
		max_value =-10E80;
		set_value =0;

		// first we find the maximum
		for (auto f :c)
		{
			std::vector<TCellID> face_nodes = (mesh->get<Face>(f)).getIDs<Node>();
			for (auto node:face_nodes)
			{
				// to avoid long syntax
				Node n = mesh->get<Node>(node);
				// working_vector represent in fact the point vector
				working_vector = Vector3d(n.X(),n.Y(),n.Z());
				/* the value of (working_vector.transpose()*	math2eigen((*edge_chosenNormal)[e])
				 will be only depending the axis of the chosen normal, and it will be positive
				 when maximum in the direction of the normal.*/
				if (working_vector.transpose()*	math2eigen((*face_chosenNormal)[f]) > max_value){
					// this value is the real value of the maximum, not only positive
					set_value = Vector3d(1.,1.,1.).transpose() * (math2eigen((*face_chosenNormal)[f])
												* (working_vector.transpose()*	math2eigen((*face_chosenNormal)[f]) ));
					max_value = working_vector.transpose()*	math2eigen((*face_chosenNormal)[f]);
				}
		}
		} // end for (auto f :c)
		// set_value /= c.size()*3;

		// we now have a maximum,
		// we set all chart's points to this maximum
		for (auto f :c)
		{
			std::vector<TCellID> face_nodes = (mesh->get<Face>(f)).getIDs<Node>();

			if((*face_chosenNormal)[f].X() != 0)
			// the normal on x is non null, we are on the direction of x
			{
				// we change the value of all the node of the face
				for (auto node:face_nodes){
					// (mesh->get<Node>(node)).setX(set_value);
					if (node_scnd_member_X.coeffRef(node)==0){
						if(chart_id==1 || chart_id==21 || chart_id==6)
							node_scnd_member_X.coeffRef(node) = set_value-5;
                        else if(chart_id==17)
                            node_scnd_member_X.coeffRef(node) = set_value+10;
                        else
							node_scnd_member_X.coeffRef(node) = set_value;
						(*node_updatedX)[node] = 1;
						tripletList_X.push_back(T(node,node,1));
					}
				}
			}
                // else
            else if ((*face_chosenNormal)[f].Y() != 0)
                // the normal on y is non null, we are on the direction of y
            {
                for (auto node:face_nodes){
                    // (mesh->get<Node>(node)).setY(set_value);
                    if (node_scnd_member_Y.coeffRef(node) ==0){
                        if(chart_id==24 || chart_id==16)
                            node_scnd_member_Y.coeffRef(node) = set_value-5;
                        else if(chart_id==15 )
                            node_scnd_member_Y.coeffRef(node) = set_value+10;
                        else
                            node_scnd_member_Y.coeffRef(node) = set_value;
                        (*node_updatedY)[node] = 1;
                        tripletList_Y.push_back(T(node,node,1));
                    }
                }
            }
            else
                // else we are on the direction of z
            {
                for (auto node:face_nodes){
                    // (mesh->get<Node>(node)).setZ(set_value);
                    if (node_scnd_member_Z.coeffRef(node)==0){
                        node_scnd_member_Z.coeffRef(node) = set_value;
                        (*node_updatedZ)[node] = 1;
                        tripletList_Z.push_back(T(node,node,1));
                    }
                }
            }

        } // end for (auto f :c)
        chart_id++;
    } // end for (auto c:chart)



    /* this last part is not optimised, as it can change the value of a same point
      Multiple time, and we could maybe move the for (I don't know about compiler optimisation)
      */



	for (auto n_id : mesh->nodes()) {

		Node n = mesh->get<Node>(n_id);
		if(node_scnd_member_X.coeffRef(n.id()) ==0 ){
			tripletList_X.push_back(T(n.id(),n.id(),n.nbEdges()));
		}
		if(node_scnd_member_Y.coeffRef(n.id())  ==0 ){
			tripletList_Y.push_back(T(n.id(),n.id(),n.nbEdges()));
		}
		if(node_scnd_member_Z.coeffRef(n.id())  ==0 ){
			tripletList_Z.push_back(T(n.id(),n.id(),n.nbEdges()));
		}

		std::vector<TCellID> edges = n.getIDs<Edge>();
		for (auto e :edges){
			TCellID other_node = ( (mesh->get<Edge>(e)).getIDs<Node>()[0]!=n.id() ? (mesh->get<Edge>(e)).getIDs<Node>()[0] : (mesh->get<Edge>(e)).getIDs<Node>()[1] );
			if(node_scnd_member_X.coeffRef(n.id())  ==0 ){
				tripletList_X.push_back(T(n.id(),other_node,-1));
			}
			if(node_scnd_member_Y.coeffRef(n.id())  ==0 ){
				tripletList_Y.push_back(T(n.id(),other_node,-1));
			}
			if(node_scnd_member_Z.coeffRef(n.id())  ==0 ){
				tripletList_Z.push_back(T(n.id(),other_node,-1));
			}
		}

	}

	SpMat laplace_matrix_X(mesh->getNbNodes(),mesh->getNbNodes());
	SpMat laplace_matrix_Y(mesh->getNbNodes(),mesh->getNbNodes());
	SpMat laplace_matrix_Z(mesh->getNbNodes(),mesh->getNbNodes());

	laplace_matrix_X.setFromTriplets(tripletList_X.begin(), tripletList_X.end());
	laplace_matrix_Y.setFromTriplets(tripletList_Y.begin(), tripletList_Y.end());
	laplace_matrix_Z.setFromTriplets(tripletList_Z.begin(), tripletList_Z.end());

	laplace_matrix_X.makeCompressed();
	laplace_matrix_Y.makeCompressed();
	laplace_matrix_Z.makeCompressed();


	std::cout << "Moving other point to have homogeneous repartition" << '\n';

	SparseLU<SparseMatrix<double>> solver_X,solver_Y,solver_Z;


	std::cout << "Computing LU decompositions ... ";
	std::cout.flush();

	const clock_t poisson_LU_time = clock();
	// solver_X.compute(laplace_matrix_X);
	// solver_Y.compute(laplace_matrix_Y);
	// solver_Z.compute(laplace_matrix_Z);

	std::cout << " ok. Computing time : " << float( clock () - poisson_LU_time ) /  CLOCKS_PER_SEC << "sec" << endl;


	std::cout << "Testing GMRES ... ";
	std::cout.flush();

	const clock_t poisson_GMRES = clock();

	SparseVector<double> X = node_scnd_member_X;
	SparseVector<double> Y = node_scnd_member_Y;
	SparseVector<double> Z = node_scnd_member_Z;

	gmresAya(laplace_matrix_X,X,node_scnd_member_X);
	gmresAya(laplace_matrix_Y,Y,node_scnd_member_Y);
	gmresAya(laplace_matrix_Z,Z,node_scnd_member_Z);

	std::cout << " ok. Computing time : " << float( clock () - poisson_GMRES ) /  CLOCKS_PER_SEC << "sec" << endl;

	std::cout << "Solving laplace parse system ..";
	std::cout.flush();
	const clock_t solver_computation_time = clock();

	// new_center_X = solver_X.solve(node_scnd_member_X);
	// new_center_Y = solver_Y.solve(node_scnd_member_Y);
	// new_center_Z = solver_Z.solve(node_scnd_member_Z);

 // solver_X.solve(node_scnd_member_X);
	//  solver_Y.solve(node_scnd_member_Y);
	//  solver_Z.solve(node_scnd_member_Z);

	std::cout << " ok. Computing time : " << float( clock () - solver_computation_time ) /  CLOCKS_PER_SEC << "sec" << endl;


	for (auto n_id : mesh->nodes()) {
		Node n = mesh->get<Node>(n_id);
		n.setX(X.coeffRef(n_id));
		n.setY(Y.coeffRef(n_id));
		n.setZ(Z.coeffRef(n_id));
	}

}

/*----------------------------------------------------------------------------*/
void PolycubeToolbox::projectionPreprocess(){
    this->gregson2011VariablesGeneration();

	std::cout << "Preparing mesh for position-driven deformation to Polycube" << '\n';
	std::cout << "Splitting multi-orientation charts ... (not implemented)" << '\n';
	std::cout << "Splitting highly non-plannar charts ... (not implemented)" << '\n';
}
/*----------------------------------------------------------------------------*/
void PolycubeToolbox::chartsGeneration(bool speaking){
	int current_chart = 0;
    chart face_to_explore;
    int marked_face = mesh->newMark<Face>();
    charts.clear();

    if (speaking)
        std::cout << "||Creating or updating charts ..." ;
    //we divide the mesh's faces in charts, depending on their chosen normal
    for(auto f_id : mesh->faces())
    {
		Face chosen_face = mesh->get<Face>(f_id);
		if (!(mesh->isMarked(chosen_face,marked_face)) && (mesh->isMarked(chosen_face,mark_bnd_faces)))
		{
			// creating the space for the new chart
			charts.push_back(chart());
			// adding the first face into the face that we will explore
			face_to_explore.push_back(chosen_face.id());

			while (!face_to_explore.empty())
			{
				// adding face to the chart
				chosen_face = mesh->get<Face>(face_to_explore.back());
				face_to_explore.pop_back();
				charts[current_chart].push_back(chosen_face.id());
				mesh->mark(chosen_face,marked_face);
				//looking for other similar faces
				std::vector<TCellID> edges = chosen_face.getIDs<Edge>();
				for (auto e:edges){// only 3 edges per face
					// if the edge is a "smooth" edge
					if (!mesh->isMarked(mesh->get<Edge>(e),mark_sharp_edges)){
							std::vector<TCellID> faces = (mesh->get<Edge>(e)).getIDs<Face>();
							for (auto f : faces){
								if ( (*face_chosenNormal)[chosen_face.id()] == (*face_chosenNormal)[f]){
									if (!mesh->isMarked(mesh->get<Face>(f),marked_face)&& mesh->isMarked(mesh->get<Face>(f),mark_bnd_faces)){
										face_to_explore.push_back(f);
										mesh->mark(mesh->get<Face>(f),marked_face);

										// we already marked the chosen_face, so it won't execute back in face_to_explore
									}
								}
							} // end for (auto f : faces)
						} // if smooth edge
				} // end for (auto e:edges)

			}// end while
			current_chart++; //we move onto the next chart
		}
	} // end for all mesh's faces
	mesh->unmarkAll<Face>(marked_face);
	mesh->freeMark<Face>(marked_face);

	// display the chart of each face
	int num = 0;
	for (auto c: charts){
		for (auto f :c){
			(*face_chart)[f] = num;
		}
		num++;
	}
	if (speaking)
	std::cout << "We have " << num << " charts" << '\n';

	if (speaking)
	std::cout << "||Detecting chart' edges and corners ..." ;
	//cleaning variable
	for (auto n_id : mesh->nodes()) {
		(*node_chartCorner)[n_id] = 0;
	}

	std::cout.flush();
	charts_edges.clear();
	marked_face = mesh->newMark<Face>();
	int marked_edge = mesh->newMark<Edge>();
	int chart_number =0;

	int chart_corner_number = 0;

	for (auto working_chart: charts)
	{
		std::vector<TCellID> edges_on_chart_bnd;

		for (auto face : working_chart){
			mesh->mark<Face>(face,marked_face);
		}
		for (auto face : working_chart){
			std::vector<TCellID> edges = (mesh->get<Face>(face)).getIDs<Edge>();
			for (auto e: edges) {
				std::vector<TCellID> bnd_faces, faces = (mesh->get<Edge>(e)).getIDs<Face>();
				for (auto f : faces){
					if (mesh->isMarked<Face>(f, mark_bnd_faces)) {
						bnd_faces.push_back(f);
					}
				}
				if (!mesh->isMarked<Face>(bnd_faces[1], marked_face) || !mesh->isMarked<Face>(bnd_faces[0], marked_face)) {
					edges_on_chart_bnd.push_back(e);
					mesh->mark<Edge>(e,marked_edge);

				}
			}
		}

		std::vector<TCellID> corner_of_this_chart;
		for (auto e: edges_on_chart_bnd){
			if (mesh->isMarked<Edge>(e,marked_edge)){
				mesh->unmark<Edge>(e,marked_edge);
				int other_chart;

				std::vector<TCellID> wanted_faces, edges_faces = (mesh->get<Edge>(e)).getIDs<Face>();
				for (auto f : edges_faces){
					if (mesh->isMarked<Face>(f,mark_bnd_faces)){
						wanted_faces.push_back(f);
					}
				}
				if ((*face_chart)[wanted_faces[1]] == chart_number){
					other_chart=(*face_chart)[wanted_faces[0]];
				} else {
					other_chart=(*face_chart)[wanted_faces[1]];
				}


				Triplet_charts_edge ce;
				ce.chart1 = chart_number;
				ce.chart2 = other_chart;
				ce.edge.push_back(e);

				std::vector<TCellID> corners_of_the_edge;
				std::vector<TCellID> node_to_explore = (mesh->get<Edge>(e)).getIDs<Node>();

				while (!node_to_explore.empty())
				{
					Node node = mesh->get<Node>(node_to_explore.back());
					node_to_explore.pop_back();
					bool is_corner = true;
					TCellID chosen_edge;
					std::vector<TCellID> faces = node.getIDs<Face>();
					bool more_than_2_chart = false;
					bool this_node_is_already_a_corner = false;
					for (auto n : corner_of_this_chart){
						if (n == node.id()){
							this_node_is_already_a_corner = true;
						}
					}
					for (auto f : faces) {
						if (mesh->isMarked<Face>(f,mark_bnd_faces) && (*face_chart)[f] != other_chart && (*face_chart)[f] != chart_number && !this_node_is_already_a_corner){
							more_than_2_chart = true;
						}
					}
					if (! more_than_2_chart){
						std::vector<TCellID> edges = node.getIDs<Edge>();
						for (auto nodes_edge : edges){
							if (mesh->isMarked<Edge>(nodes_edge,marked_edge))
							{
								wanted_faces.clear();
								edges_faces = (mesh->get<Edge>(nodes_edge)).getIDs<Face>();
								for (auto f : edges_faces){
									if (mesh->isMarked<Face>(f,mark_bnd_faces)){
										wanted_faces.push_back(f);
										chosen_edge = nodes_edge; is_corner =false;
									}
								}

								if ((*face_chart)[wanted_faces[1]] == chart_number){
									other_chart=(*face_chart)[wanted_faces[0]];
								} else {
									other_chart=(*face_chart)[wanted_faces[1]];
								}
							}
						}
					}
					if (is_corner){
						if ((*node_chartCorner)[node.id()] !=1) {
							(*node_chartCorner)[node.id()]=1; chart_corner_number++;
							corner_of_this_chart.push_back(node.id());
						}
						corners_of_the_edge.push_back(node.id());
					} else if (other_chart == ce.chart2){
						ce.edge.push_back(chosen_edge);
						mesh->unmark<Edge>(chosen_edge,marked_edge);
						std::vector<TCellID> nodes = (mesh->get<Edge>(chosen_edge)).getIDs<Node>();
						if (nodes[0] == node.id()) {
							node_to_explore.push_back(nodes[1]);
						}else {
							node_to_explore.push_back(nodes[0]);
						}
					}
					else {
						if ((*node_chartCorner)[node.id()] !=1) {
							(*node_chartCorner)[node.id()]=1; chart_corner_number++;
						}
						corners_of_the_edge.push_back(node.id());
					}

				}
				ce.corner1 = corners_of_the_edge[0];
				ce.corner2 = corners_of_the_edge[1];
				charts_edges.push_back(ce);
			}


		}


		for (auto face : working_chart){
			mesh->unmark<Face>(face,marked_face);
		}

		chart_number++;

	}




	// being here, all marking should be already erased, we don't need to unmarkAll (time consuming)
	mesh->freeMark<Face>(marked_face);
	mesh->freeMark<Edge>(marked_edge);

	// erasing duplicated edges
	std::vector<Triplet_charts_edge> new_charts_edges;
	std::vector<bool> already_taken(charts_edges.size(),false);
	for (int i=0;i<charts_edges.size();i++){
		if (!already_taken[i]){
			already_taken[i] = true;
			for (int j=0;j<charts_edges.size();j++){
				bool out = false;
				if (!already_taken[j] ){
					for (auto e : charts_edges[j].edge){
						if (e == charts_edges[i].edge[0]){
							new_charts_edges.push_back(charts_edges[i]);
							already_taken[j] = true;
							out = true;
							break;
						}
					}
				}
				if (out) { break;}
			}
		}
	}



	charts_edges = new_charts_edges;

	if (speaking)
	std::cout << "We have " << charts_edges.size() << " edges of chart "
	<< "and " << chart_corner_number << " corners." << '\n';

	if (speaking)
	std::cout << "|| ordering edges of charts..." << '\n';
	// we arrange the edge to simplify further computation

	gmds::Variable<math::Vector>* node_edge_normals;
	try {
		node_edge_normals = mesh->newVariable<math::Vector,GMDS_NODE>("edge_normal");
	} catch (const GMDSException& e) {node_edge_normals= mesh->getVariable<math::Vector,GMDS_NODE>("edge_normal");}

	for (auto& TEDGE : charts_edges){
		chart_edge new_edge;
		std::vector<math::Vector> ns_node;
		std::vector<TCellID> nodeID;
		TCellID node = TEDGE.corner1;
		ns_node.push_back(math::vec(mesh->get<Node>(node).point()));
		long long last_edge=-1;
		while (new_edge.size() !=  TEDGE.edge.size())
		{
			for (auto e : TEDGE.edge){

				if (last_edge!= e){
					std::vector<TCellID> nodes = (mesh->get<Edge>(e)).getIDs<Node>();
					if (nodes[1]==node){
						Node n1 = mesh->get<Node>(node);
						Node n2 = mesh->get<Node>(nodes[0]);
						// (*edge_normal)[e] = math::Vector(n1.X() - n2.X(),n1.Y() - n2.Y(),n1.Z() - n2.Z());

						ns_node.push_back(math::Vector({n2.X(), n2.Y(), n2.Z()})); nodeID.push_back(nodes[0]);
						node = nodes[0];
						new_edge.push_back(e); last_edge = e;
					} else if (nodes[0] == node){
						Node n1 = mesh->get<Node>(node);
						Node n2 = mesh->get<Node>(nodes[1]);
						ns_node.push_back(math::Vector({n2.X(), n2.Y(), n2.Z()})); nodeID.push_back(nodes[1]);
						node = nodes[1];
						new_edge.push_back(e); last_edge =e;
					}
				}

			}
		}

		std::vector<math::Vector> smoothed_node;
		smoothed_node.push_back(ns_node.front());
		for (size_t i = 1; i < ns_node.size()-1; i++) {
			smoothed_node.push_back((ns_node[i-1] + 2* ns_node[i] + ns_node[i+1])/4);
		}
		smoothed_node.push_back(ns_node.back());

		for (size_t i = 0; i < new_edge.size(); i++) {
			(*edge_normal)[new_edge[i]] = smoothed_node[i+1] - smoothed_node[i];
			(*edge_normal)[new_edge[i]].normalize();
			(*node_edge_normals)[nodeID[i]] = (*edge_normal)[new_edge[i]];
		}

		 TEDGE.edge = new_edge;


	}

}
/*----------------------------------------------------------------------------*/
void PolycubeToolbox::smoothing(){
	bool normal_was_changed=true;
	int iteration =0;
	int nbChange=0;
	std::cout << "-Beginning smoothing :" << '\n';
	std::cout << "|Face smoothing ..." << '\n';
	while (normal_was_changed){
		normal_was_changed = false;
		for(auto f_id : mesh->faces()){
			Face f = mesh->get<Face>(f_id);
			if (mesh->isMarked(f,mark_bnd_faces)){
				std::vector<TCellID> linked_face;
				std::vector<TCellID> via_smooth_edge;
				std::vector<TCellID> edges = f.getIDs<Edge>();
				for (auto e: edges){
					std::vector<TCellID> faces=(mesh->get<Edge>(e)).getIDs<Face>();
					for (auto face: faces) {
						if (face != f.id()&& mesh->isMarked(mesh->get<Face>(face),mark_bnd_faces)){
							linked_face.push_back(face);

							via_smooth_edge.push_back(!mesh->isMarked(mesh->get<Edge>(e),mark_sharp_edges));
						}
					}
				}
				if ((*face_chosenNormal_ID)[linked_face[0]] == (*face_chosenNormal_ID)[linked_face[1]]
				&& (*face_chosenNormal_ID)[f.id()] != (*face_chosenNormal_ID)[linked_face[1]]
				&& via_smooth_edge[0] && via_smooth_edge[1]){
					(*face_chosenNormal_ID)[f.id()] = (*face_chosenNormal_ID)[linked_face[1]];
					normal_was_changed = true; nbChange++;
				} else
				if ((*face_chosenNormal_ID)[linked_face[1]] == (*face_chosenNormal_ID)[linked_face[2]]
				&& (*face_chosenNormal_ID)[f.id()] != (*face_chosenNormal_ID)[linked_face[1]]
				&& via_smooth_edge[1] && via_smooth_edge[2]){
					(*face_chosenNormal_ID)[f.id()] = (*face_chosenNormal_ID)[linked_face[1]];
					normal_was_changed = true; nbChange++;
				} else
				if ((*face_chosenNormal_ID)[linked_face[0]] == (*face_chosenNormal_ID)[linked_face[2]]
				&& (*face_chosenNormal_ID)[f.id()] != (*face_chosenNormal_ID)[linked_face[0]]
				&& via_smooth_edge[0] && via_smooth_edge[2]){
					(*face_chosenNormal_ID)[f.id()] = (*face_chosenNormal_ID)[linked_face[0]];
					normal_was_changed = true; nbChange++;
				}
			}
		}
		iteration++;
	}
	std::cout << "|nb of iteration : "<<iteration << " for "<< nbChange << " changes" <<'\n';

	std::cout << "|[chart generation]" << '\n';
    this->chartsGeneration();
	std::cout << "|[end chart generation]" << '\n';


	std::cout << "|chart smoothing ..." << '\n';
	for (auto c : charts) {
		if (c.size() == 1){
			TCellID face = c[0];
			std::vector<TCellID> edges=(mesh->get<Face>(face)).getIDs<Edge>();
			for (auto e : edges) {
				std::vector<TCellID> faces=(mesh->get<Edge>(e)).getIDs<Face>();
				if (!mesh->isMarked(mesh->get<Edge>(e),mark_sharp_edges)){
					for (auto f : faces){
						if (f != face && mesh->isMarked(mesh->get<Face>(f),mark_bnd_faces)){

							(*face_chosenNormal_ID)[face] =(*face_chosenNormal_ID)[f];
							(*face_chosenNormal)[face] = (*face_chosenNormal)[f];
						}
					}
				} // if smooth edge
			} // for (auto e : egdes)
		} // if (c.size() == 1)
	}// for (auto c : charts)
	std::cout << "|[chart generation]" << '\n';
    this->chartsGeneration();
	std::cout << "|[end chart generation]" << '\n';
	std::cout << "-Finished smoothing." << '\n';
}


/*----------------------------------------------------------------------------*/
labelling PolycubeToolbox::closestNormal(math::Vector v){ //verifier optimalité
	double angle= v.dot(vect[0]);

	labelling label = Xp;

	for (int i = 1; i < 6; i++){
		if (v.dot(vect[i])>angle+ANGLE_DETERMINENT){
			angle =v.dot(vect[i]);
            label = (labelling)i;
		}
	}
	return label;
}


/*----------------------------------------------------------------------------*/
inline int oppositeNormale(int ANormal){
	switch (ANormal) {
		case 0: return 1;
		case 1: return 0;
		case 2: return 3;
		case 3: return 2;
		case 4: return 5;
		case 5: return 4;
	}
	return 8;
}

/*----------------------------------------------------------------------------*/
void PolycubeToolbox::
graphCutHillClimbing(double ADeltaStart, double ADeltaIncr, double ABoundingBoxRatio)
{

    this->chartsGeneration();
	this->compute_turning_point();


	int marked_face = mesh->newMark<Face>();
	bool There_is_still_turning_points = false;
	int writing_advencement= 0;
	int nb_turning_point = 0;

	std::vector<std::vector<Triplet_charts_edge>> edges_of_each_charts;
	for (auto c : charts){
		edges_of_each_charts.push_back(std::vector<Triplet_charts_edge>());
	}
	std::vector<bool> non_monotone(charts.size(),false);
	// assigning charts_edges to charts
	for (auto TEDGE : charts_edges){
		edges_of_each_charts[TEDGE.chart1].push_back(TEDGE);
		edges_of_each_charts[TEDGE.chart2].push_back(TEDGE);
		nb_turning_point += TEDGE.turningPoints.size();
		if (!TEDGE.turningPoints.empty()){
			non_monotone[TEDGE.chart1] = true;
			non_monotone[TEDGE.chart2] = true;
			There_is_still_turning_points = true;
		}
	}

	std::cout << "====" <<nb_turning_point << " turning points to take into account ====" << '\n';

	while (There_is_still_turning_points) {

		for (size_t i = 0; i < mesh->getNbFaces(); i++) {
			(*face_constrained)[i] = 0;
		}

		for (size_t c = 0; c < charts.size(); c++) {
			if (non_monotone[c]){
				for (auto f : charts[c]){
					mesh->mark<Face>(f,marked_face);
				}
			}
		}

		std::vector<TCellID> non_monotone_faces;
		map<TCellID,int> linker_ID2faces; // useful to know to which id was assigned a face.
		mesh->deleteVariable(GMDS_FACE,"non_monotone_faces");
		gmds::Variable<int>* 	var_non_monotone_faces = mesh->newVariable<int,GMDS_FACE>("non_monotone_faces");
		// putting non monotone faces into a vector
		non_monotone_faces.clear();
		linker_ID2faces.clear();
		int num=0;
		for(auto f_id:mesh->faces()){
			Face f = mesh->get<Face>(f_id);
			if (mesh->isMarked(f,marked_face)) {
				non_monotone_faces.push_back(f.id());
				(*var_non_monotone_faces)[f.id()]=1;
				linker_ID2faces.insert( std::pair<TCellID,int>(f.id(),num) ); num++;
			}
		}


		int nb_faces = non_monotone_faces.size();

		GCoptimizationGeneralGraph graphCut(nb_faces,6);
		//setting the Neighbors
		for (size_t i = 0; i < nb_faces; i++) {
			std::vector<TCellID> edges = ( mesh->get<Face>(non_monotone_faces[i]) ).getIDs<Edge>();
			for (auto e : edges) {
				std::vector<TCellID> faces = ( mesh->get<Edge>(e) ).getIDs<Face>();
				for (auto f : faces) {
					if (mesh->isMarked<Face>(f,marked_face) && non_monotone_faces[i] != f) {
						graphCut.setNeighbors(i, linker_ID2faces[f]);
					}
				}
			}
		}
		// setting adjacency cost
		std::vector<math::Vector> extraData;
		for (size_t i = 0; i < nb_faces; i++) {
			extraData.push_back( (mesh->get<Face>(non_monotone_faces[i])).normal() );
		}

		graphCut.setSmoothCost(&PolycubeToolbox::graphCutCpq, (void *)&extraData);

		// original fidelity cost for faces
		std::vector<double> dataCost_static;
		for (auto f : non_monotone_faces){
			for (size_t i = 0; i < 6; i++) {
				double nts=  (( mesh->get<Face>(f) ).normal().dot( vect[i] ) - 1) /0.2;
				dataCost_static.push_back(1-exp(-0.5*nts*nts));
			}
		}

		mesh->deleteVariable(GMDS_FACE,"non_monotone_edge");
		gmds::Variable<int>* 	non_monotone_edge = mesh->newVariable<int,GMDS_FACE>("non_monotone_edge");
		mesh->deleteVariable(GMDS_FACE,"disp_face");
		gmds::Variable<int>* 	disp_face = mesh->newVariable<int,GMDS_FACE>("disp_face");
		mesh->deleteVariable(GMDS_FACE,"part_disp_face");
		gmds::Variable<int>* 	part_disp_face = mesh->newVariable<int,GMDS_FACE>("part_disp_face");

		// // constraining face between a sharp edge in the non monotone area
		// for(IGMesh::edge_iterator ite= mesh->edges_begin(); !ite.isDone(); ite.next()){
		// 	Edge e = ite.value();
		// 	if (mesh->isMarked(e,mark_sharp_edges)){
		// 		std::vector<TCellID> wanted_faces;
		// 		for (auto f : e.getIDs<Face>()){
		// 			if (mesh->isMarked<Face>(f,marked_face)){
		// 				wanted_faces.push_back(f);
		// 			}
		// 		}
		// 		if (wanted_faces.size()==2){
		// 			dataCost_static[linker_ID2faces[wanted_faces[0]]*6+(*face_chosenNormal_ID)[wanted_faces[1]]] = 30;
		// 			dataCost_static[linker_ID2faces[wanted_faces[1]]*6+(*face_chosenNormal_ID)[wanted_faces[0]]] = 30;
		// 			(*part_disp_face)[wanted_faces[0]] = 1;
		// 			(*part_disp_face)[wanted_faces[1]] = 1;
		// 		}
		// 	}
		// }


		//constraining faces on boundary
		for (size_t c = 0; c < charts.size(); c++) {
			if (non_monotone[c]){
				for (auto TEDGE : edges_of_each_charts[c]){
					if (TEDGE.chart1 == c && !non_monotone[TEDGE.chart2])
					{

						for (auto e : TEDGE.edge){
							std::vector<TCellID> faces = (mesh->get<Edge>(e)).getIDs<Face>();
							for (auto f : faces) {
								if (mesh->isMarked<Face>(f,marked_face)){
									dataCost_static[linker_ID2faces[f]*6 +
                                            oppositeNormale((*face_chosenNormal_ID)[charts[TEDGE.chart2][0]])] = 50;
									if (mesh->isMarked<Edge>(e, mark_sharp_edges)){
										dataCost_static[linker_ID2faces[f]*6+(*face_chosenNormal_ID)[charts[TEDGE.chart2][0]]] = 30;
										(*disp_face)[f] = 1;
									}
									(*non_monotone_edge)[f]=1;
									(*face_constrained)[f]=1;

								}
							}
						}

					} else if (TEDGE.chart2 == c && !non_monotone[TEDGE.chart1])
					{

						for (auto e : TEDGE.edge){
							std::vector<TCellID> faces = (mesh->get<Edge>(e)).getIDs<Face>();
							for (auto f : faces) {
								if (mesh->isMarked<Face>(f,marked_face)){
									dataCost_static[linker_ID2faces[f]*6 +
                                            oppositeNormale((*face_chosenNormal_ID)[charts[TEDGE.chart1][0]])] = 50;
									if (mesh->isMarked<Edge>(e, mark_sharp_edges)){
										dataCost_static[linker_ID2faces[f]*6+(*face_chosenNormal_ID)[charts[TEDGE.chart1][0]]] = 30;
									}
									(*non_monotone_edge)[f]=1;
									(*face_constrained)[f]=1;
								}
							}
						}

					}
				}
			}
		}


		std::vector<int> labelling_save;
		labelling_save.reserve(nb_faces);
		double delta = ADeltaStart;
		double incr = ADeltaIncr;
		int new_nb_turning_point = nb_turning_point;
		int best_so_far = nb_turning_point;
		std::cout << "STATUS : BSF = " << best_so_far << ", NTP = " << nb_turning_point << '\n';
		while (best_so_far >= nb_turning_point){
			if (delta > 100 * incr + ADeltaStart){ incr*=2; std::cout << "BOOSTING INCREMENT" << '\n';}
			delta += incr;
			std::cout << "detla = "<< delta  << " : " << '\n';
			for (size_t label = 0; label < 6; label=label+2) {
				std::vector<double> dataCost;
				for (auto cost : dataCost_static){
					dataCost.push_back(cost);
				}

				// setting cost around turning point depending on label
				for(auto n_id :mesh->nodes()){
					Node n = mesh->get<Node>(n_id);

					if ((*node_turningPoint)[n_id] == 1){
						math::Point turning_point = n.center();
						std::vector<TCellID> faces = n.getIDs<Face>();
						for (auto f : non_monotone_faces){
							math::Point face_center = (mesh->get<Face>(f)).center();
							if (turning_point.distance(face_center) < ABoundingBoxRatio*BOUNDING_BOX_LENTH){
								dataCost[linker_ID2faces[f]*6 + label] = dataCost_static[linker_ID2faces[f]*6 + label] + delta;
								dataCost[linker_ID2faces[f]*6 + label +1] = dataCost_static[linker_ID2faces[f]*6 + label + 1] + delta;
								(*face_constrained)[f]=1;
							}
						}

					}
				}

				graphCut.setDataCost(dataCost.data());

				// running gc
				try { //graphCut.expansion(300);
					graphCut.swap(300);
				} catch (const GCException& e) {std::cout << e.message << '\n'; }

				// checking variables

				for(auto f_id : mesh->faces()){
					if (mesh->isMarked<Face>(f_id,marked_face)) {
						(*face_chosenNormal_ID)[f_id] = graphCut.whatLabel(linker_ID2faces[f_id]);
						(*face_chosenNormal)[f_id] = vect[(*face_chosenNormal_ID)[f_id]];
					}
				}
                this->chartsGeneration(false);
				this->compute_turning_point(false);

				edges_of_each_charts = std::vector<std::vector<Triplet_charts_edge>>();
				for (auto c : charts){
					edges_of_each_charts.push_back(std::vector<Triplet_charts_edge>());
				}
				non_monotone = std::vector<bool>(charts.size(),false);
				new_nb_turning_point = 0;

				for (auto TEDGE : charts_edges){
					new_nb_turning_point += TEDGE.turningPoints.size();
				}
				std::cout << new_nb_turning_point << ", ";

				if (best_so_far > new_nb_turning_point){
					best_so_far = new_nb_turning_point;
					for(auto f_id : mesh->faces()){
						if (mesh->isMarked<Face>(f_id,marked_face)) {
							labelling_save[linker_ID2faces[f_id]] = (*face_chosenNormal_ID)[f_id];
						}
					}
				}
			}
			std::cout << '\n';
		}

		for(auto f_id : mesh->faces()){
		    if (mesh->isMarked<Face>(f_id,marked_face)) {
				 (*face_chosenNormal_ID)[f_id] = labelling_save[linker_ID2faces[f_id]];
				 (*face_chosenNormal)[f_id] = vect[(*face_chosenNormal_ID)[f_id]];
			}
		}

		// recalculating edge and turning points
        this->chartsGeneration();
		this->compute_turning_point();
		// assigning charts_edges to charts
		There_is_still_turning_points = false;
		edges_of_each_charts = std::vector<std::vector<Triplet_charts_edge>>();
		for (auto c : charts){
			edges_of_each_charts.push_back(std::vector<Triplet_charts_edge>());
		}
		non_monotone = std::vector<bool>(charts.size(),false);
		nb_turning_point = 0;
		// assigning charts_edges to charts
		for (auto TEDGE : charts_edges){
			edges_of_each_charts[TEDGE.chart1].push_back(TEDGE);
			edges_of_each_charts[TEDGE.chart2].push_back(TEDGE);
			nb_turning_point += TEDGE.turningPoints.size();
			if (!TEDGE.turningPoints.empty()){
				non_monotone[TEDGE.chart1] = true;
				non_monotone[TEDGE.chart2] = true;
				There_is_still_turning_points = true;
			}
		}
		mesh->unmarkAll<Face>(marked_face);

		IGMeshIOService service(mesh);
		VTKWriter w(&service);
        w.setCellOptions(N|F);
        w.setDataOptions(N|F);
		w.write("turningPoint"  + std::to_string(writing_advencement++));
		std::cout << "writing : " << writing_advencement << '\n';
	}
	
	mesh->freeMark<Face>(marked_face);
}
/*----------------------------------------------------------------------------*/
void PolycubeToolbox::graphCutRun(){
    
    this->graphCutBasedFaceSegmentation();

    this->graphCutHillClimbing(0, 0.005, 0.15);

    this->gregson2011Init();
    
    this->gregson2011SetVariableComputationInIteration(false);
    
	for (size_t i = 0; i < 3; i++) {
		for(auto f_id: mesh->faces()){
			Face f = mesh->get<Face>(f_id);
			if (mesh->isMarked(f,mark_bnd_faces)) {
				(*face_normals)[f.id()] = f.normal();
			}
		}
        this->gregson2011Iteration();
	}
	this->compute_turning_point();
    this->projectOnPolycube();
}
/*----------------------------------------------------------------------------*/
Mesh PolycubeToolbox::getMesh(){
	return *(this->mesh);
}
/*----------------------------------------------------------------------------*/
// use as get_3rd_vector[i*6+j] for the third axis of i and j
const math::Vector get3DVector[36] = {{0, 0, 0},
	{0,0,0},
	{0,0,1},
	{0,0,1},
	{0,1,0},
	{0,1,0},
	{0,0,0},
	{0,0,0},
	{0,0,1},
	{0,0,1},
	{0,1,0},
	{0,1,0},
	{0,0,1},
	{0,0,1},
	{0,0,0},
	{0,0,0},
	{1,0,0},
	{1,0,0},
	{0,0,1},
	{0,0,1},
	{0,0,0},
	{0,0,0},
	{1,0,0},
	{1,0,0},
	{0,1,0},
	{0,1,0},
	{1,0,0},
	{1,0,0},
	{0,0,0},
	{0,0,0},
	{0,1,0},
	{0,1,0},
	{1,0,0},
	{1,0,0},
	{0,0,0},
	{0,0,0}
};

struct extra_data_cpq {
	std::vector<int> convert_to_actual_label;
	std::vector<std::map<int, math::Vector>> edge_vector;
};

double Cpq_iterative_smoothing(int p, int q, int label_p, int label_q,void *extraData){

	double sigma = 1.2;
	double cons =3;
	if (label_q == label_p){
		return 0;
	} else{

		double npnq = (( (* (extra_data_cpq *) extraData).edge_vector[p][q]
		).dot(
                get3DVector[(* (extra_data_cpq *) extraData).convert_to_actual_label[label_p] * 6
                            + (* (extra_data_cpq *) extraData).convert_to_actual_label[label_q]]
		) -1)/sigma;
		return cons*exp(-0.5*npnq*npnq);
	}
}

inline math::Vector PolycubeToolbox::convert_edge_to_vector(Edge e){
	std::vector<TCellID> nodes = e.getIDs<Node>();
	Node n1 = mesh->get<Node>(nodes[0]), n2 = mesh->get<Node>(nodes[1]);
	return math::Vector({n1.X() - n2.X(), n1.Y() - n2.Y(), n1.Z() - n2.Z()});
}


/*----------------------------------------------------------------------------*/
void PolycubeToolbox::graph_cut_labelling_update(){

	//updating around edges
	std::cout << charts_edges.size() << '\n';
	for (auto TEDGE : charts_edges){
		extra_data_cpq *extraData = new extra_data_cpq();
		int label1 = (*face_chosenNormal_ID)[charts[TEDGE.chart1][0]];
		int label2 = (*face_chosenNormal_ID)[charts[TEDGE.chart2][0]];
		(* extraData).convert_to_actual_label.push_back(label1);
		(* extraData).convert_to_actual_label.push_back(label2);

		std::vector<double> dataCost;
		std::vector<TCellID> gc2ID;
		std::map<TCellID, int> ID2gc;
		int num=0;

		for (auto f : charts[TEDGE.chart1]){

			gc2ID.push_back(f);
			ID2gc.insert(std::pair<TCellID, int>(f,num)); num++;
			(* extraData).edge_vector.push_back( std::map<int, math::Vector>() );

			bool is_on_part_bnd = false;
			for (auto e : (mesh->get<Face>(f)).getIDs<Edge>() ){
				for (auto other_f : (mesh->get<Edge>(e)).getIDs<Face>() ){
					if (other_f !=f && mesh->isMarked<Face>(other_f,mark_bnd_faces)){
						is_on_part_bnd = is_on_part_bnd || !((*face_chosenNormal_ID)[other_f] == label1 || (*face_chosenNormal_ID)[other_f] == label2 );
					}
				}
			}

			if (is_on_part_bnd) {
				dataCost.push_back(0);
				dataCost.push_back(50);
			} else {
				double nts=  (( mesh->get<Face>(f) ).normal().dot( vect[label1] ) - 1) /0.2;
				dataCost.push_back(1-exp(-0.5*nts*nts));
				nts=  (( mesh->get<Face>(f) ).normal().dot( vect[label2] ) - 1) /0.2;
				dataCost.push_back(1-exp(-0.5*nts*nts));
			}
		}
		for (auto f : charts[TEDGE.chart2]){
			gc2ID.push_back(f);
			ID2gc.insert(std::pair<TCellID, int>(f,num)); num++;
			(* extraData).edge_vector.push_back( std::map<int, math::Vector>() );

			bool is_on_part_bnd = false;
			for (auto e : (mesh->get<Face>(f)).getIDs<Edge>() ){
				for (auto other_f : (mesh->get<Edge>(e)).getIDs<Face>() ){
					if (other_f !=f && mesh->isMarked<Face>(other_f,mark_bnd_faces)){
						is_on_part_bnd = is_on_part_bnd || !((*face_chosenNormal_ID)[other_f] == label1 || (*face_chosenNormal_ID)[other_f] == label2 );
					}
				}
			}

			if (is_on_part_bnd) {
				dataCost.push_back(50);
				dataCost.push_back(0);
			} else {
				double nts=  (( mesh->get<Face>(f) ).normal().dot( vect[label1] ) - 1) /0.2;
				dataCost.push_back(1-exp(-0.5*nts*nts));
				nts=  (( mesh->get<Face>(f) ).normal().dot( vect[label2] ) - 1) /0.2;
				dataCost.push_back(1-exp(-0.5*nts*nts));
			}
		}

		GCoptimizationGeneralGraph graphCut(gc2ID.size(),2);


		for (size_t i = 0; i < gc2ID.size(); i++) {
			std::vector<TCellID> edges = ( mesh->get<Face>(gc2ID[i]) ).getIDs<Edge>();
			for (auto e : edges) {
				Edge edge = mesh->get<Edge>(e);
				std::vector<TCellID> faces = ( edge ).getIDs<Face>();
				for (auto f : faces) {
					if (((*face_chosenNormal_ID)[f] == label1 || (*face_chosenNormal_ID)[f] == label2) && gc2ID[i] != f) {
						graphCut.setNeighbors(i, ID2gc[f]);
						(* extraData).edge_vector[i].insert( std::pair<int,math::Vector>(ID2gc[f],convert_edge_to_vector(edge)) );
					}
				}
			}
		}

		graphCut.setDataCost(dataCost.data());

		graphCut.setSmoothCost(&Cpq_iterative_smoothing, (void *)extraData);

		try { //graphCut.expansion(300);

			graphCut.swap(1000);
		} catch (const GCException& e) {std::cout << e.message << '\n'; }

		for (size_t i = 0; i < gc2ID.size(); i++) {
			(*face_chosenNormal_ID)[gc2ID[i]] = graphCut.whatLabel(i);
			(*face_chosenNormal)[gc2ID[i]] = vect[(*face_chosenNormal_ID)[gc2ID[i]]];
		}

	}



	// updating around corners
	std::map<TCellID, bool> corner_already_modified;
	for (size_t i = 0; i < mesh->getNbNodes(); i++) {
		if ((*node_chartCorner)[i] == 1 ){
			corner_already_modified.insert(std::pair<TCellID,bool>(i,false));
		}
	}

	for (auto TEDGE : charts_edges){
		if (!corner_already_modified[TEDGE.corner1]){

			corner_already_modified[TEDGE.corner1] = true;
		}
		if (!corner_already_modified[TEDGE.corner2]){

			corner_already_modified[TEDGE.corner2] = true;
		}
	}


}


/*----------------------------------------------------------------------------*/
void PolycubeToolbox::set_angle_determinent(double det){
	this->ANGLE_DETERMINENT = det;
}

/*----------------------------------------------------------------------------*/
void PolycubeToolbox::test_fonction(){
	mesh->deleteVariable(GMDS_NODE,"sharp_edge");
	gmds::Variable<int>* 	sharp_edge = mesh->newVariable<int,GMDS_NODE>("sharp_edge");
	for(auto e_id : mesh->edges()){
		Edge e  = mesh->get<Edge>(e_id);
		if (mesh->isMarked(e,mark_sharp_edges)){
			std::vector<TCellID> nodes = e.getIDs<Node>();
			for (auto n : nodes){
				(*sharp_edge)[n]=1;
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<TCellID>> PolycubeToolbox::get_faces_of_block(){
	std::vector<std::vector<TCellID>> corner_of_faces;
	for (auto c : charts){
		corner_of_faces.push_back(std::vector<TCellID>());
	}
	for (auto TEDGE : charts_edges){
		corner_of_faces[TEDGE.chart1].push_back(TEDGE.corner1);
		corner_of_faces[TEDGE.chart1].push_back(TEDGE.corner2);
		corner_of_faces[TEDGE.chart2].push_back(TEDGE.corner1);
		corner_of_faces[TEDGE.chart2].push_back(TEDGE.corner2);
	}
	return corner_of_faces;
}
