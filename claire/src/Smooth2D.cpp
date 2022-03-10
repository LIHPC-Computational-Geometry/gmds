/*------------------------------------------------------------------------*/
//
// Created by Claire Roche on 21/10/2021.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/Smooth2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
Smooth2D::Smooth2D(Mesh *AMesh,
                   const Variable<int>* AVarBnd,
                   const int ANbIterations) {
	m_mesh = AMesh;
	m_node_constrained = AVarBnd;
	m_nb_max_iterations = ANbIterations;
}
/*------------------------------------------------------------------------*/
void Smooth2D::setNbIterations(const int ANbIterations)
{
	m_nb_max_iterations=ANbIterations;

}
/*------------------------------------------------------------------------*/
Smooth2D::STATUS Smooth2D::execute()
{
    //we traverse all the nodes of the whole mesh to get the nodes that
    //are not constrained, i.e are so free
	for(auto n_id:m_mesh->nodes()){
		if(m_node_constrained->value(n_id)!=1){
			m_free_nodes.push_back(n_id);
		}
	}
    //for all free nodes, we build and store stencils locally
    buildStencils();

	return Smooth2D::SUCCESS;
}
/*------------------------------------------------------------------------*/
void Smooth2D::buildStencils() {
    for(auto n_id:m_free_nodes){
        Node current_node = m_mesh->get<Node>(n_id);
        std::vector<Face> current_faces = current_node.get<Face>();
        // FIRST FACE
        Face f0 = current_faces[0];
        Node f0_n0;
        Node f0_n1;
        f0.getAdjacentNodes(current_node,f0_n0,f0_n1);
        stencil current_stencil;
        current_stencil.val[1][1] = current_node.id();
        current_stencil.val[1][0] = f0_n0.id();
        current_stencil.val[0][1] = f0_n1.id();
        std::vector<TCellID> f_nodes=f0.getIDs<Node>();
        for(auto f_n_id:f_nodes){
            if(f_n_id!=current_node.id() &&
               f_n_id!=f0_n0.id() &&
               f_n_id!=f0_n1.id()) {
                current_stencil.val[0][0] = f_n_id;
            }
        }
        //SECOND FACE
        std::vector<TCellID> faces_11_10 = m_mesh->getCommonFaces(current_node,f0_n0);
        TCellID f1_id = NullID;
        if(faces_11_10[0]==f0.id()){
            f1_id=faces_11_10[1];
        }
        else{
            f1_id=faces_11_10[0];
        }
        Face f1 = m_mesh->get<Face>(f1_id);
        Node f1_n0;
        Node f1_n1;
        f1.getAdjacentNodes(current_node,f1_n0,f1_n1);
        if(f1_n0.id()==current_stencil.val[1][0])
            current_stencil.val[2][1]=f1_n1.id();
        else
            current_stencil.val[2][1]=f1_n0.id();

        f_nodes=f1.getIDs<Node>();
        for(auto f_n_id:f_nodes){
            if(f_n_id!=current_node.id() &&
               f_n_id!=f1_n0.id() &&
               f_n_id!=f1_n1.id()) {
                current_stencil.val[2][0] = f_n_id;
            }
        }
        //THIRD FACE
        std::vector<TCellID> faces_11_21 =
                m_mesh->getCommonFaces(current_node,
                                       m_mesh->get<Node>(current_stencil.val[2][1]));
        TCellID f2_id = NullID;
        if(faces_11_21[0]==f1.id()){
            f2_id=faces_11_21[1];
        }
        else{
            f2_id=faces_11_21[0];
        }
        Face f2 = m_mesh->get<Face>(f2_id);
        Node f2_n0;
        Node f2_n1;
        f2.getAdjacentNodes(current_node,f2_n0,f2_n1);

        if(f2_n0.id()==current_stencil.val[2][1])
            current_stencil.val[1][2]=f2_n1.id();
        else
            current_stencil.val[1][2]=f2_n0.id();

        f_nodes=f2.getIDs<Node>();
        for(auto f_n_id:f_nodes){
            if(f_n_id!=current_node.id() &&
               f_n_id!=f2_n0.id() &&
               f_n_id!=f2_n1.id()) {
                current_stencil.val[2][2] = f_n_id;
            }
        }

        //FOURTH FACE
        std::vector<TCellID> faces_11_12 =
                m_mesh->getCommonFaces(current_node,
                                       m_mesh->get<Node>(current_stencil.val[1][2]));
        TCellID f3_id = NullID;
        if(faces_11_12[0]==f2.id()){
            f3_id=faces_11_12[1];
        }
        else{
            f3_id=faces_11_12[0];
        }
        Face f3 = m_mesh->get<Face>(f3_id);

        f_nodes=f3.getIDs<Node>();
        for(auto f_n_id:f_nodes){
            if(f_n_id!=current_node.id() &&
               f_n_id!=current_stencil.val[1][2] &&
               f_n_id!=current_stencil.val[0][1]) {
                current_stencil.val[0][2] = f_n_id;
            }
        }

        m_stencil.insert(std::make_pair(n_id, current_stencil));
        std::cout<<"Node "<<n_id<<std::endl;
        for(auto i=0;i<3;i++){
            for(auto j=0;j<3;j++){
                std::cout << m_stencil[n_id].val[i][j] << " ";
            }
            std::cout<<std::endl;
        }
    }
}
/*------------------------------------------------------------------------*/
