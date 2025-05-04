//
// Created by calderans on 02/03/23.
//
/*----------------------------------------------------------------------------*/
#include "gmds/io/MdlReader.h"
#include <regex>
#include <utility>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
MdlReader::MdlReader(Mesh &AMesh,std::string AString)
  :m_mesh(AMesh), m_lengthUnit(1.), m_stream(nullptr), m_name2find(std::move(AString))
{}
/*----------------------------------------------------------------------------*/
MdlReader::~MdlReader()
{
	delete m_stream;
}
/*----------------------------------------------------------------------------*/
void MdlReader::read(const std::string &AFileName){

	std::cout<<"Start reading"<<std::endl;
	//===============================================================
	// First, we check if we can read the file
	m_stream= new std::ifstream(AFileName.c_str(),std::ios::in);
	if(!m_stream->is_open()){
		std::string mess ="Impossible to read file "+AFileName;
		throw GMDSException(mess);
	}

	std::string supp_name;
	while (moveStreamOntoSupport(supp_name)) {
		CellGroup<Node> *g;
		try {
			g = m_mesh.getGroup<Node>(supp_name);
		}
		catch (GMDSException &) {
			g = m_mesh.newGroup<Node>(supp_name);
		}

		std::string current_word = "";
		*m_stream >> current_word;
		*m_stream >> current_word;
		double X, Y;

		// nodes of the support are stocked in a vector instead of the group
		// in case a group is composed of multiple supports not linked
		std::vector<TCellID> nodes;

		while (std::string::npos == current_word.find(')')) {
			X = std::stod(current_word);
			*m_stream >> current_word;
			Y = std::stod(current_word);
			*m_stream >> current_word;

			Node n = m_mesh.newNode(X, Y);
			nodes.push_back(n.id());
		}

		g->add(nodes[0]);
		for (int i = 1; i < nodes.size(); i++) {
			m_mesh.newEdge(nodes[i - 1], nodes[i]);
			g->add(nodes[i]);
		}
	}
}
/*----------------------------------------------------------------------------*/
bool MdlReader::moveStreamOntoSupport(std::string &AString){

	bool found = false;
	std::string line;
	while (!found && std::getline(*m_stream,line))  {
		// some of our input files contain a whitespace at the end of the line after Support
		found = std::regex_match(line, std::regex(m_name2find+"(.*)(Support)([' ']?)"));
	}

	AString = line.substr(0, line.find_first_of("_ "));

	return found;
}
/*----------------------------------------------------------------------------*/
Mesh* MdlReader::getMesh(){
	return &m_mesh;
}
/*----------------------------------------------------------------------------*/
void MdlReader::createVariablesFromGroups(){
	for(int i=0; i<m_mesh.getNbGroups<Node>(); i++){
		CellGroup<Node>* group = m_mesh.getGroup<Node>(i);
		Variable<int>* group_var = m_mesh.newVariable<int,GMDS_NODE>(group->name());
		for(auto n : group->cells()){
			group_var->set(n, 1);
		}
	}
}