/*----------------------------------------------------------------------------*/
/*
 * VariableManager.cpp
 *
 *  Created on: 26 juil. 2010
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/utils/VariableManager.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
VariableManager::VariableManager()
{
    m_variables.clear();
}
/*----------------------------------------------------------------------------*/
VariableManager::~VariableManager()
{
	for(unsigned int k=0;k<m_variables.size();k++)
	{
		VariableItf *v =m_variables[k];
		delete v;
	}
}
/*----------------------------------------------------------------------------*/
void VariableManager::deleteVariable(const std::string& AName){
	for(unsigned int k=0;k<m_variables.size();k++)
	{
		if (m_variables[k]->getName()==AName){
			VariableItf *v =m_variables[k];
			if(k!=m_variables.size()-1)
				m_variables[k] = m_variables.back();
			m_variables.pop_back();
			delete v;
		}
	}
}
/*----------------------------------------------------------------------------*/
void VariableManager::deleteVariable(VariableItf* AVar){
    for(unsigned int k=0;k<m_variables.size();k++)
	{
		if (m_variables[k]==AVar){
			VariableItf *v =m_variables[k];
			if(k!=m_variables.size()-1)
				m_variables[k] = m_variables.back();
			m_variables.pop_back();
			delete v;
		}
	}
}
/*----------------------------------------------------------------------------*/
void VariableManager::addEntry(const int i){
	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->addEntry(i);
}
/*----------------------------------------------------------------------------*/
void VariableManager::removeEntry(const int i){
	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->removeEntry(i);
}
/*----------------------------------------------------------------------------*/
void VariableManager::setDomain(const int size){
	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->setDomain(size);
}
/*----------------------------------------------------------------------------*/
void VariableManager::clearVariables(){
	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->clear();
}
/*----------------------------------------------------------------------------*/
void VariableManager::compact(){
	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->compact();
}
/*----------------------------------------------------------------------------*/
void VariableManager::serialize(std::ostream& AStr) {
	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->serialize(AStr);
}
/*----------------------------------------------------------------------------*/
void VariableManager::unserialize(std::istream& AStr){

	for(unsigned int k=0;k<m_variables.size();k++)
		m_variables[k]->unserialize(AStr);
}
/*----------------------------------------------------------------------------*/
int VariableManager::getNbVariables() const{
	return m_variables.size();
}
/*----------------------------------------------------------------------------*/
VariableItf* VariableManager::getVariable(const TInt i){
	if(i<0 || i>=(int)m_variables.size())
		throw GMDSException("VariableManager::getVariable - Out of the bounds");
	return m_variables[i];
}
/*----------------------------------------------------------------------------*/
bool VariableManager::doesVariableExist(const std::string& AName){
	unsigned int k=0;
	for(;k<m_variables.size();k++)
		if (m_variables[k]->getName()==AName)
			return true;

	return false;
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
