/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
 * VariableManager.cpp
 *
 *  Created on: 26 juil. 2010
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <KM/Utils/VariableManager.h>
/*----------------------------------------------------------------------------*/
using namespace kmds;
/*----------------------------------------------------------------------------*/
VariableManager::VariableManager()
{
        ;
}
/*----------------------------------------------------------------------------*/
VariableManager::~VariableManager()
{
        for (unsigned int k = 0; k < m_variables.size(); k++) {
                VariableItf* v = m_variables[k];
                delete v;
        }
}
/*----------------------------------------------------------------------------*/
void
VariableManager::deleteVariable(const std::string& AName)
{
        for (auto k = 0; k < m_variables.size(); k++) {
                VariableItf* v = m_variables[k];
                if (v->getName() == AName) {
                        if (k != m_variables.size() - 1)
                                m_variables[k] = m_variables.back();
                        m_variables.pop_back();
                        delete v;
                        return;
                }
        }

        throw KException("VariableManager::deleteVariable - variable not found");
}
/*----------------------------------------------------------------------------*/
void
VariableManager::deleteVariable(VariableItf* AVar)
{
        for (auto k = 0; k < m_variables.size(); k++) {
                VariableItf* v = m_variables[k];
                if (v == AVar) {
                        if (k != m_variables.size() - 1)
                                m_variables[k] = m_variables.back();
                        m_variables.pop_back();
                        delete v;
                        return;
                }
        }

        throw KException("VariableManager::deleteVariable - variable not found");
}
/*----------------------------------------------------------------------------*/
void
VariableManager::resize(const int ASize)
{
        for (auto v : m_variables)
                v->resize(ASize);
}
/*----------------------------------------------------------------------------*/
void
VariableManager::initializeVariables(const TCellID AID)
{
        for (auto v : m_variables)
                v->initialize(AID);
}
/*----------------------------------------------------------------------------*/
void
VariableManager::initializeVariables(const TCellID AID, const TInt32 ANb)
{
        for (auto v : m_variables)
                v->initialize(AID, ANb);
}
/*----------------------------------------------------------------------------*/
void
VariableManager::compact()
{
        for (auto v : m_variables)
                v->compact();
}
/*----------------------------------------------------------------------------*/
void
VariableManager::serialize(std::ostream& AStr)
{
        for (auto v : m_variables)
                v->serialize(AStr);
}
/*----------------------------------------------------------------------------*/
void
VariableManager::unserialize(std::istream& AStr)
{
        for (auto v : m_variables)
                v->unserialize(AStr);
}
/*----------------------------------------------------------------------------*/
int
VariableManager::getNbVariables() const
{
        return m_variables.size();
}
/*----------------------------------------------------------------------------*/
std::vector<VariableItf*>
VariableManager::getAllVariables()
{
        return m_variables;
}
/*----------------------------------------------------------------------------*/
VariableItf*
VariableManager::getVariable(const TInt32 i)
{
        if (i < 0 || i >= m_variables.size())
                throw KException("VariableManager::getVariable - Out of the bounds");
        return m_variables[i];
}
/*----------------------------------------------------------------------------*/
bool
VariableManager::doesVariableExist(const std::string& AName)
{
        unsigned int k = 0;
        for (; k < m_variables.size(); k++)
                if (m_variables[k]->getName() == AName)
                        return true;

        return false;
}
/*----------------------------------------------------------------------------*/
