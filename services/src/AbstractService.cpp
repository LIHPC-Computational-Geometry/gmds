/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/9/19.
//
/*----------------------------------------------------------------------------*/
#include <gmds/utils/Exception.h>
#include "gmds/services/AbstractService.h"
#include "gmds/services/Property.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
void AbstractService::executeAfterChecking() {
    if(!checkInput())
        throw GMDSException("Invalid input data");

    execute();
}
/*----------------------------------------------------------------------------*/
void AbstractService::addInput(const AbstractData *AData) {
    m_input.insert(AData);
}
/*----------------------------------------------------------------------------*/
bool AbstractService::addConstraint(const gmds::AbstractData *AData,
                                    const gmds::Property *AProp)
{
    //if AData is not an input or an output, we don't add this constraint
    if(m_input.find(AData) == m_input.end() &&
       m_output.find(AData)== m_output.end())
        return false;

    m_constraints[AData].push_back(AProp);
    return true;
}
/*----------------------------------------------------------------------------*/
bool AbstractService::checkInput() {

    for(auto d:m_input){
        std::vector<const Property*> props = m_constraints[d];
        for(auto p:props){
            if(!p->isValid(d)){
                return false;
            }
        }
    }
    return true;
}