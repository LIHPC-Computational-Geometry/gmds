/*----------------------------------------------------------------------------*/
/*
 *  Parameters.cpp
 *
 *  Created on: April 12, 2016
 *  Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <dictionary.h>
#include <iniparser.h>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/Parameters.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Parameters::Parameters()
:m_initialized(false)
{}
/*----------------------------------------------------------------------------*/
Parameters::~Parameters()
{}
/*----------------------------------------------------------------------------*/
bool Parameters::add(const std::string& ASection,
                     const std::string& AName,
                     const ETypeParam AType)
{
    for(auto e:m_entries){
        if(e.section==ASection && e.name==AName){
            return false;
        }
    }
    
    Entry e;
    e.section = ASection;
    e.name    = AName;
    e.type    = AType;
    
    m_entries.push_back(e);
    return true;
}
/*----------------------------------------------------------------------------*/
std::vector<Parameters::Entry> Parameters::getEntries() const
{
    return m_entries;
}
/*----------------------------------------------------------------------------*/
int Parameters::
find(const std::string& ASection, const std::string& AName)
{
    for(auto i=0; i<m_entries.size();i++){
        Entry e = m_entries[i];
        if(e.section==ASection && e.name==AName){
            //we found the right entry
            return i;
        }
    }
    return -1;
}
/*----------------------------------------------------------------------------*/
bool Parameters::
get(const std::string& ASection, const std::string& AName, double& AOut)
{
    if(!m_initialized){
        //the paramaters have no values
        return false;
    }
    int i = find(ASection,AName);
    if(i==-1){
        //there is no parameter with this section and name
        return false;
    }
    if(m_entries[i].type!=DOUBLE_P){
        //Wrong type
        return false;
    }
    //right type so we can return the right value
    AOut=std::stod(m_values[i]);
    return true;
}/*----------------------------------------------------------------------------*/
bool Parameters::
get(const std::string& ASection, const std::string& AName, int& AOut)
{
    if(!m_initialized){
        //the paramaters have no values
        return false;
    }
    int i = find(ASection,AName);
    if(i==-1){
        //there is no parameter with this section and name
        return false;
    }
    if(m_entries[i].type!=INT_P){
        //Wrong type
        return false;
    }
    //right type so we can return the right value
    AOut=std::stoi(m_values[i]);
    return true;
}
/*----------------------------------------------------------------------------*/
bool Parameters::
get(const std::string& ASection, const std::string& AName, std::string& AOut)
{
    if(!m_initialized){
        //the paramaters have no values
        return false;
    }
    int i = find(ASection,AName);
    if(i==-1){
        //there is no parameter with this section and name
        return false;
    }
    if(m_entries[i].type!=STRING_P){
        //Wrong type
        return false;
    }
    //right type so we can return the right value
    AOut=m_values[i];
    return true;
}
/*----------------------------------------------------------------------------*/
bool Parameters::
get(const std::string& ASection, const std::string& AName, bool& AOut)
{
    if(!m_initialized){
        //the paramaters have no values
        return false;
    }
    int i = find(ASection,AName);
    if(i==-1){
        //there is no parameter with this section and name
        return false;
    }
    if(m_entries[i].type!=BOOL_P){
        //Wrong type
        return false;
    }
    //right type so we can return the right value
    std::string str_val = m_values[i];
    if(str_val=="true"){
        AOut=true;
    }
    else if(str_val=="false"){
        AOut=false;
    }
    else{
        //Wrong value
        return false;
    }
    return true;
}
/*----------------------------------------------------------------------------*/
std::vector<std::string> Parameters::parseIni(const std::string& AFileName)
{
    std::vector<std::string> wrong;
    if(m_entries.empty()){
        //nothing to do
        return wrong;
    }
    m_values.resize(m_entries.size());
    
    //We create a dictionnary object to parse the read file
    dictionary* dico = iniparser_load(AFileName.c_str());
    
    //For each entry, we look a right type value in the file
    for(auto i=0; i<m_entries.size();i++){
        Entry ei = m_entries[i];
        std::string param_s = ei.section+":"+ei.name;
        const char* param_ch = param_s.c_str();
        // ----------- INT VALUE -----------------------------
        if(ei.type==INT_P){
            int not_found = -1111111;
            int v = iniparser_getint(dico,(char *)param_ch,not_found);
            if(v==not_found){
                wrong.push_back(param_s);
            }
            else{
                m_values[i]=std::to_string(v);
            }
        }
        // ----------- DOUBLE VALUE -----------------------------
        else if(ei.type==DOUBLE_P){
            double not_found = -1111.11111;
            double v = iniparser_getdouble(dico,(char *)param_ch,not_found);
            if(v==not_found){
                wrong.push_back(param_s);
            }
            else{
                m_values[i]=std::to_string(v);
            }
        }  // ----------- BOOL VALUE -----------------------------
        else if(ei.type==BOOL_P){
            int not_found = 10;
            int v = iniparser_getboolean(dico,(char *)param_ch,not_found);
            if(v==not_found){
                wrong.push_back(param_s);
            }
            else if (v==1){
                m_values[i]="true";
            }
            else {
                m_values[i]="false";
            }
        }
        // ----------- STRING VALUE -----------------------------
        else {
            const char* v = iniparser_getstring(dico,param_ch, NULL);
            if(v==NULL){
                wrong.push_back(param_s);
            }
            else{
                std::string v_s(v);
                m_values[i]=v_s;
            }
        }
    }
    
    //All the value have been found and initialized
    m_initialized=true;
    //We free the dictionary
    iniparser_freedict(dico);
    
    return wrong;
}
/*----------------------------------------------------------------------------*/
