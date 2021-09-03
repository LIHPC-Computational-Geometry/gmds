/*----------------------------------------------------------------------------*/
/*
 * VariableManager.h
 *
 *  Created on: 26 juil. 2010
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_VARIABLEMANAGER_H_
#define GMDS_VARIABLEMANAGER_H_
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include "Variable.h"
#include "CommonTypes.h"
#include "Exception.h"
#include "GMDSUtils_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
	/*----------------------------------------------------------------------------*/
	/** \class VariableManager
	 *  \brief Handle the creation and update of a collection of variables. A
	 *  	   variable is defined as a set of discrete values associated to a key.
	 *  	   Few holes are in the key numerotation.
	 */
	class GMDSUtils_API VariableManager{

	public:
		/*------------------------------------------------------------------------*/
		/** \brief  Constructor.
		*/
		VariableManager();

		/*------------------------------------------------------------------------*/
		/** \brief  Destructor.
		*/
		~VariableManager();

		/*------------------------------------------------------------------------*/
		/** \brief  creation of a variable allocated in the stack. The domain of the
		 * 			variable is initialized to [0,initSize].
		 */
		template<typename T> Variable<T>* newVariable(const std::string& AName,
			const int initSize = 2, const std::vector<int>* ref = 0);

		/*------------------------------------------------------------------------*/
		/** \brief  Returns whether a variable exists.
		*/
		bool doesVariableExist(const std::string& AName);

		/*------------------------------------------------------------------------*/
		/** \brief  Access to a variable.
		*/
		template<typename T> Variable<T>* getVariable(const std::string& AName);

		/*------------------------------------------------------------------------*/
		/** \brief  suppression of a variable. the memory used in the stack is free.
		*/
		void deleteVariable(const std::string& AName);

		/*------------------------------------------------------------------------*/
		/** \brief  suppression of a variable. the memory used in the stack is free.
		*/
		void deleteVariable(VariableItf* AVar);

		/*------------------------------------------------------------------------*/
		/** \brief  Add a default entry to all the variables. Each variable has the
		 * 			responsability to initialize the corresponding value.
		 */
		void addEntry(const int i);

		/*------------------------------------------------------------------------*/
		/** \brief  Remove the i th entry of all the managed variables.
		 */
		void removeEntry(const int i);

		/*------------------------------------------------------------------------*/
		/** \brief  set the domain of all the variables to size.
		 */
		void setDomain(const int size);

		/*------------------------------------------------------------------------*/
		/** \brief  set the domain of all the variables to size.
		 */
		int getNbVariables() const;

		/*------------------------------------------------------------------------*/
		/** \brief  get access to the i^th variable in a abstract form
		 */
		VariableItf* getVariable(const TInt i);

		/*------------------------------------------------------------------------*/
		/** \brief indicates if variables are attached
		 */
		bool empty() const;

		/*------------------------------------------------------------------------*/
		/** \brief clear all the variable
		 */
		void clearVariables();

		/*------------------------------------------------------------------------*/
		/** \brief compact all the variables
		 */
		void compact();

		/*------------------------------------------------------------------------*/
		/** \brief serialize (*this) in AStr
		 *
		 * \param AStr an output streammap
		 */
		void serialize(std::ostream& AStr);

		/*------------------------------------------------------------------------*/
		/** \brief unserialize (*this) from AStr
		 *
		 * \param AStr an input stream
		 */
		void unserialize(std::istream& AStr);

		/*------------------------------------------------------------------------*/
		/** \brief  get the list of variables
		*/
		std::vector<VariableItf*> getAllVariables(){ return m_variables; }

	private:

		std::vector<VariableItf*>  m_variables;
	};
	/*----------------------------------------------------------------------------*/
	template<typename T>
	Variable<T>* VariableManager::newVariable(const std::string& AName,
		const int initSize, const std::vector<int>* ref){

		for (unsigned int k = 0; k < m_variables.size(); k++)
            if (m_variables[k]->getName() == AName){
                std::string mess= "Impossible to create a variable "+m_variables[k]->getName()+": name already used";
			throw GMDSException(mess);
            }

		Variable<T>* v = new Variable<T>(AName);
		m_variables.push_back(v);

		if (ref == 0)
			v->setDomain(initSize + 1);
		else
			v->setDomainWithDefault(initSize + 1, *ref);

		return v;
	}
	/*----------------------------------------------------------------------------*/
	template<typename T>
	Variable<T>* VariableManager::getVariable(const std::string& AName){
		unsigned int k = 0;
		for (; k < m_variables.size(); k++)
		if (m_variables[k]->getName() == AName)
			return dynamic_cast<Variable<T>*>(m_variables[k]);
        std::string mess= "No variable named "+AName;
		throw GMDSException(mess);
	}
	/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_VARIABLEMANAGER_H_ */
/*----------------------------------------------------------------------------*/
