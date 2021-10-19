/*----------------------------------------------------------------------------*/
/*
 * VariableManager.h
 *
 *  Created on: 03/10/2017
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_VARIABLEMANAGER_H_
#define KMDS_VARIABLEMANAGER_H_
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <KM/Utils/Exception.h>
#include <KM/Utils/KTypes.h>
#include <KM/Utils/Variable.h>
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
/** \class VariableManager
 *  \brief Handle the creation and update of a collection of variables. A
 *  	   variable is defined as a set of discrete values associated to a key.
 *  	   Few holes are in the key numerotation.
 */
class EXPORT_KMDS VariableManager
{
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
        template <typename T>
        Variable<T>* createVariable(const T ADefaultValue, const std::string& AName, const int initSize = 2);

        /*------------------------------------------------------------------------*/
        /** \brief  Returns whether a variable exists.
        */
        bool doesVariableExist(const std::string& AName);

        /*------------------------------------------------------------------------*/
        /** \brief  Access to a variable.
        */
        template <typename T>
        EXPORT_KMDS Variable<T>* getVariable(const std::string& AName);

        /*------------------------------------------------------------------------*/
        /** \brief  suppression of a variable. the memory used in the stack is free.
        */
        void deleteVariable(const std::string& AName);

        /*------------------------------------------------------------------------*/
        /** \brief  suppression of a variable. the memory used in the stack is free.
        */
        void deleteVariable(VariableItf* AVar);

        /*------------------------------------------------------------------------*/
        /** \brief  set the domain of all the variables to \p ASize.
         */
        void resize(const int ASize);

        /*------------------------------------------------------------------------*/
        /** \brief  get the domain size of all the variables
                 */
        int getSize() const;

        /*------------------------------------------------------------------------*/
        /** \brief  get access to the i^th variable in a abstract form
         */
        VariableItf* getVariable(const TInt32 i);

        /*------------------------------------------------------------------------*/
        /** \brief  Initialize all the variables related to index \p AI. Default
         *          value is known by every single variable
         */
        void initializeVariables(const TCellID AI);
        /*------------------------------------------------------------------------*/
        /** \brief  Initialize all the variables related to indices \p AI to
         *          \p AI + \p ANb. Default value is known by every single variable
         */
        void initializeVariables(const TCellID AI, const TInt32 ANb);

        /*------------------------------------------------------------------------*/
        /** \brief indicates if variables are attached
         */
        bool empty() const;

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
        /** \brief  get the the number of variables
        */
        int getNbVariables() const;
        /*------------------------------------------------------------------------*/
        /** \brief  get the list of variables
        */
        std::vector<VariableItf*> getAllVariables();

 private:
        std::vector<VariableItf*> m_variables;
};
/*----------------------------------------------------------------------------*/
template <typename T>
Variable<T>*
VariableManager::createVariable(const T ADefaultValue, const std::string& AName, const int initSize)
{
        for (auto v : m_variables)
                if (v->getName() == AName) {
                        std::string mess = "Impossible to create a variable " + AName + ": name already used";
                        throw KException(mess);
                }

        Variable<T>* v = new Variable<T>(ADefaultValue, AName);
        m_variables.push_back(v);

        v->resize(initSize + 1);

        return v;
}
/*----------------------------------------------------------------------------*/
template <typename T>
Variable<T>*
VariableManager::getVariable(const std::string& AName)
{
        for (auto v : m_variables) {
                if (v->getName() == AName)
                        return dynamic_cast<Variable<T>*>(v);
        }
        std::string mess = "No variable named " + AName;
        throw KException(mess);
}
/*----------------------------------------------------------------------------*/
}  // namespace kmds
/*----------------------------------------------------------------------------*/
#endif /* KMDS_VARIABLEMANAGER_H_ */
/*----------------------------------------------------------------------------*/
