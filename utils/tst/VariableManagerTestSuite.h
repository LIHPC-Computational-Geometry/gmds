#ifndef VARIABLE_MANAGER_TESTSUITE_H
#define VARIABLE_MANAGER_TESTSUITE_H

#include "gmds/utils/VariableManager.h"
#include "gtest/gtest.h"
#include <unit_test_config.h>

using namespace gmds;

TEST(VariableManagerTestSuite, NewVariable)
{
	VariableManager variableManager;

	Variable<int> *intVariable = variableManager.newVariable<int>("intVariable", 2);
	ASSERT_TRUE(intVariable != nullptr);
	ASSERT_EQ(intVariable->getName(), "intVariable");
	ASSERT_EQ(intVariable->getDomainSize(), 3);

	Variable<double> *doubleVariable = variableManager.newVariable<double>("doubleVariable", 3);
	ASSERT_TRUE(doubleVariable != nullptr);
	ASSERT_EQ(doubleVariable->getName(), "doubleVariable");
	ASSERT_EQ(doubleVariable->getDomainSize(), 4);
}

TEST(VariableManagerTestSuite, GetVariable)
{
	VariableManager variableManager;

	Variable<int> *intVariable = variableManager.newVariable<int>("intVariable", 2);
	Variable<double> *doubleVariable = variableManager.newVariable<double>("doubleVariable", 3);

	ASSERT_EQ(variableManager.getVariable<int>("intVariable"), intVariable);
	ASSERT_EQ(variableManager.getVariable<double>("doubleVariable"), doubleVariable);

	ASSERT_THROW(variableManager.getVariable<int>("nonExistentVariable"), GMDSException);
}

TEST(VariableManagerTestSuite, DeleteVariable)
{
	VariableManager variableManager;

	Variable<int> *intVariable = variableManager.newVariable<int>("intVariable", 2);
	Variable<double> *doubleVariable = variableManager.newVariable<double>("doubleVariable", 3);

	ASSERT_TRUE(variableManager.doesVariableExist("intVariable"));
	ASSERT_TRUE(variableManager.doesVariableExist("doubleVariable"));

	variableManager.deleteVariable("intVariable");
	ASSERT_FALSE(variableManager.doesVariableExist("intVariable"));

	variableManager.deleteVariable(doubleVariable);
	ASSERT_FALSE(variableManager.doesVariableExist("doubleVariable"));
}

TEST(VariableManagerTestSuite, AddRemoveEntry)
{
	VariableManager variableManager;

	Variable<int> *intVariable = variableManager.newVariable<int>("intVariable", 2);
	Variable<double> *doubleVariable = variableManager.newVariable<double>("doubleVariable", 3);

	variableManager.addEntry(1);
	ASSERT_EQ(intVariable->getEntriesCount(), 1);
	ASSERT_EQ(doubleVariable->getEntriesCount(), 1);

	variableManager.removeEntry(1);
	ASSERT_EQ(intVariable->getEntriesCount(), 0);
	ASSERT_EQ(doubleVariable->getEntriesCount(), 0);
}

TEST(VariableManagerTestSuite, SetDomain)
{
	VariableManager variableManager;

	Variable<int> *intVariable = variableManager.newVariable<int>("intVariable", 2);
	Variable<double> *doubleVariable = variableManager.newVariable<double>("doubleVariable", 3);

	variableManager.setDomain(5);
	ASSERT_EQ(intVariable->getDomainSize(), 6);
	ASSERT_EQ(doubleVariable->getDomainSize(), 6);
}

TEST(VariableManagerTestSuite, GetNbVariables)
{
	VariableManager variableManager;

	ASSERT_EQ(variableManager.getNbVariables(), 0);

	Variable<int> *intVariable = variableManager.newVariable<int>("intVariable", 2);
	Variable<double> *doubleVariable = variableManager.newVariable<double>("doubleVariable", 3);

	ASSERT_EQ(variableManager.getNbVariables(), 2);

	variableManager.deleteVariable(intVariable);
	ASSERT_EQ(variableManager.getNbVariables(), 1);
}

TEST(VariableManagerTestSuite, Empty)
{
	VariableManager variableManager;

	ASSERT_TRUE(variableManager.empty());

	Variable<int> *intVariable = variableManager.newVariable<int>("intVariable", 2);

	ASSERT_FALSE(variableManager.empty());

	variableManager.deleteVariable(intVariable);
	ASSERT_TRUE(variableManager.empty());
}

TEST(VariableManagerTestSuite, ClearVariables)
{
	VariableManager variableManager;

	Variable<int> *intVariable = variableManager.newVariable<int>("intVariable", 2);
	Variable<double> *doubleVariable = variableManager.newVariable<double>("doubleVariable", 3);

	ASSERT_EQ(variableManager.getNbVariables(), 2);

	variableManager.clearVariables();

	ASSERT_EQ(variableManager.getNbVariables(), 0);
}

TEST(VariableManagerTestSuite, Compact)
{
	VariableManager variableManager;

	Variable<int> *intVariable = variableManager.newVariable<int>("intVariable", 2);
	Variable<double> *doubleVariable = variableManager.newVariable<double>("doubleVariable", 3);

	variableManager.deleteVariable(intVariable);
	variableManager.compact();

	ASSERT_EQ(variableManager.getNbVariables(), 1);
}

#endif     // VARIABLE_MANAGER_TEST_H
