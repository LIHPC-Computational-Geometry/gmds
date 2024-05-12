#include <iostream>
#include <memory>

#include <ortools/linear_solver/linear_solver.h>

int main() {
	std::cout<<"hello world ortools"<<std::endl;

// Create the mip solver with the SCIP backend.
	std::unique_ptr<operations_research::MPSolver> solver(operations_research::MPSolver::CreateSolver("SCIP"));
	if (!solver) {
		LOG(WARNING) << "SCIP solver unavailable.";
		return -1;
	}

	const double infinity = solver->infinity();
// x and y are integer non-negative variables.
	operations_research::MPVariable* const x = solver->MakeIntVar(0.0, infinity, "x");
	operations_research::MPVariable* const y = solver->MakeIntVar(0.0, infinity, "y");

	LOG(INFO) << "Number of variables = " << solver->NumVariables();

	// x + 7 * y <= 17.5.
	operations_research::MPConstraint* const c0 = solver->MakeRowConstraint(-infinity, 17.5, "c0");
	c0->SetCoefficient(x, 1);
	c0->SetCoefficient(y, 7);

// x <= 3.5.
	operations_research::MPConstraint* const c1 = solver->MakeRowConstraint(-infinity, 3.5, "c1");
	c1->SetCoefficient(x, 1);
	c1->SetCoefficient(y, 0);

	LOG(INFO) << "Number of constraints = " << solver->NumConstraints();

	// Maximize x + 10 * y.
	operations_research::MPObjective* const objective = solver->MutableObjective();
	objective->SetCoefficient(x, 1);
	objective->SetCoefficient(y, 10);
	objective->SetMaximization();

	const operations_research::MPSolver::ResultStatus result_status = solver->Solve();
// Check that the problem has an optimal solution.
	if (result_status != operations_research::MPSolver::OPTIMAL) {
		LOG(FATAL) << "The problem does not have an optimal solution!";
	}

	LOG(INFO) << "Solution:";
	LOG(INFO) << "Objective value = " << objective->Value();
	LOG(INFO) << "x = " << x->solution_value();
	LOG(INFO) << "y = " << y->solution_value();
}