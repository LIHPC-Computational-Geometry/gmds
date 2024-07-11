#include <iostream>
#include <memory>

#include <ortools/linear_solver/linear_solver.h>

namespace operations_research {

int SolverMIP()
{
	// Déclaration du solveur OR-TOOLS
	// Création du solveur MILP avec le backend SCIP.
	std::unique_ptr<MPSolver> solver(operations_research::MPSolver::CreateSolver("SCIP"));
	if (!solver) {
		LOG(WARNING) << "Le solveur SCIP n'est pas disponible.";
		return -1;
	}

	// Définition des variables
	const double infinity = solver->infinity();
	operations_research::MPVariable *const x = solver->MakeIntVar(0.0, infinity, "x");
	operations_research::MPVariable *const y = solver->MakeIntVar(0.0, infinity, "y");

	// Définition des contraintes
	operations_research::MPConstraint *const constraint0 = solver->MakeRowConstraint(-infinity, 17.5, "constraint0");
	constraint0->SetCoefficient(x, 1);
	constraint0->SetCoefficient(y, 7);

	operations_research::MPConstraint *const constraint1 = solver->MakeRowConstraint(-infinity, 3.5, "constraint1");
	constraint1->SetCoefficient(x, 1);
	constraint1->SetCoefficient(y, 0);

	// Définition de la fonction objectif
	operations_research::MPObjective *const objective = solver->MutableObjective();
	objective->SetCoefficient(x, 1);
	objective->SetCoefficient(y, 10);
	objective->SetMaximization();

	// Résolution du problème
	const operations_research::MPSolver::ResultStatus result_status = solver->Solve();

	// Vérification de la solution optimale
	if (result_status != operations_research::MPSolver::OPTIMAL) {
		std::cerr << "Le problème n'a pas de solution optimale." << std::endl;
		return -1;
	}

	// Affichage de la solution
	std::cout << "Solution optimale trouvée:" << std::endl;
	std::cout << "Valeur de l'objectif = " << objective->Value() << std::endl;
	std::cout << "x = " << x->solution_value() << std::endl;
	std::cout << "y = " << y->solution_value() << std::endl;

	// Affichage des informations avancées
	std::cout << "Problème résolu en " << solver->wall_time() << " millisecondes" << std::endl;
	std::cout << "Nombre d'itérations : " << solver->iterations() << std::endl;
	std::cout << "Nombre de nœuds explorés : " << solver->nodes() << std::endl;

}
}

int main() {
	operations_research::SolverMIP();
	return 0;
}




