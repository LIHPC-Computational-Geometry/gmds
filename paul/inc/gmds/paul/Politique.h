//
// Created by bourmaudp on 31/05/22.
//

#ifndef GMDS_POLITIQUE_H
#define GMDS_POLITIQUE_H

/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
#include "GMDSIgAlgo_export.h"
#include "gmds/paul/Grid.h"
#include "gmds/paul/Tools.h"
#include "gmds/paul/Actions_Agent.h"
#include "gmds/paul/Environment.h"
/*----------------------------------------------------------------------------*/
#include <time.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
class Politique
{
 public:
	Politique(Environment *AEnv);

	void initQTable();

	int getNextAction(int intervalIoU);

	int getNextActionQLearning(int intervalIoU);

	void updateQTable(int intervalIoU, int actionSelect, double newQValue);

	int getInterval(double localIoU);

	double maxQValue(int intervalIoU);


	std::vector<std::vector<double>> getQTable();

	std::vector<std::vector<double>> Q_Table;

	gmds::Environment env;
};

}

#endif     // GMDS_POLITIQUE_H
