//
// Created by bourmaudp on 30/05/22.
//

#ifndef GMDS_Q_LEARNING_ALGO_H
#define GMDS_Q_LEARNING_ALGO_H

/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
#include "GMDSIgAlgo_export.h"
#include "gmds/paul/Grid.h"
#include "gmds/paul/Tools.h"
#include "gmds/paul/Actions_Agent.h"
#include "gmds/paul/Environment.h"
/*----------------------------------------------------------------------------*/

namespace gmds{

void executeAlgoQlearning(Environment environment);

/** \brief verify if is the terminal state (global IoU > 0.9 )
 *
 * @param environment the environment
 * @return false if not the terminal state, else true
 */
bool isTerminalState(Environment environment);

auto getQTable();

int getNextAction();

void initQTable();




}

#endif     // GMDS_Q_LEARNING_ALGO_H
