//
// Created by bourmaudp on 30/05/22.
//

#ifndef GMDS_Q_LEARNING_ALGO_H
#define GMDS_Q_LEARNING_ALGO_H

/*----------------------------------------------------------------------------*/
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"
#include "gmds/io/VTKReader.h"

#include <iostream>
#include <fstream>
#include <string>


#include "gmds/ig/Mesh.h"
#include "GMDSIgAlgo_export.h"
#include "gmds/paul/Grid.h"
#include "gmds/paul/Tools.h"
#include "gmds/paul/Actions_Agent.h"
#include "gmds/paul/Environment.h"
#include"gmds/paul/Politique.h"
/*----------------------------------------------------------------------------*/

namespace gmds{

void executeTrainQlearning(Environment environment);


void executeQLearning(Environment environment);

/** \brief verify if is the terminal state (global IoU > 0.9 )
 *
 * @param environment the environment
 * @return false if not the terminal state, else true
 */
bool isTerminalState(Environment environment);

auto getQTable();


void initQTable();

void saveQTable(std::vector<std::vector<double>> theQTable,std::string NameFile);

std::vector<std::vector<double>> readQTable(std::vector<std::vector<double>> theQTable,std::string NameFile);




}

#endif     // GMDS_Q_LEARNING_ALGO_H
