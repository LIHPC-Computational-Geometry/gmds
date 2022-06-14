import os
from pathlib import Path

# set working directory
directory = str(Path(__file__).resolve().parents[2])
directory += "/cmake-build-debug/baptiste"
os.chdir(directory)
print("Current working directory: {0}".format(os.getcwd()))



from environment import RLBlockSet, Capsule, cloneBlockSet



def main():
	blockSet = RLBlockSet()
	xMin = 0.
	yMin = 0.
	xMax = 9.
	yMax = 9.
	blockSet.setFrame(xMin, yMin, xMax, yMax, 3, 3)

	nbFaces = blockSet.countBlocks()
	print("Number of blocks : " + str(nbFaces) + "\n")

	print(type(blockSet.getAllFaces()))
	for faceID in blockSet.getAllFaces():
		print(str(faceID) + "\n")
		if faceID > 4:
			blockSet.deleteBlock(faceID)
			print("Face " + str(faceID) + " deleted" + "\n")

	nbFaces = blockSet.countBlocks()
	print("Number of blocks : " + str(nbFaces) + "\n")


	for faceID in blockSet.getAllFaces():
		blockSet.editCorner(faceID, False, "y", -3)


	print("Corner editing done")
	blockSet.saveMesh("MyBlockSetPython2")
	print("Figure saved")
	
	print("End of function")

#main()


def testReward():
	c = Capsule()
	c.readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk")
	targetMesh = c.m_mesh
	blockSet = RLBlockSet()
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 3, 3)
	print("Set from file ok")
	AMesh = blockSet.m_mesh
	print("The files have been read successfully" + "\n")
	x = blockSet.getReward(targetMesh)
	print(x)

#testReward()

def executeAction(action, blockSet, faceID):
	if action["type"] == "delete":
		blockSet.deleteBlock(faceID)
	elif action["type"] == "edit":
		blockSet.editCorner(faceID, action["v"], action["axis"], action["range"])

def getAllActions():
	actions = []
	actions.append({"type" : "stay"})
	actions.append({"type" : "delete"})
	for v in [False, True]:
		for axis in ["x", "y"]:
			for arange in [-2, -1, 1, 2]:
				action = {"type" : "edit",
						  "v" : v,
						  "axis" : axis,
						  "range" : arange}
				actions.append(action)
	return actions

def getAllActionsBis(blockSet):
	id = 0
	actions = []
	actions.append({"id" : id, "type" : "stay", "face" : None})
	id += 1
	for faceID in blockSet.getAllFaces():
		actions.append({"id" : id, "type" : "delete", "face" : faceID})
		for v in [False, True]:
			for axis in ["x", "y"]:
				for arange in [-2, -1, 1, 2]:
					id += 1
					action = {"id" : id,
							  "type" : "edit",
							  "v" : v,
							  "axis" : axis,
							  "range" : arange,
							  "face" : faceID}
					actions.append(action)
	return actions


def virtualExpert(blockSet, targetMesh, nMax):
	res = []
	actions = getAllActions()
	for step in range(nMax):
		print("Step n°" + str(step))
		for faceID in blockSet.getAllFaces():
			print("Face id = " + str(faceID))
			maxReward = 0
			bestAction = actions[0]
			for action in actions:
				blockSetBis = RLBlockSet()
				cloneBlockSet(blockSet, blockSetBis)
				if step < nMax/2 and action["type"] == "edit":
					executeAction(action, blockSetBis, faceID)
				else:
					executeAction(action, blockSetBis, faceID)
				if blockSetBis.isValid():
					reward = blockSetBis.getReward(targetMesh);
					if (reward > maxReward):
						maxReward = reward
						bestAction = action
				else:
					print("Bad cell detected")
			if bestAction["type"] == "delete":
				print("A delete action has been selected")
			executeAction(bestAction, blockSet, faceID)
			res.append((bestAction, faceID))
	return res


def executeActionBis(action, blockSet):
	if action["type"] == "stay":
		return True
	elif action["face"] in blockSet.getAllFaces():
		if action["type"] == "delete":
			blockSet.deleteBlock(action["face"])
		elif action["type"] == "edit":
			blockSet.editCorner(action["face"], action["v"], action["axis"], action["range"])
		return True
	else:
		return False

def virtualExpertBis(blockSet, targetMesh, nMax):
	res = []
	actions = getAllActionsBis(blockSet)
	#for step in range(nMax):
	stop = False
	step = 0
	while not stop:
		step += 1
		print("Step n°" + str(step))
		maxReward = 0
		bestAction = actions[0]
		for action in actions:
			print(action)
			blockSetBis = RLBlockSet()
			cloneBlockSet(blockSet, blockSetBis)
			if step < nMax/2 and action["type"] == "edit":
				executeActionBis(action, blockSetBis)
			else:
				executeActionBis(action, blockSetBis)
			if blockSetBis.isValid():
				reward = blockSetBis.getReward(targetMesh)
				if (reward > maxReward):
					maxReward = reward
					bestAction = action
			else:
				print("Action can't be executedBad cell detected")
		if bestAction["type"] == "delete":
			print("A delete action has been selected")
			actions = list(filter(lambda a: a["face"] != bestAction["face"], actions))
		if bestAction["type"] == "stay":
			stop = True
		executeActionBis(bestAction, blockSet)
		res.append(bestAction)
	return res



def testDeepCopy():
	blockSet = RLBlockSet()

	xMin, yMin, xMax, yMax = 0., 0., 9., 9.
	blockSet.setFrame(xMin, yMin, xMax, yMax, 3, 3)
	
	print("Number of blocks : " + str(blockSet.countBlocks()) + "\n")

	blockSet2 = RLBlockSet(2)
	cloneBlockSet(blockSet, blockSet2)

	print("Number of blocks 2 : " + str(blockSet2.countBlocks()) + "\n")

	blockSet2.deleteBlock(0)

	print("Number of blocks : " + str(blockSet.countBlocks()) + "\n")
	print("Number of blocks 2 : " + str(blockSet2.countBlocks()) + "\n")


def testVirtualExpert(n):
	c = Capsule()
	c.readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk")
	#c.readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/HolesInSquare0.vtk")
	targetMesh = c.m_mesh
	print(type(targetMesh))
	blockSet = RLBlockSet()
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 3, 3)
	#blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/HolesInSquare0.vtk", 3, 3)
	res = virtualExpertBis(blockSet, targetMesh, n)
	print(res)
	print(len(res))
	blockSet.saveMesh("/home/bonyb/Documents/virtualExpertPyCurvedShape.vtk")
	exit(0)


testVirtualExpert(2)




# TO DO : Q-LEARNING



#import numpy as np
import random
#from IPython.display import clear_output


alpha = 0.1
gamma = 0.6
epsilon = 0.5



def getTargetShape(filename):
	c = Capsule()
	c.readMesh(filename)
	targetShape = c.m_mesh
	return targetShape



def step(state, action, targetMesh):
	nextState = RLBlockSet()
	cloneBlockSet(state, nextState)
	possible = executeActionBis(action, nextState)
	if not possible:
		reward = -100
	else:
		if nextState.isValid():
			reward = nextState.getReward(targetMesh)
		else:
			reward = -100
	print("Reward : " + str(reward))
	terminated = (reward > 0.95)
	if len(state.getAllFaces()) == 0:
		terminated = True
	return nextState, reward, terminated

def generateID(text):
  hash = 0
  for ch in text:
    hash = ( hash*281  ^ ord(ch)*997) & 0xFFFFFFFF
  return hash


def readTable(qTable, state, action):
	value = 0
	if (generateID(state.getStateID()), action["id"]) in qTable:
		value = qTable[generateID(state.getStateID()), action["id"]]
	return value

def updateTable(qTable, state, action, value):
	qTable[generateID(state.getStateID()), action["id"]] = value

import math

def train(filename, n = 1):
	qTable = dict()
	c = Capsule()
	c.readMesh(filename)
	targetShape = c.m_mesh
	print("Number of nodes in target mesh : " + str(targetShape.getNbNodes()))
	blockSet = RLBlockSet()
	blockSet.setFromFile(filename, 3, 3)
	actions = getAllActionsBis(blockSet)
	for episode in range(n):
		state = RLBlockSet()
		cloneBlockSet(blockSet, state)

		reward = 0
		terminated = False
		k = 0

		while not terminated:
			print("Iteration n°" + str(k))
			print(qTable)
			if random.uniform(0, 1) < epsilon:
				action = random.choice(actions)
			else:
				try:
					actionID = max(qTable, key = qTable.get)[1]
					action = list(filter(lambda a: a["id"] == actionID, actions))[0]
				except:
					action = random.choice(actions)

			nextState, reward, terminated = step(state, action, targetShape)
			oldValue = readTable(qTable, state, action)
			interTable = { key: value for key, value in qTable.items() if key[0] == nextState }
			try:
				nextMax = qTable[max(interTable, key = interTable.get)]
			except:
				nextMax = 0
			newValue = (1 - alpha) * oldValue + alpha * (reward + gamma * nextMax)
			"""
			if math.isnan(newValue):
				print("reward : " + str(reward))
				print("action : " + str(action))
				print("nextMax : " + str(nextMax))
				print("oldValue : " + str(oldValue))
				exit(5)
			"""

			updateTable(qTable, state, action, newValue)
			
			state = nextState
			
			k += 1
		
		if (episode + 1) % 100 == 0:
			clear_output(wait=True)
			print("Episode: {}".format(episode + 1))
			enviroment.render()

	print("**********************************")
	print("Training is done!\n")
	print("**********************************")
	state.saveMesh("/home/bonyb/Documents/resultOfQLearning.vtk")


#train("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 10)



def testSerialize():
	c = Capsule()
	c.readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk")
	targetMesh = c.m_mesh
	#print(type(targetMesh))
	blockSet = RLBlockSet()
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 3, 3)

	blockSetBis = RLBlockSet()
	cloneBlockSet(blockSet, blockSetBis)


	x = blockSet.getStateID()
	print(x)

	id_ = generateID(x)
	print(id_)
	print(generateID(blockSetBis.getStateID()))


def testIoU():
	c = Capsule()
	
	c.readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk")
	targetMesh = c.m_mesh
	print(type(targetMesh))
	blockSet = RLBlockSet()
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 3, 3)
	print("Alles in Ordnung")
	blockSet.getReward(targetMesh)
	"""
	blockSet.saveMesh("/home/bonyb/Documents/virtualExpertPy.vtk")
	"""
	#c.saveMesh("/home/bonyb/Documents/targetMeshTestIoU.vtk")
	exit(0)

#testIoU()

def testOverlay():
	c = Capsule()
	c.readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk")
	targetMesh = c.m_mesh
	blockSet = RLBlockSet()
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 3, 3)
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())
	blockSet.editCorner(0, False, "x", 2)
	print(blockSet.overlay())


#testOverlay()


def metropolis(blockSet, targetMesh):
	reward = blockSet.getReward(targetMesh)
	actions = getAllActionsBis(blockSet)

	while reward < 0.95:
		print(reward)

		action = random.choice(actions)

		blockSetBis = RLBlockSet()
		cloneBlockSet(blockSet, blockSetBis)

		executeActionBis(action, blockSetBis)

		if blockSetBis.isValid():
			newReward = blockSetBis.getReward(targetMesh)

			ratio = newReward/reward

			if ratio > 1:
				executeActionBis(action, blockSet)
				reward = newReward

			else:
				eps = random.uniform(0.9, 1)

				if ratio > eps:
					executeActionBis(action, blockSet)
					reward = newReward
				else:
					pass
		else:
			pass



def testMetropolis():
	c = Capsule()
	c.readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk")
	targetMesh = c.m_mesh
	blockSet = RLBlockSet()
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 3, 3)
	metropolis(blockSet, targetMesh)
	blockSet.saveMesh("/home/bonyb/Documents/metropolis.vtk")
	exit(0)

#testMetropolis()