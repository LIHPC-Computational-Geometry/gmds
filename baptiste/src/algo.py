from environment import RLBlockSet, Capsule, cloneBlockSet
import random
import numpy as np

alpha = 0.1
gamma = 0.6
epsilon = 0.4


def executeAction(action, blockSet, faceID):
	print("Trying to execute action ...")
	if action["type"] == "delete":
		print("Trying to delete face ...")
		blockSet.deleteBlock(faceID)
	elif action["type"] == "edit":
		print("Trying to edit face ...")
		blockSet.editCorner(faceID, action["v"], action["axis"], action["range"])



def getAllActions():
	actions = []
	actions.append({"id" : 0, "type" : "stay"})
	actions.append({"id" : 1, "type" : "delete"})
	id = 2
	for v in [False, True]:
		for axis in ["x", "y"]:
			for arange in [-2, -1, 1, 2]:
				action = {"id" : id,
								  "type" : "edit",
						  	  "v" : v,
						  	  "axis" : axis,
						  	  "range" : arange}
				actions.append(action)
				id += 1
	return actions



def step(state, action, targetMesh, faceID):
	nextState = RLBlockSet()
	cloneBlockSet(state, nextState)
	if faceID in nextState.getAllFaces():
		executeAction(action, nextState, faceID)
	else:
		print("Could not execute action on deleted face")
	if nextState.isValid() and nextState.countBlocks() > 1:
		reward = nextState.getReward(targetMesh)
	else:
		reward = -100
	print("Reward : " + str(reward))
	terminated = (reward > 1.1)
	if state.countBlocks() == 0:
		terminated = True
	return nextState, reward, terminated
	
def getCategory(state, faceID, n):
	local = state.getLocalIou(faceID)
	for i, x in enumerate(list(np.linspace(0, 1, n))):
		if local < x:
			return i-1
	return n


def training(filename, n = 1):
	print("Training !!!")
	actions = getAllActions()

	c = Capsule()
	c.readMesh(filename)
	targetShape = c.m_mesh
	blockSet = RLBlockSet()
	blockSet.setFromFile(filename, 3, 2)

	qTables = dict()	
	for faceID in blockSet.getAllFaces():
		qTables[faceID] = np.zeros([6, len(actions)])
	
	

	for episode in range(n):
		print("Episode " + str(episode))

		print(qTables)

		state = RLBlockSet()
		cloneBlockSet(blockSet, state)

		reward = 0
		terminated = False
		k = 0
		reset = False
		while not terminated:
			if reset:
				break
			print("**********************************")
			print("Episode n° " + str(k))
			print("**********************************")
			for faceID in blockSet.getAllFaces():
				qTable = qTables[faceID]
				if random.uniform(0, 1) < epsilon:
					action = random.choice(actions)
				else:
					try:
						actionID = np.argmax(qTable[getCategory(state, faceID, 5)])
						action = list(filter(lambda a: a["id"] == actionID, actions))[0]
					except:
						action = random.choice(actions)

				if state.isValid():
					state.getReward(targetShape)
				else:
					print("Could not compute reward on this state because of bad cell")
					reset = True
					break

				nextState, reward, terminated = step(state, action, targetShape, faceID)
				
				oldValue = qTable[getCategory(state, faceID, 5), action["id"]]

				nextMax = np.max(qTable[getCategory(state, faceID, 5)])

				newValue = (1 - alpha) * oldValue + alpha * (reward + gamma * nextMax)

				qTable[getCategory(state, faceID, 5), action["id"]] = newValue


				state = nextState

				qTables[faceID] = qTable
			
				k += 1
				print("LALALAALLAA")
		print("One episode just finished")
		state.saveMesh("/home/bonyb/Images/figures/newqlearning/train" + str(episode) + ".vtk")
	print("**********************************")
	print("Training is done!\n")
	print("**********************************")
	return qTables


def testing(filename):
	actions = getAllActions()
	c = Capsule()
	c.readMesh(filename)
	targetShape = c.m_mesh
	blockSet = RLBlockSet()
	blockSet.setFromFile(filename, 3, 2)
	for i in range(30):
		for faceID in blockSet.getAllFaces():
			qTable = qTables[faceID]
			if blockSet.isValid():
				blockSet.getReward(targetShape)
			actionID = np.argmax(qTable[getCategory(blockSet, faceID, 5)])
			action = list(filter(lambda a: a["id"] == actionID, actions))[0]
			executeAction(action, blockSet, faceID)
		blockSet.saveMesh("/home/bonyb/Images/figures/newqlearning/result" + str(i) + ".vtk")



qTables = training("/home/bonyb/Documents/GitHub/gmds/test_samples/HolesInSquare0.vtk", 500)

print(qTables)

testing("/home/bonyb/Documents/GitHub/gmds/test_samples/HolesInSquare0.vtk")

print("Number of episodes whihc converge : " + str(counter))