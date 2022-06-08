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
	actions = []
	actions.append({"type" : "stay"})
	for faceID in blockSet.getAllFaces():
		actions.append({"type" : "delete", "face" : faceID})
		for v in [False, True]:
			for axis in ["x", "y"]:
				for arange in [-2, -1, 1, 2]:
					action = {"type" : "edit",
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
	if action["type"] == "delete":
		blockSet.deleteBlock(action["face"])
	elif action["type"] == "edit":
		blockSet.editCorner(action["face"], action["v"], action["axis"], action["range"])

def virtualExpertBis(blockSet, targetMesh, nMax):
	res = []
	actions = getAllActionsBis(blockSet)
	for step in range(nMax):
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


def testVirtualExpert():
	c = Capsule()
	c.readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk")
	targetMesh = c.m_mesh
	blockSet = RLBlockSet()
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 3, 3)
	res = virtualExpertBis(blockSet, targetMesh, 10)
	print(res)
	print(len(res))
	blockSet.saveMesh("/home/bonyb/Documents/virtualExpertPy.vtk")
	exit(0)



#testVirtualExpert()




# TO DO : Q-LEARNING



#import numpy as np
import random
#from IPython.display import clear_output


alpha = 0.1
gamma = 0.6
epsilon = 0.1
q_table = dict()



def getTargetShape(filename):
	c = Capsule()
	c.readMesh(filename)
	targetShape = c.m_mesh
	return targetShape


def training(m, a, alpha = 0.1, gamma = 0.6, epsilon = 0.1, n = 5):
    actions = ['E', 'W', 'N', 'S']
    q_table = create_table(m, actions)
    for i in range(n):
        a.position = (m.rows, m.cols)
        # a.position = (random.randint(1, m.rows), random.randint(1, m.cols))
        state = a.position
        epochs, penalties, reward, = 0, 0, 0
        done = False
        while not done:
            if random.uniform(0, 1) < epsilon:
                action = random.choice(actions)
            else:
                action = max(q_table, key = q_table.get)[1]
                
            next_state, reward, done = step(m, a, action)
            old_value = q_table[state, action]
            
            # next_max = np.max(q_table[next_state])
            inter_table = { key: value for key, value in q_table.items() if key[0] == next_state }
            next_max = q_table[max(inter_table, key = inter_table.get)]

            new_value = (1 - alpha) * old_value + alpha * (reward + gamma * next_max)
            q_table[state, action] = new_value
            
            state = next_state
            epochs += 1
            
            clear_output(wait=True)
            print(f"Episode: {i}")
            print(f"Epoch: {epochs}")
            # print(q_table)
    return q_table


def step(state, action, targetMesh):
	nextState = RLBlockSet()
	cloneBlockSet(state, nextState)
	executeActionBis(action, nextState)
	if nextState.isValid():
		reward = nextState.getReward(nextState, targetMesh)
	else:
		reward = -100
	terminated = reward > 0.95
	return nextState, reward, terminated



def train(filename, n = 10):
	targetShape = getTargetShape(filename)
	blockSet = RLBlockSet()
	blockSet.setFromFile(filename, 3, 3)
	actions = getAllActionsBis(blockSet)
	for episode in range(n):
		# Reset the enviroment
		state = RLBlockSet()
		cloneBlockSet(blockSet, state)

		# Initialize variables
		reward = 0
		terminated = False
		epochs = 0

		while not terminated:
			print(epochs)
			# Take learned path or explore new actions based on the epsilon
			if random.uniform(0, 1) < epsilon:
				action = random.choice(actions)
			else:
				# action = np.argmax(q_table[state])
				action = max(q_table, key = q_table.get)[1]

			nextState, reward, terminated = step(state, action, targetShape)
			old_value = q_table[state, action]

			inter_table = { key: value for key, value in q_table.items() if key[0] == nextState }
			nextMax = q_table[max(inter_table, key = inter_table.get)]

			new_value = (1 - alpha) * old_value + alpha * (reward + gamma * nextMax)
			q_table[state, action] = new_value
			
			state = nextState
			epochs += 1
		
		if (episode + 1) % 100 == 0:
			clear_output(wait=True)
			print("Episode: {}".format(episode + 1))
			enviroment.render()

	print("**********************************")
	print("Training is done!\n")
	print("**********************************")


train("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk")