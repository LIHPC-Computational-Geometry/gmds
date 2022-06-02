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


def executeAction(action, blockSet, faceID):
	if action["type"] == "delete":
		blockSet.deleteBlock(faceID)
	elif action["type"] == "edit":
		blockSet.editCorner(faceID, action["v"], action["axis"], action["range"])

def getAllActions():
	actions = []
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
				reward = blockSetBis.getReward(targetMesh);
				if (reward > maxReward):
					maxReward = reward
					bestAction = action
			executeAction(bestAction, blockSet, faceID)
			res.append((bestAction, faceID))
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

	virtualExpert(blockSet, targetMesh, 10)


testVirtualExpert()