import os
from pathlib import Path

# set working directory
directory = str(Path(__file__).resolve().parents[2])
directory += "/cmake-build-debug/baptiste"
os.chdir(directory)
print("Current working directory: {0}".format(os.getcwd()))



from environment import RLBlockSet, getVolFrac, Tools



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

	nbFaces = blockSet.countBlocks();
	print("Number of blocks : " + str(nbFaces) + "\n")


	for faceID in blockSet.getAllFaces():
		blockSet.editCorner(faceID, False, "y", -3)


	print("Corner editing done")
	blockSet.saveMesh("MyBlockSetPython2")
	print("Figure saved")
	
	print("End of function")

#main()

def testReadMesh():
	m = readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk")
	print(m.getNbNodes())

#testReadMesh()

def testVolFrac():
	t = Tools()
	t.readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk")
	AImprintMesh = t.m_mesh
	print(AImprintMesh.getNbNodes())
	print("Read Ok")
	blockSet = RLBlockSet(2)
	print("Dim2 Ok")
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 8, 4)
	print("Set from file ok")
	AMesh = blockSet.m_mesh
	print("The files have been read successfully" + "\n")

	"""
	for faceID in blockSet.getAllFaces():
		print(str(faceID) + "\n")
		if faceID > 25:
			blockSet.deleteBlock(faceID);
			print("Face " + str(faceID) + " deleted" + "\n")
	
	"""

	#t.computeVolFrac(AMesh, AImprintMesh)
	#blockSet.getReward(AImprintMesh)
	print("Computation done")

	"""
	#res = getVolFrac(AMesh)

	#print(res)
	#print(type(res))

	for faceID in blockSet.getAllFaces():
		blockSet.editCorner(faceID, false, "y", -3)
		blockSet.editCorner(faceID, true, "x", 2)
	"""

	#volfraccomputation_2d(&AMesh, &AImprintMesh, AMesh.getVariable<double, GMDS_FACE>("volFrac"));

	#blockSet.saveMesh("MyBlockSet")
	print("End of function")


testVolFrac()
