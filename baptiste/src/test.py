import os
from pathlib import Path

# set working directory
directory = str(Path(__file__).resolve().parents[2])
directory += "/cmake-build-debug/baptiste"
os.chdir(directory)
print("Current working directory: {0}".format(os.getcwd()))



from environment import RLBlockSet

blockSet = RLBlockSet()

ar = blockSet.LinearSpacedArray(0, 5, 6)

print(type(ar))
print(ar)
def main():
	double xMin = 0
	double yMin = 0
	double xMax = 9
	double yMax = 9
	blockSet.setFrame(xMin, yMin, xMax, yMax)

	nbFaces = blockSet.countBlocks()
	print("Number of blocks : " + nbFaces + "\n")

	for faceID in blockSet.m_mesh.faces():
		print(faceID + "\n")
		if faceID > 4:
			blockSet.deleteBlock(faceID)
			print("Face " + faceID + " deleted" + "\n")

	nbFaces = blockSet.countBlocks();
	print("Number of blocks : " + nbFaces + "\n")

	for faceID in blockSet.m_mesh.faces():
		blockSet.editCorner(faceID, False, "y", -3)

	blockSet.saveMesh("MyBlockSetPython");