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