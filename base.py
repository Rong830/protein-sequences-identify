#!/usr/local/bin/python3

import sys
import os, subprocess, shutil


# Check if there are file with same name, if there are, remove the file
path = os.getcwd()
if os.path.isfile("aim.uid"):
    os.remove("aim.uid")

# Make the bash script excutable
subprocess.call("chmod 700 esearch.sh",shell=True)
# Ask the user to enter the aim protein
subprocess.call("./esearch.sh ",shell=True)

# Count the number of relative proteins
f'Counting the number of relative protein'
count = len(open(r"aim.uid").readlines())
if count >= 1000:
    print("WARNINGS! There are more than 1000 relative protein! Please check if you have enter the right protein name!")
    print("If you want to start over please type 'y', if you want to continue with this data please type 'c'")
elif count <= 25:
    print("WARNINGS! There are less than 25 relative protein! Please check if you have enter the right protein name!")





