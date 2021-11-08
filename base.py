#!/usr/local/bin/python3

import sys
import os, subprocess, shutil




# Ask the user to enter the aim protein
def ask_the_aim_pro():
    # Check if there are file with same name, if there are, remove the file
    path = os.getcwd()
    if os.path.isfile("aim.uid"):
        os.remove("aim.uid")

    pro = input("Please tell me the protein name you are interested in: " )
    a = input("Do you know which taxonomic group it is?(y/n) " )
    if a == 'y':
        group = input("Pleas tell me the name of its taxonomic group: " )
    else:
        print("It's ok that you don't know the taxonomic group, we can continue searching")
        group=True

    print("Finding the uids that related to " + str(pro))
    # Use esearch then efetch to get the UID values
    if group:
        subprocess.call("esearch -db protein -query " + str(pro) + " | efetch -format uid > aim.uid",shell=True)
    else:
        subprocess.call("esearch -db protein -query " + str(pro) + " AND " + str(group) + " | efetch -format uid > aim.uid",shell=True)

# Count the number of relative proteins
def count_relative_protein():
    f'Counting the number of relative protein'
    count = len(open(r"aim.uid").readlines())
    if count >= 1000:
        print("WARNINGS! There are more than 1000 relative protein! Please check if you have enter the right protein name!")
        i = input("If you want to choose another protein please type 'y', if you want to continue with this data please type 'c'")
        if i != "c":
            print("Sorry. I don't understand what you mean. Let's start over.")
            ask_the_aim_pro()
    elif count <= 25:
        print("WARNINGS! There are less than 25 relative protein! Please check if you have enter the right protein name!")
        i = input("If you want to start over please type 'y', if you want to continue with this data please type 'c'")
        if i != "c":
            print("Sorry. I don't understand what you mean. Please start over.")
            ask_the_aim_pro()

# Use esummary to bireflly analyze the variety.
def species_size():
    if os.path.isfile("type.txt"):
        os.remove("type.txt") # If there are a file named "type.txt", then delete it
    count = 0
    a = int(input("How many kinds of species do you want to have in your dataset? "))
    while True:
        if a <= 0:
            a = int(input("I don't think you enter the correct number of species, please enter again."))
        else:
            break
    # sore the uids in a list
    uids = list(open("aim.uid").read().rstrip().split())
    # print(uids)
    for i in uids:
        species = []
        count += 1
        print('Now what we were finding the species for uid: ',i)
        subprocess.call('esummary -db protein -id '+str(i)+' | grep "<SubName>" >> type.txt', shell=True)
        open_file = open("type.txt").read().split("\n ") # Count the number of species bt counting the line in file type, then use 'set' to make each species unique.
        species = list(open_file)
#        print(species)
#        print(len(set(species)))
        if len(set(species)) >= a: # When there are more than a different types of species, then stop summarizing
            print("There are "+str(count)+" proteins in your dataset with "+str(len(set(species))) + " species in it, it contains the number of species you wanted or reached the maximun number of species in the total dataset.")
            break
    return uids[:count], set(species)       








ask_the_aim_pro()
count_relative_protein()
data_uids, species = species_size()
