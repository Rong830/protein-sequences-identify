#!/usr/local/bin/python3

import sys
import os, subprocess, shutil


# Check if there are file with same name, if there are, remove the file
def file_check(file_name):
    path = os.getcwd()
    if os.path.isfile(file_name):
        os.remove(file_name)


# Ask the user to enter the aim protein
def ask_the_aim_pro():
    file_check("aim.uids")
    pro = input("Please tell me the protein name you are interested in: " )
    a = input("Do you know which taxonomic group it is?(y/n) " )
    if a == 'y':
        group = input("Pleas tell me the name of its taxonomic group: " )
        for character in group:
            if character not in string.ascii_letters :
                print("\nI don't think number should be in taxonomic group name, we will search without group name.")
                group=True
    else:
        print("It's ok that you don't know the taxonomic group, we can continue searching.")
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
    count = 0
    a = input("How many kinds of species do you want to have in your dataset? ")
    while type(a) != int or a <= 0:
        a = input("I don't think you enter the correct number of species, please enter again.")
    # Store the uids into a list
    uids = list(open("aim.uid").read().rstrip().split())
#    print(uids)
    file_check("type.txt")
    file_check("all_summary.txt")
    for i in uids:
        species = []
        count += 1
        print('Now what we were finding the species for uid: ',i)
        subprocess.call('esummary -db protein -id '+str(i)+' >> all_summary.txt', shell=True)
        subprocess.call('cat all_summary.txt | grep "<Organism>" >> type.txt', shell=True)
        open_file = open("type.txt").read().split("\n ") # Count the number of species bt counting the line in file type, then use 'set' to make each species unique.
        species = list(open_file)
#        print(species)
#        print(len(set(species)))
        if len(set(species)) >= a: # When there are more than a different types of species, then stop summarizing
            print("There are "+str(count)+" proteins in your dataset with "+str(len(set(species))) + " species in it, it contains the number of species you wanted or reached the maximun number of species in the total dataset.")
            break
    file_check("species.txt")
    # Write the total species names in a file called "species.txt"
    file_out = open("species.txt", "w")
    for i in set(species):   
        file_out.write(str(i) + "\n")
    file_out.close()
    return uids[:count], set(species)    


def efetch_fasta_file(data_uids):
    # Use edirect to get the sequence as fasta format, and save it to a file
    file_check("data.fa")
    for i in data_uids:
        print("I'm downloading the fasta file about " + str(i))
        subprocess.call("esearch -db protein -query "+i+" | efetch -format fasta >> data.fa", shell=True)


def clustao_align():   
    file_check("align_output.txt")
    print("I'm currently aligning the sequences.")
    subprocess.call("clustalo -i data.fa -o align_output.txt -v") # Using clustao to align the conservation interval between the sequences.


ask_the_aim_pro()
os.system('clear')# Clear the screen

count_relative_protein()
os.system('clear')# Clear the screen

data_uids, species = species_size()
os.system('clear')# Clear the screen

efetch_fasta_file(data_uids)
os.system('clear')# Clear the screen

clustao_align()
os.system('clear')# Clear the screen


