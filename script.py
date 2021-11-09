#!/usr/local/bin/python3

################################################################
# Before Start this Assignment

# Things I need to check
# 1. check the inputs is correct (without special characers)
# 2. check the number of relevent protein before download it all (more than 3000 or equal to 0 will start over, less than 25 will give warning)
# 3. check species number (too small?)(defalt)




################################################################
# Import section
import string
import sys
import os, subprocess, shutil
import re

################################################################
#Function defination area.

# Check if there are file with same name, if there are, remove the file
def file_check(file_name):
    path = os.getcwd()
    if os.path.isfile(file_name):
        os.remove(file_name)


# Check if there are special character in the string
def test_special(s):
    regular = re.compile(r'[a-zA-Z0-9_-]')
    if len(s) == len(regular.findall(s)):
        return True
    else :
        return False


# A decoration of printing        
def print_deco(str):
    os.system('clear')# Clear the screen
    print("*********************************************************\n\n\n")
    print(str+"\n\n\n\n*********************************************************")
    

# A decoration of printing 
def input_deco(str):
    os.system('clear')# Clear the screen
    print("*********************************************************\n\n\n")
    a = input(str+"\n\n\n\n*********************************************************\n\n")
    return a


#check if there is protein inputted with organism inputted in the database.
def pro_tax_check(pro,tax):
	num = int(subprocess.check_output("esearch -db protein -query '{0}[prot]' -organism '{1}' | xtract -pattern ENTREZ_DIRECT -element Count".format(pro,tax),shell=True))
	if num == 0:
		return False
	return True


# Ask the user to enter the aim protein
def ask_the_aim_pro():
    file_check("aim.uids")
    pro = input_deco("Please tell me the protein name you are interested in: " )
    while test_special(pro) != True:
        pro = input_deco("You in put included special character which is not allowed. \nPlease tell me the protein name you are interested in: " )
    a = input_deco("Do you know which taxonomic group it is?(y/n) " )
    if a == 'y':
        tax = input_deco("Pleas tell me the name of its taxonomic group: " )
        while test_special(tax) != True:
            tax = input_deco("You in put included special character which is not allowed. \nPlease tell me the taxonomic group name: " )
        for character in tax:
            if character not in string.ascii_letters :
                print_deco("I don't think number should be in taxonomic group name, we will search without group name.")
                tax = False
    else:
        print_deco("It's ok that you don't know the taxonomic group, we can continue searching.")
        tax = False
    print_deco("Finding the uids that related to " + str(pro))
    # Use esearch then efetch to get the UID values
    if tax != False:
        return pro, tax
    else:
        return pro, "No taxonomic"
    

def mkdir(folder_name):
    # Make a folder in current working space.
    dir = os.path.exists(folder_name)
    # Check if there is an existing directory.
    if not dir:
        os.makedirs(folder_name)
    else :
        shutil.rmtree(folder_name)
        os.makedirs(folder_name) 

# Count the number of relative protein in the database
def pro_count(pro):
	count = int(subprocess.check_output("esearch -db protein -query '{0}' | xtract -pattern ENTREZ_DIRECT -element Count".format(pro),shell=True))
	if count == 0:
		return False
	return True 


def pro_tax_count(pro,tax):
	count = int(subprocess.check_output("esearch -db protein -query '{0}[prot]' -organism '{1}' | xtract -pattern ENTREZ_DIRECT -element Count".format(pro,tax),shell=True))
	if count == 0:
		return False
	return True 


# Count the number of relative proteins
def count_relative_protein(pro, tax):
    print_deco('Counting the number of relative protein')
    if tax == "No taxonomic":
        count = pro_count(pro)
    else :
        count =  pro_tax_count(pro,tax)
    while count >= 30000 :
        pro = input_deco("WARNINGS! There are more than 30000 relative protein! Please check if you have enter the right protein name! \nPlease enter the right name of the protein family:\n")
        ask_the_aim_pro()
    while count == 0 :
        pro = input_deco("WARNINGS! \nThere is no such "+str(pro)+" protein in the database, Please start over!\nPlease enter the right name of the protein family:\n")
        ask_the_aim_pro()
    if count <= 25:
        i = input_deco("WARNINGS! There are less than 25 relative protein! Please check if you have enter the right protein name! \nIf you want to start over please type 'y', if you want to continue with this data please type 'c'")
        if i != "c":
            print_deco("Sorry. I don't understand what you mean. Let's start over.")
            ask_the_aim_pro()


def test_number(n):
    while True:
        try :
            int(n)
            break
        except :
            n = input("I don't think you enter the correct number, please enter again. ")
    return int(n)

def seq_num_check(num):
    # If the number of sequences is more than 2500, we can't continue due to the huge data.
    if num >= 2500:
        print_deco("Opps! There are more than 2500 sequences which is too large for us to handle. \nPlease restart this program with other protein.")
        exit()
    

def clustao_align():   
    file_check("align_output.txt")
    print_deco("I'm currently aligning the sequences.")
    subprocess.call("clustalo -i data.fa -o align_output.txt -v --force", shell=True) # Using clustao to align the conservation interval between the sequences.


def plotcon():
    print_deco("Now, let's plot consercation of a sequence alignment.")
    file_check("plotcon.svg")
    subprocess.call("plotcon -sequences align_output.txt -winsize 4 -graph svg", shell = True)




################################################################
# Actually run the program
### Step 1. Find out want the user want, and download the data
# Ask what the user want and check if the input is correct.
pro, tax = ask_the_aim_pro()
count_relative_protein(pro, tax)

# Make a directory to start analyzing.
print_deco("Thanks for your cooperation, now I made a directory for you named " + str(pro) + str(tax))
mkdir(str(pro) + "_" + str(tax))
os.chdir('{0}/{1}_{2}'.format(os.getcwd(),pro,tax))

# Download the data and count the number of sequences and species. 
file_check(str(pro) + "_" + str(tax) + ".txt")
subprocess.call('esearch -db protein -query "{0}[prot]" -organism "{1}" | efetch -format docsum | xtract -pattern DocumentSummary -element Organism > {0}_{1}.txt'.format(pro,tax), shell=True)

# Count the number of sequences and species.
seq_num = len(open('{0}_{1}.txt'.format(pro,tax)).readlines())
species_num = len(set(open('{0}_{1}.txt'.format(pro,tax)).readlines()))
# Check if the number of sequences is too large.
seq_num_check(seq_num)

# Tell the user what we have found and whether the user want to continue.
ans = input_deco("I have found " + str(species_num) + "species and " + str(pro) + "from what you ask for.\nDo you want to continue analyze? \nPress any key to continue, type 'q' to exit this program.")
if ans == 'q':
    print_deco("Have a nice day!/nBye~")
    exit()


