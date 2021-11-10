#!/usr/local/bin/python3

################################################################
# Before Start this Assignment, what I have thought

# 1. Check the inputs is correct (without special characers)
# 2. Check the number of relevent protein before download it all (more than 3000 or equal to 0 will start over, less than 25 will give warning)
# 3. Check species number (too small?)(defalt)
# 4. Check if the file or folder is existed before creating it
# 5. Check the conservation level and plot them correctly
# 6. How to multiple scan the the motif by prosite database
# 7. Conclude the results about motif

################################################################
# Import section
import string
import sys
import os
import subprocess
import re
import shutil

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
    tax = input_deco("Pleas tell me the name of its taxonomic group: " )
    while test_special(tax) != True:
        tax = input_deco("You in put included special character which is not allowed. \nPlease tell me the taxonomic group name: " )
    for character in tax:
        if character not in string.ascii_letters :
            print_deco("I don't think number should be in taxonomic group name, we will search without group name.")
            tax = False
    print_deco("Finding the uids that related to " + str(pro))
    # Use esearch then efetch to get the UID values
    return pro, tax
    

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
    count = pro_count(pro)
    while count >= 1000 :
        pro = input_deco("WARNINGS! There are more than 1000 relative protein! Please check if you have enter the right protein name! \nPlease enter the right name of the protein family:\n")
        ask_the_aim_pro()
    while count == 0 :
        pro = input_deco("WARNINGS! \nThere is no such "+str(pro)+" protein in the database, Please start over!\nPlease enter the right name of the protein family:\n")
        ask_the_aim_pro()
    if count <= 5:
        i = input_deco("WARNINGS! There are less than 5 relative protein! Please check if you have enter the right protein name! \nIf you want to start over please type 'y', if you want to continue with this data please type 'c'")
        if i.lower() != "c":
            print_deco("Sorry. I don't understand what you mean. Let's start over.")
            ask_the_aim_pro()


def test_number(n):
    while True:
        try :
            int(n)
            break
        except :
            n = input("I don't think you entered a number, please enter again. ")
    return int(n)

def seq_num_check(num):
    # If the number of sequences is more than 1000, we can't continue due to the huge data.
    if num >= 1000:
        print_deco("Opps! There are more than 1000 sequences which is too large for us to handle. \nPlease restart this program with other protein.")
        exit()
    

def clustalo_align(pro, tax):   
    print_deco("I'm currently aligning the sequences.")
    # Align all the protein sequences using clustalo
    # With '--force' we don't need to check whether the ouput file is exsited.
    subprocess.call('clustalo -i ' + str(pro) + '_' + str(tax) + '.fasta -o ' + str(pro) + '_' + str(tax) + '_align.fasta -v --force', shell =True)


def plotcon(filename):
    # Plot consercation of a sequence alignment.
    print_deco("Now, let's plot consercation of a sequence alignment.")
    file_check("plotcon.svg")
    subprocess.call('plotcon -sequence ' + str(filename) + ' -winsize 4 -graph svg', shell=True)
    # Ask the user if he/she wants to see the result.
    ans = input_deco("Do you want to see the graph of the conservaton? (y/n)")
    if ans.lower() == 'y':
        # Use firefox to open the graph in order to let the uesr see the graph.
	    subprocess.call('firefox plotcon.svg&', shell=True)
    print_deco("This graph is saved, you could check it later if you want to.")


def clustalo_percent_change(pro, tax):
    # Ask the user if he/she want to calculate the percetage changes by using infoalign.
    ans = input_deco("Do you want to calculate the percetage changes by using infoalign? (y/n) ")
    if ans.lower() == 'y':
        #download name part and percetage change part, and make a dictionary after sort.
        print_deco("I'm calculating now, just be patient.")
        file_check('output1.infoalign')
        file_check('output2.infoalign')
        subprocess.call('infoalign -sequence ' + str(pro) + '_' + str(tax) + '_align.fasta -only -name -outfile output1.infoalign', shell=True)
        subprocess.call('infoalign -sequence ' + str(pro) + '_' + str(tax) + '_align.fasta -only -change -outfile output2.infoalign', shell=True)
        f1 = open("output1.infoalign").read().rstrip("\n") 
        keys = list(f1.split("\n")) 
        f2 = open("output2.infoalign").read().rstrip("\n") 
        list2 = f2.split("\n") 
        values = list(map(float, list2))
        dic = dict(zip(keys, values))
        dic = sorted(dic.items(), key=lambda values:values[1])
        ans = input_deco("The percentage changes are sorted from low to hign.\nDo you want to store them? (y/n) ")
        print(dic)
        if ans == 'y':
            # Store the sequeces to the file called acc_bitscore_sorted.txt
            file_check("clustalo_percent_change.txt")
            file = open("clustalo_percent_change.txt", "w")
            for i in dic:
                file.write("\n" + str(i)) 
            file.close()


def create_blastdb(pro, tax):
    # Make a blast database from fasta files and save them to an output for user when doing blastp processing. 
    print_deco('I am creating a blastdb from fasta files')
    file_check(str(pro) + "_" + str(tax) + "_database")
    subprocess.call("makeblastdb -in " + str(pro) + "_" + str(tax) + ".fasta -dbtype prot -out " + str(pro) + "_" + str(tax) + "_database", shell = True)
    print_deco('I have now created blastdb.')








################################################################
# Actually running part
### Part 1. Find out want the user want, and download the data
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
seq_num = len(open(str(pro) + '_' + str(tax) + '.txt').readlines())
species_num = len(set(open('{0}_{1}.txt'.format(pro,tax)).readlines()))
# Check if the number of sequences is too large.
seq_num_check(seq_num)

# Tell the user what we have found and whether the user want to continue.
ans = input_deco("I have found " + str(species_num) + " species and " + str(pro) + " from what you ask for.\nDo you want to continue analyze? \nPress any key to continue, type 'q' to exit this program.")
if ans == 'q':
    print_deco("Have a nice day!/nBye~")
    exit()



### Part 2. Analyze the conservation between the sequences we found.
# Download the fasta files to the current directory.
print_deco("Sure! Now I will download the fasta files that related to what you searching for.")
subprocess.call('esearch -db protein -query "{0}[prot]" -organism "{1}" | efetch -format fasta > {0}_{1}.fasta'.format(pro,tax), shell=True)

# Using clustalo to align the consevation area beatween all the sequences.
clustalo_align(pro, tax)

# Plot consercation of a sequence alignment.
plotcon(str(pro) + '_' + str(tax) + '_align.fasta')

# Calculate the percentage changes.
clustalo_percent_change(pro, tax)

# Starting to do blastdb
create_blastdb(pro, tax)

# Using the EMBOSS Cons to get the consensus sequence from Clustalo output.
print_deco('Using the EMBOSS Cons to get the consensus sequence from Clustalo output.')
file_check(str(pro) + '_' + str(tax) + '_consensus.fasta')
subprocess.call('cons -plurality 0.8 -sequence ' + str(pro) + '_' + str(tax) + '_align.fasta -outseq ' + str(pro) + '_' + str(tax) + '_consensus.fasta', shell=True) 
print_deco('Consensus has been created now.')

# Using blastp to blast the protein that the user is interested in.
print_deco('Using blastp to search this file: ' + str(pro) + '_' + str(tax) + '_consensus.fasta in blastdb')
subprocess.call('blastp -db ' + str(pro) + '_' + str(tax) + '_database -query ' + str(pro) + '_' + str(tax) + '_consensus.fasta -outfmt 7 -out ' + str(pro) + '_' + str(tax) + '_blast.txt', shell = True)

# Creat a dictionary that contains the acc and the bit score.
blastp_result = {}
line = []
for i in open(str(pro) + '_' + str(tax) + '_blast.txt', 'r'):
    if '#' not in i:
        line.append(i)
print(line)
with open(str(pro) + '_' + str(tax) + '_blast_lines.txt', 'w', encoding="utf-8") as my_file:
    my_file.writelines(line)
for i in line :
    piece = i.split('\t')
    blastp_result[piece[1]] = piece[-1]

# Ask the user to enter a number of most similar sequence.
ans = input_deco("I have sorted the related sequences by their bit score.\nHow many sequence you want to store? \nPlease enter an number.")
while True :
    ans = test_number(ans)
    if ans <= 0:
        ans = input_deco("The number is negtive. \nPlease enter again.")
    else :
        break
print_deco("Now I'm storing the sequence for you.")

# Store the sequeces to the file called acc_bitscore_sorted.txt
file_check("acc_bitscore_sorted.txt")
file = open("acc_bitscore_sorted.txt", "w")
result = {}
for i in sorted(blastp_result.items(),key=lambda item:item[1],reverse=True)[:ans]:
    result[i[0]] = i[1]
file.write("\n".join(result.keys())) 
file.close()

# Use the provided 'pullseq' programme to extract the sequences from the fasta file.
file_check("acc_bitscore_sorted.fasta")
print_deco("Now I'm ecxtracting sequencs from a fasta file")
subprocess.call('/localdisk/data/BPSM/Assignment2/pullseq -i ' + str(pro) + '_' + str(tax) + '_align.fasta -n subject_acc.txt > acc_bitscore_sorted.fasta -v', shell=True)

# Use plotcon to plot the similar sequences alignment result 
ans = input_deco("Do you want to see the plot of blast result by plotcon? (y/n)")
if ans.lower() == 'y':
    plotcon("acc_bitscore_sorted.fasta")










