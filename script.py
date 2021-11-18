#!/usr/local/bin/python3

################################################################
# Before Start this Assignment, what I have thought

# 1. Check the inputs is correct (without special characers)
# 2. Check the number of relevent protein before download it all (more than 1000 or equal to 0 will start over, less than 5 will give warning)
# 3. Check species number (too small?)(defalt)
# 4. Check if the file or folder is existed before creating it
# 5. Check the conservation level and plot them correctly
# 6. How to multiple scan the the motif by prosite database
# 7. Conclude the results about motif and send summary email to the user
# 8. Check if the user have the imported package and edirect and PROSITE database

# EMBOSS softwares I used
# plotcon
# infoalign
# Cons
# seqretsplit
# patmatmotif
# pepstats
# pepwindowall
# distmat

################################################################
# Import section

import sys
import os
import subprocess
import re
import shutil
import glob
import string

subprocess.call("pip3 install pandas", shell = True)
subprocess.call("pip3 install matplotlib", shell = True)
subprocess.call("pip3 install sklearn", shell = True)
subprocess.call("pip3 install seaborn", shell = True)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
from pandas.api.types import CategoricalDtype
import seaborn as sns
# subprocess.call("pip3 install " + str(package), shell = True)
################################################################
# Function defination area.

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
    input(str+"\n\n\n\n***************(Press any key to continue)***************\n")
    

# A decoration of printing 
def input_deco(str):
    os.system('clear')# Clear the screen
    print("*********************************************************\n\n\n")
    a = input(str+"\n\n\n\n*********************************************************\n\n")
    return a


# Check if the user have edirect in his/her home directory.
def edirect_check():
    print_deco('Checking if you have edirect in you home directory')
    if os.path.isdir(os.path.expanduser('~/edirect')):
        print_deco('You have edirect in you home space, we can continue now.')
    else: 
        print_deco('You do not have installed edirect in your home directory, I will install it for you.')
        subprocess.call('sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"', shell = True)
        print_deco('You have edirect in you home space, we can continue now.')


  # Ask the user where to put the output files, the defual setting is the user home directory.
def output_path():
    ans = 'A'
    while ans.lower() != 'y' and ans.lower() != 'n':
        ans = input_deco("Do you want to store all the output file in your current working directory? (y/n)")
    if ans.lower() == 'y':
        print_deco("Sure! The output files will be stored in your current working directory.")
        return os. getcwd()
    if ans.lower() == 'n':
        ans = input_deco("Sure! Please tell me the path name of where you want to store your output files.")
        if not os.path.exists(ans):
            os.mkdir(ans)
            print_deco('Now we will store everything in '+ str(ans))
        else:
            print('Now we will store everything in '+ str(ans))
        return ans 


#check if there is protein inputted with organism inputted in the database.
def pro_tax_check(pro,tax):
	num = int(subprocess.check_output("esearch -db protein -query '{0}[prot]' -organism '{1}' | xtract -pattern ENTREZ_DIRECT -element Count".format(pro,tax),shell=True))
	if num == 0:
		return False
	return True


# Ask the user to enter the aim protein
def ask_the_aim_pro():
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
    print_deco("Finding the protein sequences that related to " + str(pro))
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
	return count 


def pro_tax_count(pro,tax):
	count = int(subprocess.check_output("esearch -db protein -query '{0}[prot]' -organism '{1}' | xtract -pattern ENTREZ_DIRECT -element Count".format(pro,tax),shell=True))
	return count


# Count the number of relative proteins
def count_relative_protein(pro, tax):
    print_deco('Counting the number of relative protein')
    count = pro_tax_count(pro,tax)
    while count >= 1000 :
        pro = input_deco("WARNINGS! There are more than 1000 relative protein! Please check if you have enter the right protein name! \nPlease enter the right name of the protein family:\n")
        ask_the_aim_pro()
    while count == 0 :
        pro = input_deco("WARNINGS! \nThere is no such "+str(pro)+" protein in the database, Please start over!\nPlease enter the right name of the protein family:\n")
        ask_the_aim_pro()


def test_number(n):
    while True:
        try :
            int(n)
            break
        except :
            n = input("I don't think you entered a number, please enter again. ")
    return int(n)


def specific_download(pro, tax):
    # Ask the user whether she/him wants to include 'partial', 'low quality', 'hypothetical', 'predicted', 'isoform' sequences.
    specific = []
    for i in ['partial', 'low quality', 'hypothetical', 'predicted', 'isoform']:
        ans = input_deco('Do you want to include '+ str(i) + ' in your target sequences? (y/n)')
        if ans.lower() == 'n':
            specific.append('NOT ' + str(i).upper() + ' ')
    print_deco('Your choice is: ' + str(specific) + ' Press any key to start downloding sequences.')
    file_check(str(pro) + '_' + str(tax) + ".fasta")
    subprocess.call("esearch -db protein -query "+str(pro)+"[prot] -organism "+str(tax)+''.join(specific)+" | efetch -format fasta > " + str(pro) + '_' + str(tax) + ".fasta", shell = True)
    # Download the data and count the number of sequences and species. 
    file_check(str(pro) + "_" + str(tax) + ".txt")
    subprocess.call('esearch -db protein -query "{0}[prot]" -organism "{1}" | efetch -format docsum | xtract -pattern DocumentSummary -element Organism > {0}_{1}.txt'.format(pro,tax), shell=True)


def seq_num_check(num):
    # If the number of sequences is more than 1000, we can't continue due to the huge data.
    if num >= 1000:
        print_deco("Opps! There are more than 1000 sequences which is too large for us to handle. \nPlease restart this program with other protein.")
        exit()
    else:
        species_num = len(set(open('{0}_{1}.txt'.format(pro,tax)).readlines()))
        # Tell the user what we have found and whether the user want to continue.
        ans = input_deco("I have found " + str(species_num) + " species from what you want to asked.\nDo you want to continue analyze? \nPress any key to continue, type 'q' to exit this program.")
        if ans.lower() == 'q':
            print_deco("Have a nice day! Bye~")
            exit()


def clustalo_align(pro, tax):   
    print_deco("I'm currently aligning the sequences.")
    # Align all the protein sequences using clustalo
    # With '--force' we don't need to check whether the ouput file is exsited.
    subprocess.call('clustalo -i ' + str(pro) + '_' + str(tax) + '.fasta -o ' + str(pro) + '_' + str(tax) + '_align.fasta -v --force', shell =True)


def hydropathy(pro,tax):
    print_deco('Generating Hydropathy plot for alignment')
    subprocess.call('pepwindowall ' + str(pro) + '_' + str(tax) + '_align.fasta -graph svg -goutfile hydropathy', shell = True)
    print_deco('The Hydropathy plot is generated, it has been stored in hydropathy')


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
        l = f2.split("\n") 
        values = list(map(float, l))
        dic = dict(zip(keys, values))
        dic = sorted(dic.items(), key=lambda values:values[1])
        ans = input_deco("The percentage changes are sorted from low to hign.\nDo you want to store them? (y/n) \n\n"+str(dic))
        if ans == 'y':
            # Store the sequeces to the file called acc_bitscore_sorted.txt
            file_check("clustalo_percent_change.txt")
            file = open("clustalo_percent_change.txt", "w")
            for i in dic:
                file.write("\n" + str(i)) 
            file.close()


def consensus_sequence(pro,tax):
    print_deco('Using the EMBOSS Cons to get the consensus sequence from Clustalo output.')
    file_check(str(pro) + '_' + str(tax) + '_consensus.fasta')
    subprocess.call('cons -plurality 0.8 -sequence ' + str(pro) + '_' + str(tax) + '_align.fasta -outseq ' + str(pro) + '_' + str(tax) + '_consensus.fasta', shell=True) 
    print_deco('Consensus has been created now.')


def create_blastdb(pro, tax):
    # Make a blast database from fasta files and save them to an output for user when doing blastp processing. 
    print_deco('I am creating a blastdb from fasta files')
    file_check(str(pro) + "_" + str(tax) + "_database")
    subprocess.call("makeblastdb -in " + str(pro) + "_" + str(tax) + ".fasta -dbtype prot -out " + str(pro) + "_" + str(tax) + "_database", shell = True)
    print_deco('I have now created blast database.')


def pepstats():
    subprocess.call("cp ../acc_bitscore_sorted.fasta .", shell=True) 
    subprocess.call("cp ../acc_bitscore_sorted.txt .", shell=True)
    # Seperated the fasta sequences by using seqresplit
    print_deco("I will get protein statistics from the fasta files.")
    subprocess.call("seqretsplit -sequence acc_bitscore_sorted.fasta -sformat fasta -osformat fasta -auto",shell=True)
    # Uses EMBOSS pepstats to get protein statistics from the fasta files and save them to an output. 
    file_name = []
    file_name = glob.glob("./*.fasta")
    for i in file_name:
        subprocess.call("pepstats -sequence " + str(i) + " -outfile " + str(i) + ".pepstats",shell = True)
    print_deco("I have finish searching ever motif in the sequences you selected! \nYou can review the result in the folder called 'motif'.")


# Check if the prosite files are in user's home directory.
def PROSITE_db_check(prosite_db_directory = '/localdisk/software/EMBOSS-6.6.0/share/EMBOSS/data/PROSITE/prosite.lines'):
    if os.path.exists(prosite_db_directory):
        print_deco("You have PROSITE database in your home directory, we can continue to the next step.")
    else: 
        print('Prosite file \"prosite.lines\" was not found. Attempting to create by downloading .dat & .doc files from EMB=EBI')
        home = os.path.expanduser('~')
        print_deco("I'm checking if you have already download the PROSITE database in your home directory.")
        if os.path.exists(home+"/prosite/prosite.doc") and os.path.exists(home+"/prosite/prosite.dat"):
            print_deco('You have all the database file we needed. We can continue to the next step.')
        else: 
            print("You don't have the files we needed.\nI'm donwloading everything you need now.")
            if not os.path.exists(home+"/prosite"):
                os.mkdir(home+"/prosite")
                subprocess.call("wget http://ftp.ebi.ac.uk/pub/databases/prosite/prosite.doc -O "+home+"/prosite/prosite.doc", shell = True)
                subprocess.call("wget http://ftp.ebi.ac.uk/pub/databases/prosite/prosite.dat -O "+home+"/prosite/prosite.dat", shell = True)
                print_deco('You have all the database file we needed. We can continue to the next step.')
                print_deco("Uses prosextract to process prosite database before we start searching.")
                subprocess.call("prosextract -prositedir PROSITE", shell = True)


# Search motifs from PROSITE database by using patmatmotifs.
def motifs_search():
    file_name = []
    file_name = glob.glob("./*.fasta")
    for i in file_name:
        file_motif = i.replace('.fasta', '.motif')
        subprocess.call("patmatmotifs -sequence " + str(i) + " -auto Y -outfile " + str(file_motif), shell =True)


# Generate a report about the motifs that have been found from PROSITE database.
# The user can choose whether or not send the summary to his/her email.
def summary_email():
    motif_files = [m for m in glob.glob("./*.motif")]
    motif_dict = {}
    HitCount = re.compile('# HitCount:(.*)\n')
    Motif = re.compile('Motif = (.*)\n')
    file_check("summary.txt")
    for m in motif_files:
        with open(m) as f:
            data = f.read()
            motif_num = HitCount.search(data).group(1)
            motifs = Motif.findall(data)
            if int(motif_num) > 0:
                motif_dict[m] = {'no_motifs': motif_num, 'motifs': motifs}
    file = open("summary.txt", 'w')
    file.write('Files Analysed:' + str(len(motif_files)) + ' motifs found in: ' + str(len(motif_files)) + str("\n"))
    for k in motif_dict.keys():
        file.write('In file: ' + str(k) + ' We have found ' + str(motif_dict[k]['no_motifs']) + ' motifs, and the motif(s) is: ' + str(motif_dict[k]['motifs']) + '\n')      
    file.close()
    ans = input_deco("Do you want a summary sent to your email? (y/n) ")
    if ans.lower() == 'y':
        ans = input_deco("Please give me you email address.")
        regular = re.compile(r'[a-zA-Z0-9_@.]')
        while len(ans) != len(regular.findall(ans)) or '@' not in ans or '.' not in ans:
            ans = input_deco("You email seems invalid, please try agin.")
        subprocess.call('mail -s "$DATE" < summary.txt ' + str(ans), shell = True)


def create_pd_from_blastp(file):
    headers=['query_acc.', 'subject_acc.', '%_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bit_score']
    blastdb = pd.read_csv(file, comment = '#', sep ='\t', names = headers)
    return blastdb


def distmat():
    ans = input_deco("Do you want to generate the distance matrix? (y/n)")
    if ans.lower() == 'y':
        subprocess.call('distmat -sequence acc_bitscore_sorted.fasta -outfile distmat.txt -protmethod 0', shell = True)
        print_deco('The distance matrix is now generated and stored in a file called "distmat.txt".')


def PCA_plot(output_path,pro,tax):
    # Open the alignment output file
    file = open(str(output_path) + "/"+str(pro)+"_"+str(tax)+"/"+str(pro)+"_"+str(tax)+"_align.fasta", 'r')
    lines = file.readlines()
    seq_id = []
    seq = []
    temp = ''
    for l in lines:
        if l[0] == '>':
            # head line
            seq_id.append(l.strip())
            seq.append(temp)
            temp = ''
        else:
            # Sequence line
            temp += l.strip()
    file.close()
    seq.pop(0)
    seq_id.pop(-1)
    org = []
    # Separate the name of organism
    for i in seq_id:
        org.append(i.split('[')[1].strip(']'))
    seq_sep = []
    for i in seq:
        seq_sep.append(list(i))
    # Store all the sequences in a dataframe
    align = pd.DataFrame.from_records(seq_sep)
    # 'org' is the label of each sequences
    align['organism'] = org
    align = align.replace('-',np.nan)
    # Convert all the aa into categorical type
    ctype = CategoricalDtype(['I', 'Y', 'G', 'Q', 'L', 'S', 'W', 'M', 'N', 'H', '-', 'V', 'D', 'E', 'F', 'A', 'R', 'K', 'C', 'P', 'T'])
    align_1hot = pd.get_dummies(align.drop('organism', axis=1).astype(ctype))
    X = align_1hot.values
    # sub_labels = list(set(org))
    # Ask if the user want to visualize the sequences variances in a 3D graph
    ans = input_deco("Do you want to visualize the sequences variances in a 3D graph? (y/n)")
    # while ans.lower() != 'y' and ans.lower() != 'n':
    #     ans = input_deco("Do you want to visualize the sequences variances in a 3D graph? (y/n)")
    if ans.lower() == 'y':
        PCA_3D_plot(X)


def PCA_3D_plot(X):
    X_3d = PCA(n_components=3).fit(X).transform(X)
    fig = plt.figure(1, figsize=(8, 6))
    ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=30, azim=300,auto_add_to_figure=False)
    fig.add_axes(ax)
    ax.scatter(X_3d[:,0], X_3d[:,1], X_3d[:,2])
    plt.legend(loc='center left', bbox_to_anchor=[1.1, .5])
    plt.show()   


def distmap_heatmap(pro, tax):
    ans = input_deco("Do you want to generate a heatmap for distant matrix? (y/n)")
    if ans.lower() == 'y':
        subprocess.call('clustalo -i acc_bitscore_sorted.fasta --distmat-out=DISTMAT --full', shell = True)
        file = open("DISTMAT",'r')
        dist = []
        lines = file.readlines()
        for l in lines:
            dist.append(l.rstrip().split(' '))
        # dist = list(filter(None, dist))
        file.close()
        dist_arr = pd.DataFrame(dist[1:])
        dist_arr = dist_arr.replace('', None)
        dist_arr = dist_arr.iloc[:,:-4]
        sns.set(font_scale=0.6)
        sns.clustermap(dist_arr.drop([0],axis = 1).astype(float), cmap='rainbow',yticklabels=dist_arr[0])
        plt.show()



################################################################
# Actually running part
### Part 1. Find out want the user want, and download the data
# Check if the user have installed edirect.
edirect_check()
# Ask the user where to store the ouput files.
output_path = output_path()
os.chdir(output_path)
# Ask what the user want and check if the input is correct.
pro, tax = ask_the_aim_pro()
count_relative_protein(pro, tax)

# Make a directory to start analyzing.
print_deco("Thanks for your cooperation, now I made a directory for you named " + str(pro) + str(tax))
mkdir(str(pro) + "_" + str(tax))
os.chdir('{0}/{1}_{2}'.format(os.getcwd(),pro,tax))

# Download the data and count the number of sequences and species. 
file_check(str(pro) + "_" + str(tax) + ".txt")
specific_download(pro,tax)
# subprocess.call('esearch -db protein -query "{0}[prot]" -organism "{1}" | efetch -format docsum | xtract -pattern DocumentSummary -element Organism > {0}_{1}.txt'.format(pro,tax), shell=True)

# Count the number of sequences and species.
seq_num = len(open(str(pro) + '_' + str(tax) + '.txt').readlines())

# Check if the number of sequences is too large.
seq_num_check(seq_num)

# Download the fasta files to the current directory.
print_deco("Sure! Now I will download the fasta files that related to what you searching for.")
subprocess.call('esearch -db protein -query "{0}[prot]" -organism "{1}" | efetch -format fasta > {0}_{1}.fasta'.format(pro,tax), shell=True)


### Part 2. Analyze the conservation between the sequences we found.

# Using clustalo to align the conservation area between all the sequences.
clustalo_align(pro, tax)
# Make a plot on hydropathy alignment.
hydropathy(pro,tax)
# Plot consecration of a sequence alignment.
plotcon(str(pro) + '_' + str(tax) + '_align.fasta')
# Calculate the percentage changes.
clustalo_percent_change(pro, tax)
# Using the EMBOSS Cons to get the consensus sequence from Clustalo output.
consensus_sequence(pro,tax)

# Starting to do blastdb
create_blastdb(pro, tax)
# Using blastp to blast the protein that the user is interested in.
print_deco('Using blastp to search this file: ' + str(pro) + '_' + str(tax) + '_consensus.fasta in blastdb')
subprocess.call('blastp -db ' + str(pro) + '_' + str(tax) + '_database -query ' + str(pro) + '_' + str(tax) + '_consensus.fasta -outfmt 7 -out ' + str(pro) + '_' + str(tax) + '_blast.txt', shell = True)

# Creat a dictionary that contains the acc and the bit score.
blastp_result = {}
line = []
line = [i for i in open(str(pro) + '_' + str(tax) + '_blast.txt', 'r') if '#' not in i]
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
subprocess.call('/localdisk/data/BPSM/Assignment2/pullseq -i ' + str(pro) + '_' + str(tax) + '_align.fasta -n acc_bitscore_sorted.txt > acc_bitscore_sorted.fasta -v', shell=True)

# Use plotcon to plot the similar sequences alignment result 
ans = input_deco("Do you want to see the plot of blast result by plotcon? (y/n)")
if ans.lower() == 'y':
    plotcon("acc_bitscore_sorted.fasta")
# Generate a distance matrix.
distmat()


### Part 3. Scan protein sequence(s) of interest with motifs from the PROSITE database
# Make a directory to store the seperate sequences files
mkdir("motif")
os.chdir("motif")
pepstats()

# Check it the user have PROSITE database in the home directory.
PROSITE_db_check()
# Search motifs from PROSITE database by using patmatmotifs.
motifs_search()
# Generate a summary result about motifs and ask if the user wants a copy send to his/her email.
summary_email()

### Part 4. Extra ploting (PCA)(Heatmap).
PCA_plot(output_path,pro,tax)
distmap_heatmap(pro,tax)

################################################################
# I tried to use different classifier to classify the protein sequences.
# However the training data is too little, thus the classifier seems not very accurate

# from sklearn.metrics import adjusted_rand_score
# from sklearn.ensemble import RandomForestClassifier
# for ii in range(2, len(X)):
#     X_pca = PCA(n_components = ii).fit_transform(X)
#     rf = RandomForestClassifier(n_estimators=100, random_state=0, oob_score=True)
#     rf.fit(X_pca, org)
#     y_pred = np.argmax(rf.oob_decision_function_, axis=1)
#     print('Random Forests on {}D PCA:\nARI: {}, Mean Accuracy: {}\n\n'.
#           format(ii, adjusted_rand_score(org, y_pred), rf.oob_score_))
# rf = RandomForestClassifier(n_estimators=100, random_state=0, oob_score=True)
# rf.fit(X, org)
# y_pred = np.argmax(rf.oob_decision_function_, axis=1)
# print('Random Forest on Full Data\nARI: {}, Mean Accuracy: {}\n\n'.
#       format(adjusted_rand_score(org, y_pred), rf.oob_score_))


# from sklearn.naive_bayes import GaussianNB
# for ii in range(2, len(X)):
#     X_pca = PCA(n_components = ii).fit_transform(X)
#     gnb = GaussianNB()
#     gnb.fit(X_pca, org)
#     y_pred = gnb.predict(X_pca)
#     print('Gaussian Naive Bayes on {}D PCA:\nARI: {}, Mean Accuracy: {}\n\n'.
#           format(ii, adjusted_rand_score(org, y_pred), gnb.score(X_pca, org)))
# gnb = GaussianNB()
# gnb.fit(X, org)
# y_pred = gnb.predict(X)
# print('Gaussian Naive Bayes on Full Data\nARI: {}, Mean Accuracy: {}\n\n'.
#       format(adjusted_rand_score(org, y_pred), gnb.score(X, org)))


# from sklearn.cluster import KMeans
# for ii in range(2, len(X)):
#     X_pca = PCA(n_components = ii).fit_transform(X)
#     km = KMeans(n_clusters=6, random_state=0)  
#     km.fit(X_pca)
#     y_pred = km.labels_
#     print('k-means with 6 clusters on {}D PCA:\nARI: {}, Inertia: {}\n\n'.
#           format(ii, adjusted_rand_score(org, y_pred), km.inertia_))


print_deco("This is the end of this programe. \nHope your experience is good. \nThanks for using this programe.\nHave a nice day~~BYE")









