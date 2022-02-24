# Read in FASTA
# file = open(str(pro)+"_"+str(tax)+"_alian.fasta", 'r')

# clustalo -i glucose-6-phosphatase_Aves.fasta -o glucose-6-phosphatase_Aves_align.fasta -v --force -seqnos
# esearch -db protein -query "glucose-6-phosphatase[prot]" -organism "Aves" | efetch -db protein -format docsum

from script.py import *


# def PCA():
#     file = open("/localdisk/home/s2160628/CW2/glucose-6-phosphatase_Aves/glucose-6-phosphatase_Aves_align.fasta", 'r')
#     lines = file.readlines()
#     seq_id = []
#     seq = []
#     temp = ''
#     for l in lines:
#         if l[0] == '>':
#             # head line
#             seq_id.append(l.strip())
#             seq.append(temp)
#             temp = ''
#         else:
#             # Sequence line
#             temp += l.strip()
#     file.close()
#     seq.pop(0)
#     seq_id.pop(-1)
#     org = []
#     for i in seq_id:
#         org.append(i.split('[')[1].strip(']'))
#     seq_sep = []
#     for i in seq:
#         seq_sep.append(list(i))
#     align = pd.DataFrame.from_records(seq_sep)
#     align['organism'] = org
#     align = align.replace('-',np.nan)

#     ctype = CategoricalDtype(['I', 'Y', 'G', 'Q', 'L', 'S', 'W', 'M', 'N', 'H', '-', 'V', 'D', 'E', 'F', 'A', 'R', 'K', 'C', 'P', 'T'])
#     align_1hot = pd.get_dummies(align.drop('organism', axis=1).astype(ctype))
#     X = align_1hot.values
#     sub_labels = list(set(org))

# def PCA_3D(X):
#     # sub_cats = [sub_labels[label] for label in sub_labels]
#     X_3d = PCA(n_components=3).fit(X).transform(X)
#     fig = plt.figure(1, figsize=(8, 6))
#     ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=30, azim=300,auto_add_to_figure=False)
#     fig.add_axes(ax)
#     ax.scatter(X_3d[:,0], X_3d[:,1], X_3d[:,2])
#     plt.legend(loc='center left', bbox_to_anchor=[1.1, .5])
#     plt.show()


# def PCA_2D():
#     pca = PCA(n_components=2)
#     X_2d = pca.fit_transform(X)
#     plt.figure(figsize=(12,8))
#     plt.scatter(X_2d[:,0], X_2d[:,1], alpha=.5)
#     plt.axis('equal')
#     plt.legend(loc='center left', scatterpoints=3, bbox_to_anchor=[1.01, 0.5])
#     plt.title('Labelled data in PCA space')
#     plt.xlabel('PC1')
#     plt.ylabel('PC2')
#     top_plot = plt.gca()
#     plt.show()

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









# pca = PCA().fit(X)

# fig, ax = plt.subplots(ncols=2, figsize=(20,6))
# ax1, ax2 = ax.ravel()

# ratio = pca.explained_variance_ratio_
# ax1.bar(range(len(ratio)), ratio, color='purple', alpha=0.8)
# ax1.set_title('Explained Variance Ratio PCA', fontsize=20)
# ax1.set_xticks(range(len(ratio)))
# ax1.set_xticklabels(['PC {}'.format(i+1) for i in range(len(ratio))])
# ax1.set_ylabel('Explained Variance Ratio')

# # ratio[0]=0
# ratio = pca.explained_variance_ratio_
# ax2.plot(np.cumsum(ratio), 'o-')

# ax2.set_title('Cumulative Sum of Explained Variance Ratio PCA', fontsize=20)

# ax2.set_ylim(0,1.1)
# ax2.set_xticks(range(len(ratio)))
# ax2.set_xticklabels(['PC {}'.format(i+1) for i in range(len(ratio))])
# ax2.set_ylabel('Cumulative Sum of Explained Variance Ratio')
# plt.show()



def PCA_plot(pro,tax):
    # Open the alignment output file
    file = open(os.getcwd() + "/"+str(pro)+"_"+str(tax)+"/"+str(pro)+"_"+str(tax)+"_align.fasta", 'r')
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
        X_3d = PCA(n_components=3).fit(X).transform(X)
        fig = plt.figure(1, figsize=(8, 6))
        ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=30, azim=300,auto_add_to_figure=False)
        fig.add_axes(ax)
        ax.scatter(X_3d[:,0], X_3d[:,1], X_3d[:,2])
        # plt.legend(loc='center left', bbox_to_anchor=[1.1, .5])
        plt.show()  


# sns.heatmap(flights, cmap='rainbow', ax=ax).set(title='Heatmap for distance matrix')

PCA_plot(os.getcwd,'glucose-6-phosphatase','Aves')


# distmat -sequence acc_bitscore_sorted.fasta -outfile distmat.txt -protmethod 0
# clustalo -i acc_bitscore_sorted.fasta --distmat-out=DISTMAT --full
# clustalo -i acc_bitscore_sorted.fasta -o my-out-seqs.fa -v --distmat-out=output.distmat
# clustalo -i acc_bitscore_sorted.fasta --distmat-out=DISTMAT --full


def distmap_heatmap(output_path,pro,tax):
    ans = input_deco("Do you want to generate a heatmap for distant matrix? (y/n)")
    if ans.lower() == 'y':
subprocess.call('clustalo -i '+str(output_path) + "/"+str(pro)+"_"+str(tax)+"/acc_bitscore_sorted.fasta --distmat-out=DISTMAT --full", shell = True)
        with open("DISTMAT").read() as file:
            dist = []
            lines = file.readlines()
            for l in lines:
                if len(l.split('\t')) != 1:
                    dist.append(l.split('\t'))
            dist_arr = pd.DataFrame(dist)
            sns.clustermap(dist_arr, row_colors='rainbow')
            plt.show()


clustalo -i acc_bitscore_sorted.fasta --distmat-out=DISTMAT --full

subprocess.call('clustalo -i '+str(output_path) + "/"+str(pro)+"_"+str(tax)+"/acc_bitscore_sorted.fasta --distmat-out=DISTMAT --full", shell = True)

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

dist_clean = []
for i in dist:
    if '' in i:
        dist_clean.append(i.remove(''))
    else:
        dist_clean.append(i)
dist_clean = [x for x in dist_clean if x is not None]
dist_arr = pd.DataFrame(dist_clean)


a = sns.clustermap(dist_arr.drop([0],axis = 1).astype(float), cmap='rainbow',yticklabels=dist_arr[0])
plt.show()


num_elements = len(dist_arr)
X_Tick_List = []
X_Tick_Label_List=[]
for item in range (0,num_elements):
    X_Tick_List.append(dist_arr.drop([0],axis = 1)[item])
    X_Tick_Label_List.append(dist_arr[0](item+1))


sns.heatmap(dist_arr.drop([0],axis = 1).astype(float), cmap='rainbow')
plt.xticks(ticks=X_Tick_List,labels=X_Tick_LabeL_List, rotation=25,fontsize=8) 
plt.show()      


a='XP_040548341.1 0.000000 0.000000 0.002315 0.002315 0.094972 0.078212 0.170391 0.173184 0.175978 0.170391 0.175978 0.243243 0.097765 0.170391 0.120112 0.178771 0.103352 0.170391 0.050926 0.108939 0.114525 0.108939 0.108939 0.106145 0.106145 0.108939 0.108939 0.111732 0.111732 0.106145 0.111732 0.108939 0.108939 0.108939 0.108939 0.108939 0.106145 0.106145 0.106145 0.106145 0.175978 0.046296 0.173184 0.285714 0.178273 0.131285 0.131285 0.181564 0.173184 0.162011 0.094972 0.097765 0.170391 0.173184 0.239808 0.167598 0.106145 0.041667 0.106145 0.131285 0.114525 0.170391 0.175978 0.178273 0.108939 0.170391 0.125698 0.168067 0.111732 0.094972 0.167464 0.120112 0.120112 0.103352 0.103352 0.194444 0.097765 0.096317 0.120112 0.163017 0.172805 0.109551 0.124646 0.123209 0.101983 0.118980 0.136872 0.107649 0.104816 0.096591 0.107649 0.099150 0.133144 0.113314 0.096317 0.127479 0.118980 0.172805 0.117318 0.112360'

dist_arr.to_csv(index=False)




# 1

def __init__(self):
    self.mean = None
    self.cov = None
@staticmethod
def reX(Ztrn):
    return Ztrn[..., np.newaxis] if Ztrn.ndim == 2 else Ztrn
def fit(self, Ztrn):
    Ztrn = self.reX(Ztrn)
    self.mean = Ztrn.mean(0)
    self.cov = np.einsum('ijk,ikj->jk', Ztrn-self.mean, Ztrn-self.mean) / (Ztrn.shape[0]-1)
def prob(self, Ztrn):
    Ztrn = self.reX(Ztrn)
    a = (2*np.pi)**(-self.mean.shape[0]/2)*np.linalg.det(self.cov)**(-1/2)
    b = np.exp((-1/2)*np.einsum('ijk,jl,ilk->ik',Ztrn-self.mean,np.linalg.inv(self.cov),Ztrn-self.mean))
    return a*b

    
    mod = Normal()
    mod.fit(Ztrn)
    print(mod.mean)
    print(mod.cov)
    xx,xy = np.meshgrid(*np.linspace(0,50,318)[np.newaxis,...].repeat(2,0))
    m = np.array([xx,xy]).transpose((1,2,0)).reshape(-1,2,1)
    p = mod.prob(m)

cs = plt.contour(xx,xy,p.reshape(318,318))
plt.clabel(cs)
plt.scatter(*Ztrn.T, alpha=0.5, s  = 5)
plt.grid()
plt.gca().set_aspect('equal')
plt.colorbar(cs) # Add a colorbar to a plot

plt.title('Contours Plot for Ztrn')
plt.xlabel=('A4')
plt.ylabel=('A7')
plt.figure(figsize=(50,50))
plt.show()