import os
import glob
import copy
import numpy as np
import random

'''
Python implementation of the Smith-Waterman algorithm for local alignment of protein sequences, specifically tailored to read
BLOSUM matricies and .fa protein amino acid sequence files.

Based largely on code by Ryan Boehning, https://gist.github.com/radaniba/11019717
'''

class similarity_matrix_class:
    '''
    Quick similarity matrix class to keep track of the scoring matrix and make searching through it to find similarities faster
    '''
    def __init__(self):
        self.index = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'B':20,'Z':21,'X':22,'*':23}
        self.matrix = np.zeros((24,24))

def read_blosum(filename):
    '''
    Method that takes in a filename and reads it into a BLOSUM matrix class
    
    Input: filename of BLOSUM matrix
    Output: similarity matrix class
    '''
    similarity_matrix = similarity_matrix_class()
    with open(filename, "r") as f:
        for i, line in enumerate(f):
            if i > 6: # This is the definition of not being robust
                int_list = line.split()
                for j, char in enumerate(int_list):
                    similarity_matrix.matrix[i - 7,j] = int(char) # as is this
    return similarity_matrix

def read_matio(filename):
    '''
    Method that takes in a filename and reads it into a BLOSUM matrix class
    
    Input: filename of BLOSUM matrix
    Output: similarity matrix class
    '''
    similarity_matrix = similarity_matrix_class()
    with open(filename, "r") as f:
        for i, line in enumerate(f):
            if i > 2: # This is the definition of not being robust
                int_list = line.split()
                for j, char in enumerate(int_list):
                    similarity_matrix.matrix[i - 3,j] = int(char) # as is this
    return similarity_matrix

def read_pam(filename):
    '''
    Method that takes in a filename and reads it into a BLOSUM matrix class
    
    Input: filename of BLOSUM matrix
    Output: similarity matrix class
    '''
    similarity_matrix = similarity_matrix_class()
    with open(filename, "r") as f:
        for i, line in enumerate(f):
            if i > 9: # This is the definition of not being robust
                int_list = line.split()
                for j, char in enumerate(int_list):
                    similarity_matrix.matrix[i - 10,j] = int(char) # as is this
    return similarity_matrix


def read_sequences(dir):
    """
    Read in all of the sequences from the given directory.

    Input: directory
    Output: list of sequences as strings 
    """

    sequences = []
    # iterate over each .fa file in the given directory
    for filepath in glob.iglob(os.path.join(dir, "*.fa")):
        with open(filepath, "r") as f:
            sequence = ''
            for i, line in enumerate(f):
                if i > 0: # Again, this is just bad code
                    sequence += line
            sequence = sequence.replace('\n','').replace('\r','')
            sequences.append(sequence)
    print("Read in %d sequences"%len(sequences))
    return sequences

def read_sequence(filename):
    '''
    Read in a specific sequence.
    
    Input: filepath
    Output: sequence as string
    '''
    with open(filename, "r") as f:
        sequence = ''
        for i, line in enumerate(f):
            if i > 0: # Again, this is just bad code
                sequence += line
        sequence = sequence.replace('\n','').replace('\r','')
    return sequence


def read_pairs(filename):
    '''
    Read pairs linewise from a txt file (assumes file contains direct path to file, including directories)
    
    Input: filename
    Output: list of paired sequences as strings
    '''
    pairs = []
    with open(filename, "r") as f:
        for i, line in enumerate(f):
            pair = []
            to_compare = line.split()
            for seq_path in to_compare:
                sequence = read_sequence(seq_path)
                pair.append(sequence)
            pairs.append(pair)
    print('successfully read all sequences in ' + str(filename) + '.')   
    return pairs



def get_similarity(a, b, matrix):
    '''
    Using the given scoring matrix, do a quick lookup to see the similarity between given input residue values (as strings)
    
    Input: Two residues as strings, a corresponding similarity matrix.
    Output: A raw integer score that signifies the similarity between the residues.
    '''
    
    # In the case of some error, assume that there is no match and give the similarity a value of -5
    a_index = 23
    b_index = 23
    
    for name, index in matrix.index.items():
        if a == name:
            a_index = index
        if b == name:           
            b_index = index
    score = matrix.matrix[a_index][b_index]
    return score

def create_score_matrix(a,b,path_matrix,matrix, gap_opening_penalty, gap_extension_penalty):
    '''
    Given your sequences, construct a matrix of score values and keep track of the path taken to get there.
    More of a wrapper function for get_score()
    
    Inputs: Sequences a & b, a blank path matrix (same size as score matrix), and a similarity matrix
    Outputs: Score matrix, tuple of highest scoring cell as well as coordinates
    '''
    rows = len(a) + 1
    cols = len(b) + 1          
    score_matrix = np.zeros((rows, cols))
    # Hiscore keeps track of the highest score, then the coordinates of that score
    hiscore = (0,0,0)
    
    # Iterate through all cells, building the score matrix
    for i in range(1, rows):
        for j in range(1, cols):
            score = get_score(score_matrix, path_matrix, i, j, a, b, matrix, gap_opening_penalty, gap_extension_penalty)
            if score > hiscore[0]:
                hiscore = (score,i,j)
            score_matrix[i][j] = score

    return score_matrix, hiscore


def create_path_matrix(a,b):
    # Pretty redundant but ehh whatever
    path_matrix = np.empty((len(a) + 1,len(b) + 1),dtype=str)
    return path_matrix
    
def get_score(score_matrix, path_matrix, i, j, a, b, matrix, gap_opening_penalty, gap_extension_penalty):
    '''
    Using the similarity matrix, get a similarity score and account for matches, gap openings & extensions.
    Fill in the score matrix, and update the path matrix with the route to get there. Lastly, return the highest score for that cell.
    
    Inputs: Score matrix, path matrix, current indicies, the sequences being compared, and the similartiy matrix
    Outputs: Updating the path matrix and returning the maximum score for that cell
    '''
    similarity = get_similarity(a[i - 1], b[j - 1], matrix)
    
    # Determine the scores resulting from taking any particular path to that cell
    # Score for a match, coming from the cell diagonally above
    diag_score = score_matrix[i - 1][j - 1] + similarity
                             
    # Coming from a gap: check to see if you're opening or extending it
    # Coming from the cell above:
    if path_matrix[i - 1][j] == 'd': # you're opening a gap
        up_score = score_matrix[i - 1][j] - gap_opening_penalty
    else: # you're extending a gap
        up_score = score_matrix[i - 1][j] - gap_extension_penalty
    # coming from the cell to the left
    if path_matrix[i][j - 1] == 'd': # you're opening a gap
        left_score = score_matrix[i][j - 1] - gap_opening_penalty
    else: # you're extending a gap
        left_score = score_matrix[i][j - 1] - gap_extension_penalty                           

    # Choose the highest score to determine which path should be taken to get to that cell                         
    max_score = max(0, diag_score, up_score, left_score)    
                   
    # Keep track of how you got to the current cell in the path matrix
    if max_score == 0:
        path_matrix[i][j] = 'n'    
    elif max_score == diag_score:
        path_matrix[i][j] = 'd'
    elif max_score == up_score:
        path_matrix[i][j] = 'u'
    elif max_score == left_score:
        path_matrix[i][j] = 'l'                             
    return max_score    





def get_pair_scores(pairs, gap_opening_penalty, gap_extension_penalty, matrix):
    pair_hiscores = [0]*len(pairs)
    for index, pair in enumerate(pairs):
        a = pair[0]
        b = pair[1]
        hiscore = (0,0,0)
        
        path_matrix = create_path_matrix(a,b)
        score_matrix, hiscore = create_score_matrix(a, b, path_matrix, matrix, gap_opening_penalty, gap_extension_penalty)
        
        pair_hiscores[index] = hiscore[0]
        #print('Index: ' + str(index) + '\n \t \t hiscore: ' + str(hiscore))
    return pair_hiscores

def get_pair_scores_adjusted(pairs, gap_opening_penalty, gap_extension_penalty, matrix):
    pair_hiscores = [0]*len(pairs)
    for index, pair in enumerate(pairs):
        a = pair[0]
        b = pair[1]
        hiscore = (0,0,0)
        
        path_matrix = create_path_matrix(a,b)
        score_matrix, hiscore = create_score_matrix(a, b, path_matrix, matrix, gap_opening_penalty, gap_extension_penalty)
        
        pair_hiscores[index] = hiscore[0] / min(len(a), len(b))
    return pair_hiscores


def get_best_alignments(pairs, gap_opening_penalty, gap_extension_penalty, matrix):
    pair_hiscores = [0]*len(pairs)
    for index, pair in enumerate(pairs):
        a = pair[0]
        b = pair[1]
        hiscore = (0,0,0)
        
        path_matrix = create_path_matrix(a,b)
        score_matrix, hiscore = create_score_matrix(a, b, path_matrix, matrix, gap_opening_penalty, gap_extension_penalty)
        
        
        pair_hiscores[index] = hiscore[0]
        print('hiscore of index ' + str(index) + ' is ' + str(hiscore))
    return pair_hiscores

def traceback(a, b, score_matrix, path_matrix, hiscore):
    done = False
    a_seq_alignment = []
    b_seq_alignment = []
    x, y = hiscore[1], hiscore[2]

    while not done:
        diag = score_matrix[x - 1][y - 1]
        up = score_matrix[x - 1][y]
        left = score_matrix[x][y - 1]
        if diag >= up and diag >= left:
            if diag != 0:
                a_seq_alignment.append(a[x - 1])
                b_seq_alignment.append(b[y - 1])
                y -= 1
                x -= 1
            else:
                done = True
        elif up > diag and up >= left:
            if up != 0:
                a_seq_alignment.append(a[x - 1])
                b_seq_alignment.append('*')
                x -= 1
            else:
                done = True
        elif left > diag and left > up:
            if left != 0:
                b_seq_alignment.append(b[y - 1])
                a_seq_alignment.append('*')
                y -= 1
            else:
                done = True
        
    return ''.join(reversed(a_seq_alignment)), ''.join(reversed(b_seq_alignment))
        
def get_smaller_pairs(pairs, gap_opening_penalty, gap_extension_penalty, matrix):
    smaller_pairs = []*len(pairs)
    for index, pair in enumerate(pairs):
        small_pair = []
        a = pair[0]
        b = pair[1]
        hiscore = (0,0,0)
        
        path_matrix = create_path_matrix(a,b)
        score_matrix, hiscore = create_score_matrix(a, b, path_matrix, matrix, gap_opening_penalty, gap_extension_penalty)
        

        a_trace, b_trace = traceback(a, b, score_matrix, path_matrix, hiscore)
        small_pair.append(a_trace)
        small_pair.append(b_trace)
        
        smaller_pairs.append(small_pair)
        print('smaller_pairs is')
        print(smaller_pairs)
    return smaller_pairs


def test_gap_penalties(pospairs_hiscores, negpairs_hiscores):
    '''
    Iterate through all this damn s*** and get lists of tuples of varying gap opening / extension penalties and corresponding scores
    '''
    # Initialize matricies to keep track of the best scores obtained for each run
    pospairs_hiscores_all = []
    negpairs_hiscores_all = [] 
    
    for gap_opening_penalty in range(1,20):
        for gap_extension_penalty in range(1,5):
            print('Testing gap opening / extension penalty of (' + str(gap_opening_penalty) + ',' + str(gap_extension_penalty) + ').')
            
            # For each positive pair, compute the highest scoring alignment
            pospairs_toadd = get_pair_scores(pospairs, gap_opening_penalty, gap_extension_penalty, matrix)
            pospairs_hiscores_all.append((gap_opening_penalty, gap_extension_penalty, pospairs_toadd))
            
            # For each negative pair, compute the highest scoring alignment
            negpairs_toadd = get_pair_scores(negpairs, gap_opening_penalty, gap_extension_penalty, matrix)
            negpairs_hiscores_all.append((gap_opening_penalty, gap_extension_penalty, negpairs_toadd))
     
    return pospairs_hiscores_all, negpairs_hiscores_all
    


def get_ROC_data(pospairs_hiscores, negpairs_hiscores):
    minscore = min(negpairs_hiscores)
    maxscore = max(pospairs_hiscores)
    ROC_data = np.zeros((len(np.arange(minscore,maxscore,maxscore / 200)),2))
    score_thresholds = np.arange(minscore,maxscore,maxscore / 200)
    for index in range(0,len(np.arange(minscore,maxscore,maxscore / 200))):
        ROC_data[index][0] = fraction_above_threshold(negpairs_hiscores,score_thresholds[index])
        ROC_data[index][1] = fraction_above_threshold(pospairs_hiscores,score_thresholds[index])
    return ROC_data     

def fraction_above_threshold(hiscores,score_threshold):
    above_threshold = 0
    for index in range(0,len(hiscores)):
        if hiscores[index] > score_threshold:
            above_threshold += 1
    above_threshold = above_threshold / len(hiscores)
    return above_threshold


def get_cutoff(pair_hiscores, cutoff):
    '''
    Determine the score threshold at which you have (x%) of your scores above the cutoff.
    '''
    # Iterate from the maximum score down to find the cutoff
    for threshold in np.arange(min(pair_hiscores), max(pair_hiscores) + 2 * (max(pair_hiscores) / 100), max(pair_hiscores) / 100)[::-1]:
        above_cutoff = 0
        for score in pair_hiscores:
            if score > threshold:
                above_cutoff += 1
        above_cutoff = above_cutoff / len(pair_hiscores)
        if above_cutoff >= cutoff:
            return threshold  
    return threshold


def get_fpr_array(negpairs_hiscores_all, cutoff_array):
    '''
    Determine the false positive rate given a list of tuples that contains: (gap opening penalty, gap extension penalty, score threshold)
    
    Input: all of the negative pair hiscores, cutoff array
    Output:
    '''
    fpr_array = np.append(cutoff_array,np.zeros([len(cutoff_array),1]),1)
    min_fpr = 1
    min_fpr_index = 0
    for index in range(0,len(fpr_array)):
        for jdex in range(0,len(negpairs_hiscores_all)):
            if fpr_array[index][0] == negpairs_hiscores_all[jdex][0] and fpr_array[index][1] == negpairs_hiscores_all[jdex][1]:
                above_threshold = 0
                for score in range(0,len(negpairs_hiscores_all[jdex][2])):
                    if negpairs_hiscores_all[jdex][2][score] > fpr_array[index][2]:
                        above_threshold += 1
                above_threshold = above_threshold / len(negpairs_hiscores_all[jdex][2])
                # print(str(above_threshold) + ' percent of the matches are above the score threshold')
                fpr_array[index][3] = above_threshold
                if above_threshold < min_fpr:
                    min_fpr = above_threshold
                    min_fpr_index = index
                    #print('found a lower minimum threshold at' + str(fpr_array[index]))
    
    print('The lowest false positive rate, ' + str(min_fpr) + ' that still achieves a 70% true positive rate occurs when the score threshold is ' 
              + str(fpr_array[min_fpr_index][2]) + ', the gap opening penalty is ' + str(fpr_array[min_fpr_index][0]) + 
                ' and the gap extension penalty is ' + str(fpr_array[min_fpr_index][1]))
    return fpr_array

def top_n(a,N):
    return np.argsort(a)[::-1][:N]



'''
And now, we actually call some of the functions...
This repo is going private as soon as the assignment is over, this code is too hideous for the outside world
'''




# Initialize your scoring matrix and path matrix
blosum_matrix = read_blosum('BLOSUM50')
blosum_matrix_50 = read_blosum('BLOSUM50')
blosum_matrix_62 = read_blosum('BLOSUM62')

matio_matrix = read_matio('MATIO')
pam100_matrix = read_pam('PAM100')
pam250_matrix = read_pam('PAM250')


# Read in all pairs of sequences to be alligned
pospairs = read_pairs('Pospairs.txt')
negpairs = read_pairs('Negpairs.txt')

#pospair_hiscores_50 = get_pair_scores_adjusted(pospairs, 9, 4, blosum_matrix_50)
#pospair_hiscores_62 = get_pair_scores_adjusted(pospairs, 9, 4, blosum_matrix_62)

#negpair_hiscores_50 = get_pair_scores_adjusted(negpairs, 9, 4, blosum_matrix_50)
#negpair_hiscores_62 = get_pair_scores_adjusted(negpairs, 9, 4, blosum_matrix_62)

#roc_50 = get_ROC_data(pospair_hiscores_50,negpair_hiscores_50)
#roc_62 = get_ROC_data(pospair_hiscores_62,negpair_hiscores_62)
    
#pospair_hiscores_matio = get_pair_scores(pospairs, 9, 4, matio_matrix)
#negpair_hiscores_matio = get_pair_scores(negpairs, 9, 4, matio_matrix)
#roc_matio = get_ROC_data(pospair_hiscores_matio, negpair_hiscores_matio, min(negpair_hiscores_matio), max(pospair_hiscores_matio))

#pospair_hiscores_pam100_adjusted = get_pair_scores_adjusted(pospairs, 9, 4, pam100_matrix)
#negpair_hiscores_pam100_adjusted = get_pair_scores_adjusted(negpairs, 9, 4, pam100_matrix)
#roc_pam100_adjusted = get_ROC_data(pospair_hiscores_pam100_adjusted, negpair_hiscores_pam100_adjusted, min(negpair_hiscores_pam100_adjusted), max(pospair_hiscores_pam100_adjusted))

#pospair_hiscores_pam250 = get_pair_scores(pospairs, 9, 4, pam250_matrix)
#negpair_hiscores_pam250 = get_pair_scores(negpairs, 9, 4, pam250_matrix)
#roc_pam250 = get_ROC_data(pospair_hiscores_pam250, negpair_hiscores_pam250, min(negpair_hiscores_pam250), max(pospair_hiscores_pam250))

#pospair_hiscores_pam100 = get_pair_scores(pospairs, 9, 4, pam100_matrix)
#negpair_hiscores_pam100 = get_pair_scores(negpairs, 9, 4, pam100_matrix)





# Make a new, mutable copy of the pam100 matrix
new_pam = copy.deepcopy(pam100_matrix)
new_pam.matrix = new_pam.matrix.astype(np.float64)


# Define number of iterations, weight at which you modify scores in the matrix
iterations = 300
alpha = 2

# When I initialized the code, the TPR was 1.98, so I use that here as my baseline
tpr_max = 1.98

# Get smaller pairs for alignment purposes
pospairs_small = get_smaller_pairs(pospairs, 7,4, new_pam)
negpairs_small = get_smaller_pairs(negpairs, 7,4, new_pam)

for i in range(0, iterations):
    print('iteration: ' + str(i + 1))
    old_pam = copy.deepcopy(new_pam)
    randint = random.randint(0,23)
    new_pam.matrix[randint][randint] = new_pam.matrix[randint][randint] + alpha * random.uniform(-3,3)
    
    pospair_small_hiscores = get_pair_scores(pospairs_small, 9, 4, new_pam)
    negpair_small_hiscores = get_pair_scores(negpairs_small, 9, 4, new_pam)
    
    # Calculate score cutoffs corresponding to FPRs of 0.0, 0.1, 0.2, and 0.3, and their corresponding TPRs
    fpr_cutoffs = [0.0] * 4
    fpr = [0.0] * 4
    tpr = [0.0] * 4       
    tpr_sum = 0
          
    for index in range(0,4):
        fpr_cutoffs[index] = get_cutoff(negpair_small_hiscores, 0.1 * index)
        fpr[index] = fraction_above_threshold(negpair_small_hiscores, fpr_cutoffs[index])
        tpr[index] = fraction_above_threshold(pospair_small_hiscores, fpr_cutoffs[index])
        tpr_sum += tpr[index]
    
    print('TPR max: ' + str(tpr_max))
    print('TPR now: ' + str(tpr_sum))
    
    if tpr_sum < tpr_max:
        new_pam = old_pam
    else:
        tpr_max = tpr_sum
    
print('After iterating ' + str(iterations) + ' times, TPR is now ' + str(tpr_max))    
    
small_roc = get_ROC_data(pospair_small_hiscores, negpair_small_hiscores)
    






