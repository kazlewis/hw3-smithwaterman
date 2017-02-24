import os
import glob
import numpy as np

'''
Python implementation of the Smith-Waterman algorithm for local alignment of protein sequences, specifically tailored to read
BLOSUM matricies and .fa protein amino acid sequence files.

Based largely on code by Ryan Boehning, https://gist.github.com/radaniba/11019717
'''

class blosum_matrix_class:
    '''
    Quick blosum matrix class to keep track of the scoring matrix and make searching through it to find similarities faster
    '''
    def __init__(self):
        self.index = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'B':20,'Z':21,'X':22,'*':23}
        self.matrix = np.zeros((24,24))

def read_blosum(filename):
    '''
    Method that takes in a filename and reads it into a BLOSUM matrix class
    
    Input: filename of BLOSUM matrix
    Output: BLOSUM matrix class
    '''
    blosum_matrix = blosum_matrix_class()
    with open(filename, "r") as f:
        for i, line in enumerate(f):
            if i > 6: # This is the definition of not being robust
                int_list = line.split()
                for j, char in enumerate(int_list):
                    blosum_matrix.matrix[i - 7,j] = int(char) # as is this
    return (blosum_matrix)

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



def get_similarity(a, b, blosum_matrix):
    '''
    Using the given BLOSUM matrix, do a quick lookup to see the similarity between given input residue values (as strings)
    
    Input: Two residues as strings, a corresponding BLOSUM matrix.
    Output: A raw integer score that signifies the similarity between the residues.
    '''
    
    # In the case of screwed up residues, assume that there is no match and give the similarity a value of -5
    a_index = 23
    b_index = 23
    
    for name, index in blosum_matrix.index.items():
        if a == name:
            a_index = index
        if b == name:           
            b_index = index
    score = blosum_matrix.matrix[a_index][b_index]
    return score

def create_score_matrix(a,b,path_matrix,blosum_matrix, gap_opening_penalty, gap_extension_penalty):
    '''
    Given your sequences, construct a matrix of score values and keep track of the path taken to get there.
    More of a wrapper function for get_score()
    
    Inputs: Sequences a & b, a blank path matrix (same size as score matrix), and a similarity BLOSUM matrix
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
            score = get_score(score_matrix, path_matrix, i, j, a, b, blosum_matrix, gap_opening_penalty, gap_extension_penalty)
            if score > hiscore[0]:
                hiscore = (score,i,j)
            score_matrix[i][j] = score

    return score_matrix, hiscore


def create_path_matrix(a,b):
    # Pretty redundant but ehh whatever
    path_matrix = np.empty((len(a) + 1,len(b) + 1),dtype=str)
    return path_matrix
    
def get_score(score_matrix, path_matrix, i, j, a, b, blosum_matrix, gap_opening_penalty, gap_extension_penalty):
    '''
    Using the BLOSUM matrix, get a similarity score and account for matches, gap openings & extensions.
    Fill in the score matrix, and update the path matrix with the route to get there. Lastly, return the highest score for that cell.
    
    Inputs: Score matrix, path matrix, current indicies, the sequences being compared, and the BLOSUM scoring matrix
    Outputs: Updating the path matrix and returning the maximum score for that cell
    '''
    similarity = get_similarity(a[i - 1], b[j - 1], blosum_matrix)
    
    # Determine the scores resulting from taking any particular path to that cell
    # Score for a match, coming from the cell diagonally above
    diag_score = score_matrix[i - 1][j - 1] + similarity
                             
    # Coming from a gap: check to see if you're opening or extending it
    # Coming from the cell above:
    if path_matrix[i][j - 1] == 'd': # you're opening a gap
        up_score = score_matrix[i][j - 1] - gap_opening_penalty
    else: # you're extending a gap
        up_score = score_matrix[i][j - 1] - gap_extension_penalty
    # coming from the cell to the left
    if path_matrix[i - 1][j] == 'd': # you're opening a gap
        left_score = score_matrix[i - 1][j] - gap_opening_penalty
    else: # you're extending a gap
        left_score = score_matrix[i - 1][j] - gap_extension_penalty                           

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





def get_pair_scores(pairs, gap_opening_penalty, gap_extension_penalty, blosum_matrix):
    pair_hiscores = [0]*len(pairs)
    for index, pair in enumerate(pairs):
        a = pair[0]
        b = pair[1]
        hiscore = (0,0,0)
        
        path_matrix = create_path_matrix(a,b)
        score_matrix, hiscore = create_score_matrix(a, b, path_matrix, blosum_matrix, gap_opening_penalty, gap_extension_penalty)
        
        pair_hiscores[index] = hiscore[0]
        print('hiscore of index ' + str(index) + ' is ' + str(hiscore))
    return pair_hiscores

def get_pair_scores_adjusted(pairs, gap_opening_penalty, gap_extension_penalty, blosum_matrix):
    pair_hiscores = [0]*len(pairs)
    for index, pair in enumerate(pairs):
        a = pair[0]
        b = pair[1]
        hiscore = (0,0,0)
        
        path_matrix = create_path_matrix(a,b)
        score_matrix, hiscore = create_score_matrix(a, b, path_matrix, blosum_matrix, gap_opening_penalty, gap_extension_penalty)
        
        pair_hiscores[index] = hiscore[0] / min(len(a), len(b))
        print('hiscore of index ' + str(index) + ' is ' + str(hiscore))
    return pair_hiscores


def get_pair_scores_traced(pairs, gap_opening_penalty, gap_extension_penalty, blosum_matrix):
    pair_hiscores = [0]*len(pairs)
    for index, pair in enumerate(pairs):
        a = pair[0]
        b = pair[1]
        hiscore = (0,0,0)
        
        path_matrix = create_path_matrix(a,b)
        score_matrix, hiscore = create_score_matrix(a, b, path_matrix, blosum_matrix, gap_opening_penalty, gap_extension_penalty)
        
        pair_hiscores[index] = hiscore[0]
        print('hiscore of index ' + str(index) + ' is ' + str(hiscore))
    return pair_hiscores


def traceback(score_matrix, start_pos):
    '''Find the optimal path through the matrix.
    This function traces a path from the bottom-right to the top-left corner of
    the scoring matrix. Each move corresponds to a match, mismatch, or gap in one
    or both of the sequences being aligned. Moves are determined by the score of
    three adjacent squares: the upper square, the left square, and the diagonal
    upper-left square.
    WHAT EACH MOVE REPRESENTS
        diagonal: match/mismatch
        up:       gap in sequence 1
        left:     gap in sequence 2
    '''

    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y         = start_pos
    move         = next_move(score_matrix, x, y)
    while move != END:
        if move == DIAG:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append(seq2[y - 1])
            x -= 1
            y -= 1
        elif move == UP:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append('-')
            x -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[y - 1])
            y -= 1

        move = next_move(score_matrix, x, y)

    aligned_seq1.append(seq1[x - 1])
    aligned_seq2.append(seq1[y - 1])

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


def next_move(score_matrix, x, y):
    diag = score_matrix[x - 1][y - 1]
    up   = score_matrix[x - 1][y]
    left = score_matrix[x][y - 1]
    if diag >= up and diag >= left:     # Tie goes to the DIAG move.
        return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
    elif up > diag and up >= left:      # Tie goes to UP move.
        return 2 if up != 0 else 0      # UP move or end.
    elif left > diag and left > up:
        return 3 if left != 0 else 0    # LEFT move or end.
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')




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
            pospairs_toadd = get_pair_scores(pospairs, gap_opening_penalty, gap_extension_penalty, blosum_matrix)
            pospairs_hiscores_all.append((gap_opening_penalty, gap_extension_penalty, pospairs_toadd))
            
            # For each negative pair, compute the highest scoring alignment
            negpairs_toadd = get_pair_scores(negpairs, gap_opening_penalty, gap_extension_penalty, blosum_matrix)
            negpairs_hiscores_all.append((gap_opening_penalty, gap_extension_penalty, negpairs_toadd))
     
    return pospairs_hiscores_all, negpairs_hiscores_all
    


def get_ROC_data(pospairs_hiscores, negpairs_hiscores):
    ROC_data = np.zeros((len(np.arange(0,3,0.01)),2))
    score_thresholds = np.arange(0,3,0.01)
    for index in range(0,len(np.arange(0,3,0.01))):
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


def get_cutoff_matrix(pospairs_hiscores_all, cutoff):
    '''
    Determine a gap opening and extension penalty and score at which you have (x%) of your scores above the cutoff.
    '''
    cutoff_array = []
    
    for score_threshold in np.arange(0,300,1):
        for index in range(0,len(pospairs_hiscores_all)):
            above_threshold = 0
            for score in range(0,len(pospairs_hiscores_all[index][2])):
                if pospairs_hiscores_all[index][2][score] > score_threshold:
                    above_threshold += 1
            above_threshold = above_threshold / len(pospairs_hiscores_all[index][2])
            # print(str(above_threshold) + ' percent of the matches are above the score threshold')
    
            if above_threshold == cutoff:
                #print('0.7 occurs at a score_threshold of ' + str(score_threshold))   
                cutoff_array.append((pospairs_hiscores_all[index][0],pospairs_hiscores_all[index][1],score_threshold)) 
                #assert False
    return cutoff_array

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

'''
And now, we actually call some of the functions
'''




# Initialize your scoring matrix and path matrix
blosum_matrix = read_blosum('BLOSUM50')
blosum_matrix_50 = read_blosum('BLOSUM50')
blosum_matrix_62 = read_blosum('BLOSUM62')

# Read in all pairs of sequences to be alligned
pospairs = read_pairs('Pospairs.txt')
negpairs = read_pairs('Negpairs.txt')

#pospair_hiscores_50 = get_pair_scores_adjusted(pospairs, 9, 4, blosum_matrix_50)
#pospair_hiscores_62 = get_pair_scores_adjusted(pospairs, 9, 4, blosum_matrix_62)

#negpair_hiscores_50 = get_pair_scores_adjusted(negpairs, 9, 4, blosum_matrix_50)
#negpair_hiscores_62 = get_pair_scores_adjusted(negpairs, 9, 4, blosum_matrix_62)

roc_50 = get_ROC_data(pospair_hiscores_50,negpair_hiscores_50)
roc_62 = get_ROC_data(pospair_hiscores_62,negpair_hiscores_62)
    

