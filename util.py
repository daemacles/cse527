import numpy as np

import pyximport; pyximport.install()
from nw_innerloop import nwil, trackback

DTYPE = np.double

# Maps numeric values to characters for bases and vice versa
numToBase = ('T', 'C', 'A', 'G')
baseToNum = np.zeros(128, np.uint8)
for b in numToBase:
    baseToNum[ord(b)] = numToBase.index(b)

START = np.array(map(numToBase.index, 'ATG'), dtype=np.uint8)
STOP_OCHRE = np.array(map(numToBase.index, 'TAA'), dtype=np.uint8)
STOP_AMBER = np.array(map(numToBase.index, 'TAG'), dtype=np.uint8)
STOP_OPAL = np.array(map(numToBase.index, 'AGA'), dtype=np.uint8)

# Maps numeric values to characters for amino acids and vice versa
numToAmino = ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
#              0    1    2    3    4    5    6    7    8    9
              'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 
#              10   11   12   13   14   15   16   17   18   19 
              '$', 'B', 'Z', 'X', '*')
#              20   21   22   23   24,   $ := Stop
aminoToNum = np.zeros(128, np.uint8)
for b in numToAmino:
    aminoToNum[ord(b)] = numToAmino.index(b)

# Create a LUT for BP triplet to amino acids
def tripletToIdx(t):
    return ((baseToNum[ord(t[0])] << 4) +
            (baseToNum[ord(t[1])] << 2) +
            (baseToNum[ord(t[2])] << 0))

lut_pairs = (                   # this is a human readable the 'master' list
    ('A', ('GCT', 'GCC', 'GCA', 'GCG')),
    ('R', ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG')),
    ('N', ('AAT', 'AAC')),
    ('D', ('GAT', 'GAC')),
    ('C', ('TGT', 'TGC')),
    ('Q', ('CAA', 'CAG')),
    ('E', ('GAA', 'GAG')),
    ('G', ('GGT', 'GGC', 'GGA', 'GGG')),
    ('H', ('CAT', 'CAC')),
    ('I', ('ATT', 'ATC', 'ATA')),
    ('L', ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG')),
    ('K', ('AAA', 'AAG')),
    ('M', ('ATG',)),
    ('F', ('TTT', 'TTC')),
    ('P', ('CCT', 'CCC', 'CCA', 'CCG')),
    ('S', ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC')),
    ('T', ('ACT', 'ACC', 'ACA', 'ACG')),
    ('W', ('TGG',)),
    ('Y', ('TAT', 'TAC')),
    ('V', ('GTT', 'GTC', 'GTA', 'GTG')),
    ('$', ('TAA', 'TGA', 'TAG'))) # $ = Stop

# Here we convert lut_paris to a LUT, where each triplet as a binary number
# specifies an index, and where the leftmost base defines the 2 most
# significant bits.
#
# For example tripletToIdx('GCT') = (3 << 4) + (1 << 2) + 0 = 52
#             tripletToAmino[52] = 0
#             numToAmino[0] = 'A'
# Thus, eventually, we can get that GCT encodes for Alanine
tripletToAmino = np.zeros(64, np.uint8) # the LUT
for amino, triplets in lut_pairs:
    for t in triplets:
        tripletToAmino[tripletToIdx(t)] = aminoToNum[ord(amino)]

##############################################################################
# Needleman-Wunsch cost matrix for base pairs
nw_basepairs_cost = np.matrix([[ 1, -1, -1, -1],
                               [-1,  1, -1, -1],
                               [-1, -1,  1, -1],
                               [-1, -1, -1,  1]], dtype=DTYPE)

##############################################################################
# Loads a cost matrix from a file and returns it as a look-up matrix indexed
# by the capitalized ascii values of fields.  Any uses of the returned matrix
# should first convert value to look up to upper case.
def loadCostFile(filename):
    lines = open(filename, 'r').readlines()
    matrix = np.zeros([25,25], dtype=DTYPE)
    indices = map(numToAmino.index, lines[0].strip().split())
    row_idx = 0
    for line in lines[1:]:
        vals = line.split()
        r = indices[row_idx]
        col_idx = 0
        for val in vals[1:]:
            c = indices[col_idx]
            matrix[r,c] = float(val)
            col_idx += 1
        row_idx += 1
    return matrix

# BLOSUM62 is our default cost matrix
blosum62 = loadCostFile('blosum.txt')


def seqToBaseString(seq):
    '''Returns a string represenation of the base pairs in a sequence.'''
    global numToBase
    return ''.join(map(lambda a: numToBase[a], seq))

def seqToAminoString(seq):
    '''Returns the amino acids encoded by a sequences.'''
    if len(seq) % 3:
        raise ValueError("The length of seq is not evenly divisible by 3.")
    global tripletToAmino
    aminos = [numToAmino[tripletToAmino[(seq[i] << 4) + (seq[i+1] << 2) + seq[i+2]]] for
              i in xrange(0, len(seq), 3)]
    string = ''.join(aminos)
    if string[-1] != '$':
        print "This sequence did not end in the stop codon"
        return string
    else:
        return string[:-1]


def print_genes(seq):
    in_gene = False
    meth_count = 0
    gene_count = 0
    start = 0
    genes = []
    i = 0
    while i < len(seq)-3:
        next3 = seq[i:i+3]
        '''
        if in_gene:
            if ((next3 == STOP_OCHRE).all() or 
                (next3 == STOP_AMBER).all() or
                (next3 == STOP_OPAL).all()):
                in_gene = False
                gene_count += 1
                length = i-start
                genes.append((start, i, length, length/3, length%3))
            if (next3 == START).all():
                meth_count += 1
            i += 3
        else:                   # not in a gene
            if (next3 == START).all():
                in_gene = True
                start = i
                i += 3
            else:
                i += 1
        '''
        if (next3 == START).all():
            meth_count += 1
            i += 2              # skip rest of this
        i += 1
    print 'Gene count', gene_count
    print 'Meth count', meth_count
    return genes


def align(seq1, seq2, gap_cost, cost_mx, v_matrix, BT=True):
    '''
    Calculates the optimal alignment of sequence SEQ2 against this
    sequence using the Needleman-Wunsch algorithm.

    If BT is False, then no backtracking is performed and just the score is
    given.
    '''

    # 0) Initialize the value matrix, if needed.
    #v_matrix = initValueMatrix(seq2)

    # 1) Fill in the value matrix.
    nwil(v_matrix, seq1, seq2, cost_mx, gap_cost)
    score = v_matrix[len(seq1), len(seq2)]

    if BT:
        # 2) Trace back along an optimal path. Set row and col to the bottom
        #    right index of v_matrix to start.
        my_result, seq2_result = trackback(v_matrix, seq1, seq2,
                                           cost_mx, gap_cost)

        return score, my_result, seq2_result
    else:
        return score, '', ''

def initValueMatrix(len_seq1, len_seq2, gap_cost):
    '''
    This simply initializes a value matrix for use in the align
    method. Also creates a track back matrix for back tracking an optimal
    alignment.

    The rows are the elements of this sequence and the columns are the
    elements of seq2.
    '''
    v_matrix = np.zeros([len_seq1 + 1, len_seq2 + 1],
                        dtype=DTYPE)
    v_matrix[0,:] = np.arange(len_seq2 + 1) * gap_cost
    v_matrix[:,0] = np.arange(len_seq1 + 1) * gap_cost
    return v_matrix
    
