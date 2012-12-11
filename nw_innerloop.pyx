##############################################################################
# CSE527 Au12
# Problem Set 4 -- auxillary file
# James Youngquist
#
# The inner function from the Needleman-Wunsch algorithm was converted to
# Cython specific code following the tutorial at
#    http://wiki.cython.org/tutorials/numpy
#
##############################################################################

import numpy as np
cimport numpy as np
cimport cython

# Here we define the types that will be used throughout.
DTYPE = np.double
ctypedef np.double_t DTYPE_t

##############################################################################
# Create custom MAX functions that don't bother with all of the regular Python
# type checking etc.
cdef inline DTYPE_t float_max(DTYPE_t a, DTYPE_t b): return a if a >= b else b
cdef inline     int   int_max(    int a,     int b): return a if a >= b else b


##############################################################################
# The Needleman-Wunsch Inner Loop
@cython.boundscheck(False)      # turn off array bounds checking
def nwil(np.ndarray[DTYPE_t,  ndim=2] v_matrix  not None, # value matrix
         np.ndarray[np.uint8_t, ndim=1] seq1    not None, # sequence 1
         np.ndarray[np.uint8_t, ndim=1] seq2 not None, # sequence 2
         np.ndarray[DTYPE_t,  ndim=2] cost_mx   not None, # cost matrix
         gap_cost):                                       # gap cost

    # Just do a double check here
    assert v_matrix.dtype == DTYPE and cost_mx.dtype == DTYPE

    # Create type definitions so that cython can optimize out a lot of the
    # background type checking and conversion that python does by default.
    cdef DTYPE_t choice1, choice2, choice3
    cdef DTYPE_t cost = gap_cost
    cdef unsigned int row, col
    cdef unsigned int my_base, other_base
    cdef unsigned int rows = seq1.shape[0] + 1
    cdef unsigned int cols = seq2.shape[0] + 1

    # Fill in the value matrix.  This is done by working across the columns of
    # each row in turn, starting with row 1 and going down.  Accessing
    # elements in row-major order is slightly faster than column-major due to
    # the layout of the Numpy arrays in memory.
    row = 1
    while row != rows:
        my_base = seq1[row-1] # sequences are offset by 1 in the matrix
        col = 1
        while col != cols:
            other_base = seq2[col-1] # sequences are offset by 1 in the
                                          # matrix
            
            # - Moving along the diagonal corresponds to matching my_base
            #   with other_base. (1)
            # - Moving down a column corresponds to matching my_base with
            #   a dash after the optimal matching up till other_base. (2)
            # - Moving across a row corresponds to matching other_base with
            #   a dash after the optimal matching up till my_base. (3)
            choice1 = v_matrix[row-1,col-1] + cost_mx[my_base, other_base]
            choice2 = v_matrix[row-1,col  ] + cost
            choice3 = v_matrix[row,  col-1] + cost
            v_matrix[row,col] = float_max(choice1, float_max(choice2, choice3))
            col += 1
        row += 1
    return


##############################################################################
# Tracking back along the Needleman-Wunsch value matrix
@cython.boundscheck(False)      # turn off array bounds checking
def trackback(np.ndarray[DTYPE_t,  ndim=2] v_matrix  not None, # value matrix
              np.ndarray[np.uint8_t, ndim=1] seq1    not None, # sequence 1
              np.ndarray[np.uint8_t, ndim=1] seq2    not None, # sequence 2
              np.ndarray[DTYPE_t,  ndim=2] cost_mx   not None, # cost matrix
              gap_cost):                                       # gap cost
    cdef unsigned int row = seq1.shape[0]
    cdef unsigned int col = seq2.shape[0]
    cdef unsigned char dash_val = 8
    
    # These are arrays of chars with a default value of ' ' (space)
    cdef unsigned int result_idx = int_max(row+1, col+1)*2
    cdef np.ndarray my_result = np.zeros(result_idx, dtype=np.uint8)
    cdef np.ndarray other_result = np.zeros(result_idx, dtype=np.uint8)

    result_idx -= 1             # set to last element in the array
    
    while row != 0 or col != 0:
        if (row > 0 and col > 0 and
            v_matrix[row,col] - cost_mx[seq1[row-1], seq2[col-1]] == v_matrix[row-1,col-1]):
            row -= 1
            col -= 1
            my_result[result_idx] = seq1[row]
            other_result[result_idx] = seq2[col]
        elif row > 0 and v_matrix[row,col] - gap_cost == v_matrix[row-1, col]: # up
            row -= 1
            my_result[result_idx] = seq1[row]
            other_result[result_idx] = dash_val
        elif col > 0 and v_matrix[row,col] - gap_cost == v_matrix[row, col-1]: # left
            col -= 1
            my_result[result_idx] = dash_val
            other_result[result_idx] = seq2[col]
        else:
            raise LookupError("Don't know how you got here...")
        result_idx -= 1         # move to the previous element

    my_result_str = ''.join(map(chr, my_result)).strip()
    other_result_str = ''.join(map(chr, other_result)).strip()
    return my_result_str, other_result_str

              

