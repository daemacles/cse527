##############################################################################
# File: inner_loops.pyx
# Desc: This file contains the cython optimized inner loops of various
#       functions.

import numpy as np
cimport numpy as np
cimport cython

import util

import pyximport; pyximport.install()
from nw_innerloop import nwil, trackback

BASETYPE = np.uint8             # type used for DNA bases
ctypedef np.uint8_t BASETYPE_t

DOUBLE = np.double
ctypedef np.double_t DOUBLE_t

##############################################################################
# Create custom MAX functions that don't bother with all of the regular Python
# type checking etc.
cdef inline DOUBLE_t float_max(DOUBLE_t a, DOUBLE_t b): return a if a >= b else b
cdef inline int int_max(int a, int b): return a if a >= b else b

##############################################################################
# TwoBit decoding
@cython.boundscheck(False)
def twobit_decode(np.ndarray[BASETYPE_t, ndim=1] data not None,
                  np.ndarray[BASETYPE_t, ndim=1] base_lut not None,
                  seq_start,
                  slice_start,
                  slice_stop):
    cdef unsigned int length = slice_stop - slice_start
    cdef unsigned int byte_quotient = length / 4
    cdef unsigned int byte_remainder = length % 4 > 0
    cdef unsigned int num = byte_quotient + byte_remainder
    cdef unsigned int idx = (slice_start/4) + seq_start
    cdef unsigned int res_idx = 0
    cdef np.ndarray result = np.zeros(num*4, dtype=BASETYPE)
    while num > 0:
        next_byte = data[idx]
        idx += 1
        num -= 1
        result[res_idx  ] = (next_byte & 0xC0) >> 6
        result[res_idx+1] = (next_byte & 0x30) >> 4
        result[res_idx+2] = (next_byte & 0x0C) >> 2
        result[res_idx+3] = (next_byte & 0x03) >> 0
        # result[res_idx  ] = base_lut[(next_byte & 0xC0) >> 6]
        # result[res_idx+1] = base_lut[(next_byte & 0x30) >> 4]
        # result[res_idx+2] = base_lut[(next_byte & 0x0C) >> 2]
        # result[res_idx+3] = base_lut[(next_byte & 0x03) >> 0]
        res_idx += 4
    result = result[:length]
    return result

##############################################################################
# Returns the Reverse Complement of a sequence
@cython.boundscheck(False)
def seqRC(np.ndarray[BASETYPE_t, ndim=1] seq not None):
    cdef unsigned int length = seq.shape[0]
    cdef unsigned int idx = 0
    cdef np.ndarray out_seq = np.zeros(length, dtype=BASETYPE)
    while idx < length:
        out_seq[length - idx - 1] = (seq[idx] + 2) % 4
        idx += 1
    return out_seq
        
##############################################################################
# Increments the count of amino acid triples
@cython.boundscheck(False)
def aaTripleCount(np.ndarray[BASETYPE_t, ndim=1] seq not None,
                  np.ndarray[DOUBLE_t, ndim=3] counts not None):
    TTA = util.tripletToAmino
    cdef unsigned int num_bases = seq.shape[0] - 6
    if TTA[(seq[-3] << 4) + (seq[-2] << 2) + (seq[-1] << 0)] > 0:
        num_bases -= 3          # ignore STOP codon
    cdef unsigned int idx = 0
    cdef unsigned int i1, i2, i3
    while idx < num_bases:
        i1 = TTA[(seq[idx  ] << 4) + (seq[idx+1] << 2) + (seq[idx+2] << 0)]
        i2 = TTA[(seq[idx+3] << 4) + (seq[idx+4] << 2) + (seq[idx+5] << 0)]
        i3 = TTA[(seq[idx+6] << 4) + (seq[idx+7] << 2) + (seq[idx+8] << 0)]
        counts[i1,i2,i3] += 1.0
        idx += 3
    return

##############################################################################
# Does a lookup of 9 basepairs
@cython.boundscheck(False)
def aaLookup(np.ndarray[BASETYPE_t, ndim=1] seq not None,
             np.ndarray[DOUBLE_t, ndim=3] counts not None):
    pass

##############################################################################
# Convolves a (short) guess sequence against a (longer) target sequence
def convolveSubsequence(np.ndarray[BASETYPE_t, ndim=1] guess not None,
                        np.ndarray[BASETYPE_t, ndim=1] target_seq not None,
                        score_deque):
    cdef unsigned int slack = 1
    cdef unsigned int guess_len = guess.shape[0]
    cdef unsigned int target_len = target_seq.shape[0]
    cdef DOUBLE_t score_max
    # Account for really short target sequences
    cdef unsigned int num_tests = int_max(target_len - guess_len - slack + 1, 1)
    cdef np.ndarray v_matrix = util.initValueMatrix(guess_len, guess_len+slack, -1)
    cdef unsigned int i
    for i in range(num_tests):
        nwil(v_matrix, guess, target_seq, util.nw_basepairs_cost, -1,
             i, i+guess_len+slack)
        score = v_matrix[guess_len, guess_len+slack]
        if score > guess_len - 2 - slack:
            score_deque.append((i, score))
        score_max = float_max(score, score_max)
    return score_max
