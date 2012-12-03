##############################################################################
# File: inner_loops.pyx
# Desc: This file contains the cython optimized inner loops of various
#       functions.

import numpy as np
cimport numpy as np
cimport cython

BASETYPE = np.uint8             # type used for DNA bases
ctypedef np.uint8_t BASETYPE_t

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
    cdef np.ndarray result = np.zeros(num*4, dtype=np.uint8)
    while num > 0:
        next_byte = data[idx]
        idx += 1
        num -= 1
        bases = (base_lut[(next_byte & 0xC0) >> 6],
                 base_lut[(next_byte & 0x30) >> 4],
                 base_lut[(next_byte & 0x0C) >> 2],
                 base_lut[(next_byte & 0x03) >> 0])
        result[res_idx:res_idx+4] = bases
        res_idx += 4
    result = result[:length]
    return result
