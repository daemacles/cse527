##############################################################################
# This module supports the .2bit genome data format
# See http://genome.ucsc.edu/FAQ/FAQformat.html#format7

import numpy as np
import struct
import os

import pyximport; pyximport.install()
from inner_loops import twobit_decode

class BinaryReader(object):
    '''
    This is a wrapper around struct for easy reading of sequential values in a
    raw binary file.
    '''
    def __init__(self, data):
        self.data = data
        self.idx = 0
        self.e = '<'
        return

    def reset(self):
        self.idx = 0
        return

    def endian(self, e):
        if e == 'little':
            self.e = '<'
        else:
            self.e = '>'
        return

    def _parse(self, num, size, fmt_char):
        fmt = self.e + fmt_char*num
        ret = struct.unpack(fmt, self.data[self.idx:self.idx + (size*num)])
        self.idx += size*num
        if num == 1:
            return ret[0]
        else:
            return ret

    def uint32(self, num=1):
        return self._parse(num, 4, 'I')
    def uint8(self, num=1):
        return self._parse(num, 1, 'B')
    def char(self, num=1):
        ret = self.data[self.idx:self.idx+num]
        self.idx += num
        return ret

class Sequence(object):
    '''
    Abstract class that represents a sequence as an array, with an underlying
    memmapped numpy array backing it.
    '''
    def __getitem__(self):
        raise NotImplementedError()
    def __setitem__(self):
        raise NotImplementedError()

class TwoBitSeq(Sequence):

    #                    T C A G
    base_lut = np.array([0,1,2,3], dtype=np.uint8)
    #base_lut = np.array(map(ord, ['T', 'C', 'A', 'G']), dtype=np.uint8)
    
    def __init__(self, data, offset, name):
        self.name = name
        self.data = data        # 
        self._rec_start = offset    # where this record starts in data
        b_read = BinaryReader(data)
        b_read.idx = offset
        self.size = b_read.uint32() # number of bases
        self.nBlockCount = b_read.uint32() # number of N blocks
        self.nBlockStarts = [b_read.uint32() for i in xrange(self.nBlockCount)]
        self.nBlockSizes = [b_read.uint32() for i in xrange(self.nBlockCount)]
        self.maskBlockCount = b_read.uint32()
        self.maskBlockStarts = [b_read.uint32() for i in xrange(self.nBlockCount)]
        self.maskBlockSizes = [b_read.uint32() for i in xrange(self.nBlockCount)]
        reserved = b_read.uint32()

        # b_read.idx should now be pointing to the start of the packedDNA record
        self._seq_start = b_read.idx
        return

    def __getitem__(self, name):
        base_lut = self.base_lut
        data = self.data
        try:
            if not name.start: start = 0
            else: start = name.start
            result = twobit_decode(self.data, self.base_lut, self._seq_start,
                                 start, name.stop)
        except:
            next_byte = data[name/4 + self._seq_start]
            bases = (base_lut[(next_byte & 0xC0) >> 6],
                     base_lut[(next_byte & 0x30) >> 4],
                     base_lut[(next_byte & 0x0C) >> 2],
                     base_lut[(next_byte & 0x03) >> 0])
            result = np.array([bases[name % 4]])
        return result
        
class TwoBit(object):
    '''
    This class lets one treat a .2bit file as an array of flat arrays, one for
    each sequence.  The metadata included in the file and each sequence record
    are available as class members.  It decompresses on the fly packed gene
    data into FASTA.  Any caching is up to the user of this class.
    '''
    def __init__(self, filename):
        '''
        filename - string
        '''

        self.filename = filename
        parts = filename.split('/')
        basename = '/'.join(parts[:-1])
        self.name = parts[-1]
        self.type = '.2bit'
        data = np.memmap(filename, dtype=np.uint8, mode='r')
        self.data = data
        
        # Parse the .2bit header - 16 bytes = 4 32-bit fields
        b_read = BinaryReader(data)
        
        signature = b_read.uint32()
        if signature != 0x1A412743: # compare against magic number
            b_read.endian('big')
            b_read.reset()
            signature = b_read.uint32()
            if signature != 0x1A412743:
                raise ValueError(".2bit signature didn't match: %s"%hex(signature))
            
        version = b_read.uint32()
        if version != 0:
            raise ValueError('Version not 0')

        self.sequenceCount = b_read.uint32()
        reserved = b_read.uint32()

        # These (implicitly) maps sequence names to their offset in the file
        self.sequences = []
        self.sequence_map = {}

        idx = 16
        for seq_idx in range(self.sequenceCount):
            nameSize = b_read.uint8()
            name = ''.join(map(chr,b_read.char(nameSize)))
            offset = b_read.uint32()
            tbs = TwoBitSeq(data, offset, name)
            self.sequences.append((name, tbs))
            self.sequence_map[name] = tbs

        self.total_size = sum(s[1].size for s in self.sequences)
        return


##############################################################################
# Tests
if __name__ == "__main__":
    data = TwoBit('sacCer3_yeast2011.2bit')
    print data.sequenceCount
    s1 = data.sequences[0][1]
    print s1[0]
    print s1[:3]
    print s1[:20]
    print ''.join(map(chr,s1[:80]))
            
