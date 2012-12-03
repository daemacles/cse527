import numpy as np

class Genome(object):
    '''
    This class represents a sequence of base pairs.  In order to deal with
    large genomes, we use memory mapping and demand loading for efficiency.
    '''

    def __init__(self, genomeFile):
        '''
        genomeFile - file-like object that represents a flat array of the
                     genome data.  It is assumed to be static and will be read
                     as needed to populate self.data
        '''
        pass
