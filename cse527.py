import numpy as np
import os
from filetypes import TwoBit
import sqlite3
import time
from collections import namedtuple

from util import *
import pyximport; pyximport.install()
from inner_loops import seqRC, aaTripleCount, convolveSubsequence

IterRecord = namedtuple('IterRecord', ['start', 'end', 'name', 'strand'])
                
class Genome(object):
    '''
    This class represents a sequence of base pairs.  In order to deal with
    large genomes, we use memory mapping and demand loading for efficiency.
    '''

    def __init__(self, dataFilename):
        '''
        seq_col - a collection of sequences
        '''

        self.name = ''.join(dataFilename.split('.')[:-1]) # strip file extension

        # Call appropriate data type handler
        seq_col = None
        if '.2bit' in dataFilename:
            seq_col = TwoBit('sacCer3_yeast2011.2bit')
            
        self.sequences = []

        # Create the backing cache file if it doesn't exist.
        self.mmap_file = self.name + '.cache'
        if not os.path.exists(self.mmap_file):
            self.data = np.memmap(self.mmap_file, dtype=np.uint8, mode='w+',
                                  shape=(seq_col.total_size,))
            # Now store all data in the cache
            print 'Creating', self.mmap_file, '...'
            print 'Please be patient, this could take a while.'
            start = 0
            end = 0
            idx = 1
            time_start = time.time()
            for name, sequence in seq_col.sequences:
                print '%3d/%3d Loading sequence %10s (%9d BPs)'%(idx, len(seq_col.sequences), name, sequence.size)
                idx += 1
                end += sequence.size
                self.data[start:end] = sequence[:sequence.size]
                start = end
            print 'Took %.2fs to cache %d sequences'%(time.time()-time_start,
                                                      len(seq_col.sequences))
        else:
            print 'Using existing cache:', self.mmap_file, '(delete to recreate it).'
            self.data = np.memmap(self.mmap_file, dtype=np.uint8, mode='r+',
                                  shape=(seq_col.total_size,))

        # Now create views onto our data
        start = 0
        end = 0
        for name, sequence in seq_col.sequences:
            end += sequence.size
            self.sequences.append((name, self.data[start:end]))
            start = end
        self.seq_map = dict(self.sequences)  # convenience attribute
        self.seq_names = self.seq_map.keys() # convenience attribute

        # Finally, try to make an annotations database
        self.conn = None
        self.loadAnnotations()
        return

    def loadAnnotations(self, geneFile=None, descFile=None):
        '''
        Creates a SQLite3 database if it doesn't exist with annotations about the
        genome of this organism

        It assumes the geneFile and descFile have tab separated fields exactly as
        those downloaded from (for example)
          http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/database/
        and renamed e.g. sgdGene.txt -> <self.name>.gene.csv
                         sgdDescription.txt -> <self.name>.desc.csv

        Returns if db/table creation was successful
        '''
        db_name = self.name + '.db'
        if os.path.exists(db_name):
            print 'Using existing db:', db_name, '(delete to recreate it).'
            self.conn = sqlite3.connect(db_name)
            return True         # assume if file exists, the database is fine
        
        if not geneFile: geneFile = self.name + '.gene.csv'
        if not descFile: descFile = self.name + '.desc.csv'
        if not os.path.exists(geneFile):
            print 'ERROR: could not find', geneFile,
            print 'please rerun loadAnnotations() directly with the right gene file name.'
            return False
        if not os.path.exists(descFile):
            print 'ERROR: could not find', descFile,
            print 'please rerun loadAnnotations() directly with the right description file name.'
            return False

        print 'Creating annotations database:', db_name
        # First set up the database
        conn = sqlite3.connect(db_name)
        self.conn = conn
        c = conn.cursor()
        c.execute('drop table if exists annotations')
        c.execute('''create table annotations (
                     name text primary key not null,
                     type text,
                     desc text,
                     seq_name text,
                     strand text,
                     txStart integer,
                     txEnd integer,
                     cdsStart integer,
                     cdsEnd integer,
                     exonCount integer,
                     exonStarts text,
                     exonEnds text,
                     proteinID text
                     )
                     ''')

        # Now read the gene and description files.
        # 1) Start with descriptions to create initial table rows
        for line in open(descFile):
            name, typ, desc = line.strip().split('\t')
            c.execute('''insert into annotations (name, type, desc) values (?, ?, ?)''',
                      (name, typ, desc))

        # 2) Now fill in the rest of the data
        for line in open(geneFile):
            (b, name, seq_name, strand, txStart, txEnd,
             cdsStart, cdsEnd, exonCount, exonStarts,
             exonEnds, proteinID) = line.strip().split('\t')
            c.execute('''update annotations set
    seq_name=?, strand=?, txStart=?, txEnd=?, cdsStart=?, cdsEnd=?,
    exonCount=?, exonStarts=?, exonEnds=?, proteinID=? where name=?
    ''',
                          (seq_name, strand, txStart, txEnd, cdsStart, cdsEnd,
                           exonCount, exonStarts, exonEnds, proteinID, name))

        # 3) And finally, drop all records described as 'dubious'
        c.execute('delete from annotations where desc like "%dubious%"')
        conn.commit()
        c.close()
        self.createNonGeneDB()
        return True

    def createNonGeneDB(self):
        '''
        This function creates a table in the DB of non-gene coding regions.
        Start and stop positions are relative to their sequence, not to the
        underlying data array.

        The column 'name' refers the gene immediately following this region.
        '''
        print 'Creating nongene table'
        old_factory = self.conn.row_factory
        self.conn.row_factory = sqlite3.Row
        c = self.conn.cursor()  # cursor for querying
        c2 = self.conn.cursor() # cursor for inserting
        c.execute('drop table if exists nongene')
        c.execute('''create table nongene (
                     name text primary key not null,
                     seq_name text,
                     strand text,
                     start integer,
                     end integer,
                     overlap integer
                     )
                     ''')
        cur_seq_name = None
        cur_name = None
        start = 0
        end = 0
        gene_end = 0
        overlaps = 0            # number of overlapping "genes"
        # Need to be careful here that overlapping "genes" are detected.
        gene_rows = c.execute('select * from annotations order by seq_name asc, txStart asc')
        for row in gene_rows:
            overlap = 0
            if cur_seq_name != row['seq_name']: # switched to a new sequence
                cur_seq_name = row['seq_name']
                start = 0
            end = row['cdsStart'] - 1
            if start > end:
                overlaps += 1
                overlap = 1
            else:
                c2.execute('insert into nongene (name, seq_name, strand, start, end, overlap) values (?, ?, "+", ?, ?, ?)',
                          (row['name'], cur_seq_name, start, end, overlap))
            start = row['cdsEnd'] + 1 # move next start to the first base after this gene
        self.conn.commit()
        c.close()
        c2.close()
        self.conn.row_factory = old_factory
        print "%d overlapping 'genes' detected"%overlaps
        return

    def _genericIter(self, cursor, seq_name):
        for row in cursor:
            start, end, name, strand = row
            if strand == '-':
                # If strand is reverse, return the reverse complement for the
                # coding strand
                seq = seqRC(self.seq_map[seq_name][start:end])
            else:
                seq = self.seq_map[seq_name][start:end]
            yield (seq, IterRecord._make(row))
        return

    def cdsIter(self, seq_name):
        '''
        Returns an iterator that yeilds tuples (sub_sequence, geneName)
        for the coding region of each gene on seq_name.
        '''
        try:
            c = self.conn.cursor()
        except AttributeError:
            raise AttributeError("Annotations database not yet created." +
                                 "Please run loadAnnotations().")
        query = 'select cdsStart, cdsEnd, name, strand from annotations where seq_name=? order by cdsStart asc'
        return self._genericIter(c.execute(query, (seq_name,)), seq_name)

    def txIter(self, seq_name):
        '''
        Returns an iterator that yeilds tuples (sub_sequence, geneName) for
        the transcription region of each gene on sequence (including UTRs).
        '''
        try:
            c = self.conn.cursor()
        except AttributeError:
            raise AttributeError("Annotations database not yet created." +
                                 "Please run loadAnnotations().")
        query = 'select txStart, txEnd, name, strand from annotations where seq_name=? order by txStart asc'
        return self._genericIter(c.execute(query, (seq_name,)), seq_name)
        
    def nongeneIter(self, seq_name):
        '''
        Returns an iterator that yeilds tuples (sub_sequence, following geneName)
        for the non-transcribed regions between genes on sequence.
        '''
        try:
            c = self.conn.cursor()
        except AttributeError:
            raise AttributeError("Database not yet created." +
                                 "Please run .loadAnnotations().")
        try:
            rows = c.execute('select start, end, name, strand from nongene where seq_name=? order by start asc',
                             (seq_name,))
        except sqlite3.OperationalError:
            raise AttributeError("Database doesn't have the nongene table. " +
                                 "Please run .createNonGeneDB()")
        return self._genericIter(rows, seq_name)

if __name__ == '__main__':
    import sys
    genome = Genome('sacCer3_yeast2011.2bit')
    chrI = genome.seq_map['chrI']

    # Calculate amino acid distribution
    '''
    try:
        gene_counts
    except NameError:
        gene_counts = np.zeros([21,21,21], dtype=np.double)
        for seq in genome.seq_names:
            print 'GC: Going over', seq
            for gene_seq, record in genome.txIter(seq):
                aaTripleCount(gene_seq, gene_counts)
        gene_counts /= sum(gene_counts.reshape(21**3))

    X = [0.0]
    Y = [0.0]
    for gene, record in genome.txIter('chrI'):
        X.extend((record.start, record.start, record.end, record.end))
        Y.extend((0.0, 1.0, 1.0, 0.0))
    '''

    '''
    TTA = tripletToAmino
    gX = np.arange(start, stop, step)

    chrI_gene_guess = np.zeros(chrI.shape[0]/3+1)
    seq = chrI
    start = 1
    stop = chrI.shape[0]-9
    step = 3
    for idx in xrange(start, stop, step):
        i1 = TTA[(seq[idx  ] << 4) + (seq[idx+1] << 2) + (seq[idx+2] << 0)]
        i2 = TTA[(seq[idx+3] << 4) + (seq[idx+4] << 2) + (seq[idx+5] << 0)]
        i3 = TTA[(seq[idx+6] << 4) + (seq[idx+7] << 2) + (seq[idx+8] << 0)]
        chrI_gene_guess[idx/3] = counts[i1,i2,i3]
    chrI_gene_guess /= max(chrI_gene_guess)
    '''

    
                                                           
                                                           
    
