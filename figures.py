#!/usr/bin/env python

from cse527 import *

##############################################################################
# This script generates figures from the data I collected.  It is intended to
# be run from an interactivei python --pylab session.

# Interesting motifs, in order from lowest to highest hits
interesting = map(lambda s: np.array(map(numToBase.index, s), dtype=np.uint8),
                  ['CGGAAAAAAG', 'CTCCTTTCTT', 'ATAGAGGAAA', 'TTCATTATTG',
                   'CATATCATAA', 'ATTCATTGTC', 'TTTCTAATTC', 'TTCTTATGTT',
                   'ACAGAAACAA', 'AATTTTCACC', 'TGGGATTCCA', 'AACTTTGAAA',
                   'AAATGATGTA', 'TAATCGTAAT', 'AGTAAATGAA', 'GGAGAAAAAT',
                   'AAAAACATCA', 'ATTATCCTTT', 'CTTATTTTTG', 'TGTAAAATTA',
                   'AAACGTAAAA', 'TGTAAGAAAA', 'TTCAATTATT', 'GAAAAGCAAA',
                   'AAAACAATTT', 'ATAAAGGAAA', 'AAAAATTCTT', 'ATATACATTA',
                   'AAGAAAAATT', 'TACTTTTTTT'])

# Ordered list of chromasomes
chroms = ['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI']

genome = Genome('sacCer3_yeast2011.2bit')

try:
    RES
except NameError:
    RES = [gwConvolveTrack(genome, seq) for seq in interesting[-10:]]

def plotGenome(genome, chroms):
    idx = 0
    for chrom in chroms:
        X, Y = geneGraph(genome, chrom)
        X = np.array(X)
        Y = np.array(Y, dtype=np.float)
        Y = Y*.2 + (len(chroms) - idx - .5)
        idx += 1
        fill(X,Y,'red',linewidth=0.2)
        text(X[-1], Y[-1], ' '+chrom, horizontalalignment='left')
        title('Non-Gene Motifs on S. cerevisiae')
        xlabel('Base Pair Position on Chromosome')
        ylabel('Chromosome')

if __name__ == "__main__":
    # Create a new figure
    fig = figure()
    plotGenome(genome, chroms)

    # Draw lines for the top 10 sequences
    idx = 0
    symbols = '.*+ov^<>sp*h+xD1234|-'*4
    j = 0
    for X,Y,P in RES[-10:-1]:
        idx = 0
        for chrom in chroms:
            X = np.array(P[chrom])
            Y = np.ones(X.shape[0])*.2 + j*0.05 + (len(chroms) - idx - 0.5 + 0.1)
            if idx < len(chroms)-1:
                for i in range(X.shape[0]):
                    plot([X[i],X[i]],[Y[i],floor(Y[i]-.5)+.5], c='%f'%(j/10.0))
            else:
                for i in range(X.shape[0]):
                    plot([X[i],X[i]],[Y[i],floor(Y[i])+.5], c='%f'%(j/10.0))
            idx += 1
        #labels.append(map(seqToBaseString, interesting)[-j-2])
        j += 1
    idx = 0
    symbols = '.*+ov^<>sp*h+xD1234|-'*4

    # Draw symbols for the top 10 sequences
    j = 0
    for X,Y,P in RES[-10:-1]:
        idx = 0
        for chrom in chroms:
            X = np.array(P[chrom])
            Y = np.ones(X.shape[0])*.2 + j*0.05 + (len(chroms) - idx - 0.5 + 0.1)
            if idx < len(chroms)-1:
                plot(X,Y,symbols[j], c='%f'%(j/10.0))
            else:
                plot(X,Y,symbols[j], c='%f'%(j/10.0),
                     label='%2d) '%(10-j) + map(seqToBaseString, interesting[-10:-1])[j])
            idx += 1
        #labels.append(map(seqToBaseString, interesting)[-j-2])
        j += 1

    idx = 0
    for chrom in chroms:
        X,Y,P = RES[-1]
        X = np.array(P[chrom])
        Y = np.ones(X.shape[0])*.2 + (len(chroms) - idx - 0.5 + 0.05)
        if idx < len(chroms)-1:
            plot(X,Y,'ro')
        else:
            plot(X,Y,'ro', label=' 1) ' + seqToBaseString(interesting[-1]))
        idx += 1
    legend(prop={'family':'monospace', 'size':'medium'}, loc='lower right')
    print ' *** Please execute show() to display the chart ***'
