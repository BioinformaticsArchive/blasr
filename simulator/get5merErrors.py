#!/usr/bin/python
import sys

from pbpy.model.SeqUtils import *
from pbpy.io.cmph5 import factory

h5 = factory.create( sys.argv[1], 'r' )

K = 5
K2 = (K-1)/2

def update_counts( aligned_query, aligned_target, counts ):
    i = 0
    bases = []
    for ch1,ch2 in zip(aligned_query,aligned_target):
        if ch2=='-':
            try:
                bases[i-1][1].append( ch1 )
            except IndexError:
                return counts
        else:
            if ch1!=ch2:
                bases.append( (ch2,[ch1.lower()]) )
            else:
                bases.append( (ch2,[ch1]) )
            i += 1

    def is_match(seq):
        if len(seq)!=1: return False
        if seq[0]=='-': return False
        if seq[0].islower(): return False
        return True

    for i in xrange(len(bases)-K):
        seq = [ v[1] for v in bases[ i:i+K ] ]
        if not( all([is_match(s) for s in seq[0:K2]]) and \
            all([is_match(s) for s in seq[K2+1:]]) ):
            continue
        ref = ''.join([v[0] for v in bases[i:i+K]])
        obs = ''.join(seq[K2]).upper()
        if len(obs)>1:
            obs = obs.replace('-','')
        j = kmerToIndex(ref)
        if obs not in counts[j]: counts[j][obs]=0
        counts[j][obs] += 1
    return counts

def dump_counts( counts ):
    for k,v in enumerate(counts):
        kmer = indexToKmer( k, K )
        sv = ';'.join(['%s=%d'%(x,y) for x,y in v.iteritems()])
        print '%s\t%s' % ( kmer, sv )

counts = [ {} for i in xrange(4**K) ]

for i, hit in enumerate(h5.alnHitIterator()):
    aq = hit.alignedQuery.upper()
    at = hit.alignedTarget.upper()
    #if hit.zScore<6:
    #    continue
    counts = update_counts(aq,at,counts)
    #if i>100:
    #    break

dump_counts( counts )
