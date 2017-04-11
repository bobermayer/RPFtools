import os
import sys
import numpy as np
import bisect
import pysam
from optparse import OptionParser

parser=OptionParser()
parser.add_option('-b','--bed',dest='bed',help="12-column bed file with ORF definitions, should be sorted (with ORFs on the same transcript in subsequent entries)")
parser.add_option('-B','--bam',dest='bam',help="comma-separated list of BAM files with mapped reads, should have indices")
parser.add_option('-n','--names',dest='names',default=None,help="header names for bam files (comma-separated)")
parser.add_option('-L','--L',dest='L',default='30',help="read lengths to use (comma-separated; separated by | for each bam file) [30]")
parser.add_option('-o','--offsets',dest='offsets',default='12',help="offsets of P-sites from 5' end of read (corresponding to read lengths) [12]")
parser.add_option('-s','--stranded',dest='stranded',default='yes',help="strand information (yes/no/reverse) [yes]")
parser.add_option('-N','--N',dest='N',default=150,type=int,help="nucleotide windows around start and stop codon [150]")

options,args=parser.parse_args()

print >> sys.stderr, 'using bam files',options.bam
try:
	bam_files=[pysam.Samfile(bam.strip(),'rb') for bam in options.bam.split(',')]
	nmapped=np.array([bam.mapped for bam in bam_files])
	nB=len(bam_files)
except:
	raise Exception ("couldn't open bam files")

LL=[map(int,options.L.split('|')[n].split(',')) for n in range(nB)]
Lmax=max(map(max,LL))
Lmin=min(map(min,LL))
offset=[dict(zip(LL[n],map(int,options.offsets.split('|')[n].split(',')))) for n in range(nB)]

if options.names is not None:
	names=dict((n,x.strip()) for n,x in enumerate(options.names.split(',')))
	if len(names)!=nB:
		raise Exception("number of header names doesn't match number of bam files")
else:
	names=dict(zip(range(nB),range(1,nB+1)))

print >> sys.stderr, 'using bed file',options.bed

sys.stdout.write('# bed file: '+options.bed+'\n# bam files:\n')
for n in range(nB):
	sys.stdout.write('#  {0}: {1} ({2} reads)\n'.format(names[n],options.bam.split(',')[n],nmapped[n]))

sys.stdout.write('# gene')
for n in range(nB):
	sys.stdout.write('\tcov_start_{0}\tcov_stop_{0}'.format(names[n]))
sys.stdout.write('\n')
				 
with open(options.bed) as inf:

	for line in inf:

		ls=line.split()
		chrom,tstart,tend,name,_,strand,cstart,cend,_,nexons,exon_size,exon_start=ls
		tstart,tend,cstart,cend,nexons=map(int,[tstart,tend,cstart,cend,nexons])
		txlen=tend-tstart
		exon_size=map(int,exon_size.strip(',').split(','))
		exon_start=map(int,exon_start.strip(',').split(','))
		fe=bisect.bisect(exon_start,cstart-tstart)-1
		le=bisect.bisect(exon_start,cend-tstart)-1
		rel_start=sum(exon_size[:fe])+(cstart-tstart-exon_start[fe])
		rel_end=sum(exon_size[:le])+(cend-tstart-exon_start[le])

		try:
			cov=np.zeros((nB,txlen),dtype=int)
		except MemoryError:
			print >> sys.stderr, 'not enough memory for {0}; skipping this transcript'.format(name)
			continue

		for n,bam in enumerate(bam_files):

			if chrom not in bam.references:
				chrom=chrom.strip('chr')
			if chrom not in bam.references:
				continue

			for read in bam.fetch(chrom,max(0,tstart-Lmax),tend+Lmax):
				if read.is_unmapped or read.is_duplicate or read.is_qcfail:
					continue
				try:
					is_multimapper=sum(k > 1 for t,k in read.tags if t=='NH')
				except OverflowError:
					is_multimapper=False
					pass
				pos=np.array(read.positions)-tstart
				L=len(pos)
				if L not in offset[n]:
					continue
				if strand=='-' and (options.stranded=='yes' and read.is_reverse or\
									options.stranded=='reverse' and not read.is_reverse or\
									options.stranded=='no'):
					psite=pos[-offset[n][L]-1]
				elif strand=='+' and (options.stranded=='yes' and not read.is_reverse or\
									  options.stranded=='reverse' and read.is_reverse or\
									  options.stranded=='no'):
					psite=pos[offset[n][L]]
				else:
					continue
				if psite >= 0 and psite < txlen:
					cov[n,psite]+=1

		cov=np.concatenate([cov[:,est:est+esi] for est,esi in zip(exon_start,exon_size)],axis=1)
		txlen=cov.shape[1]

		# normalize to coverage in enterior of CDS
		# mean_cov=np.mean(cov[:,rel_start+options.N/2:rel_end-options.N/2],axis=1)
		# normalize to RPKM of ribo reads over transcript
		# mean_cov=np.sum(cov,axis=1)/(nmapped/1.e6)/(txlen/1.e3)
		# normalize to mean number of reads per transcript
		mean_cov=np.mean(cov,axis=1)
		mean_cov[mean_cov==0]=np.nan

		if np.sum(np.isfinite(mean_cov))==0:
			continue

		# get 2*N nt normalized coverage around start & stop, fill with NaN in intergenic space
		cov_start=np.concatenate([np.nan*np.ones((nB,max(0,options.N-rel_start))),\
								  cov[:,max(rel_start-options.N,0):min(txlen,rel_start+options.N)],\
								  np.nan*np.ones((nB,max(0,rel_start+options.N-txlen)))],axis=1)/mean_cov[:,None]
		cov_stop=np.concatenate([np.nan*np.ones((nB,max(0,options.N-rel_end))),\
								 cov[:,max(rel_end-options.N,0):min(txlen,rel_end+options.N)],\
								 np.nan*np.ones((nB,max(0,rel_end+options.N-txlen)))],axis=1)/mean_cov[:,None]

		sys.stdout.write("{0}".format(name))

		for n in range(nB):

			if strand=='+':
				sys.stdout.write('\t{0}\t{1}'.format(','.join('{0:.4e}'.format(c) for c in cov_start[n]),
													 ','.join('{0:.4e}'.format(c) for c in cov_stop[n])))
			else:
				sys.stdout.write('\t{0}\t{1}'.format(','.join('{0:.4e}'.format(c) for c in cov_stop[n][::-1]),
													 ','.join('{0:.4e}'.format(c) for c in cov_start[n][::-1])))

		sys.stdout.write("\n")

print >> sys.stderr, 'done'
