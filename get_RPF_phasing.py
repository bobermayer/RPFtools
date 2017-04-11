import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import bisect
import pysam
from matplotlib import pyplot as plt
from optparse import OptionParser

matplotlib.rcParams['font.size']=6
matplotlib.rcParams['axes.linewidth']=.5

parser=OptionParser()
parser.add_option('-b','--bed',dest='bed',help='12-column bed file with ORF definitions')
parser.add_option('-B','--bam',dest='bam',help='comma-separated list of BAM files with mapped reads, need to have indices')
parser.add_option('-l','--Lmin',dest='Lmin',default=27,type=int,help='min read length to use [27]')
parser.add_option('-L','--Lmax',dest='Lmax',default=31,type=int,help='max read length ot use [31]')
parser.add_option('-o','--out',dest='out',help='name of output figure file')
parser.add_option('-s','--stranded',dest='stranded',default='yes',help="strand information (yes/no/reverse) [yes]")
parser.add_option('','--offset',dest='offset',default=12,type=int,help="P-site offset [12]")

options,args=parser.parse_args()

Lmin=options.Lmin
Lmax=options.Lmax
nL=Lmax-Lmin+1

print >> sys.stderr, 'using bam files',options.bam
bam_files=[pysam.Samfile(bam.strip(),'rb') for bam in options.bam.split(',')]

print >> sys.stderr, 'using bed file',options.bed

cov_start=np.zeros((nL,108),dtype=int)
cov_stop=np.zeros((nL,108),dtype=int)
cov_tot=np.zeros((nL,3))
n=0
nskipped=0

with open(options.bed) as inf:

	for line in inf:

		ls=line.strip('\n').split('\t')
		chrom,tstart,tend,name,_,strand,cstart,cend,_,nexons,exon_size,exon_start=ls
		tx='|'.join([chrom,tstart,tend,strand,nexons,exon_size,exon_start])
		tstart,tend,cstart,cend,nexons=map(int,[tstart,tend,cstart,cend,nexons])
		tstart-=9
		tend+=9
		exon_size=map(int,exon_size.strip(',').split(','))
		exon_start=map(int,exon_start.strip(',').split(','))
		txlen=tend-tstart
		if cstart==cend:
			continue
		fe=bisect.bisect(exon_start,cstart-tstart)-1
		le=bisect.bisect(exon_start,cend-tstart)-1
		rel_start=sum(exon_size[:fe])+(cstart-tstart-exon_start[fe])
		rel_end=sum(exon_size[:le])+(cend-tstart-exon_start[le])

		if (rel_end-rel_start)%3!=0:
			nskipped+=1
			continue

		try:
			cov=np.zeros((nL,txlen),dtype=np.uint16)
		except MemoryError:
			print >> sys.stderr, 'not enough memory for {0}; skipping this transcript'.format(name)
			continue

		for bam in bam_files:
			if chrom not in bam.references:
				chrom=chrom.strip('chr')
			if chrom not in bam.references:
				continue
			for read in bam.fetch(chrom,max(0,tstart-Lmax),tend+Lmax):
				if read.is_unmapped or read.is_duplicate or read.is_qcfail:
					continue
				# ignore reads with soft-clipped nucleotides
				if read.cigar[0][0]==4 or read.cigar[-1][0]==4:
					continue
				pos=read.positions
				rlen=len(pos)
				if rlen < Lmin or rlen > Lmax:
					continue
				if strand=='-' and (options.stranded=='yes' and read.is_reverse or\
									options.stranded=='reverse' and not read.is_reverse or\
									options.stranded=='no'):
					psite=pos[-options.offset-1]-tstart
				elif strand=='+' and (options.stranded=='yes' and not read.is_reverse or\
									  options.stranded=='reverse' and read.is_reverse or\
									  options.stranded=='no'):
					psite=pos[options.offset]-tstart
				else:
					continue
				if psite >= 0 and psite < txlen:
					cov[rlen-Lmin,psite]+=1

		cov_here=np.concatenate([cov[:,est:est+esi] for est,esi in zip(exon_start,exon_size)],axis=1)

		cov_orf=cov_here[:,rel_start-9:rel_end+9]
		if strand=='-':
			cov_orf=cov_orf[:,::-1]

		if cov_orf.shape[1] >= 108:
			print >> sys.stderr, '\r{1:20s}\t({0} ORFs)'.format(n,name),
			cov_start+=cov_orf[:,:108]
			cov_stop+=cov_orf[:,-108:]
			cov_tot+=np.array([np.sum(cov_orf[:,(12+k):(-12+k):3],axis=1) for k in range(3)]).T
			n+=1

print >> sys.stderr, 'done ({0} skipped)'.format(nskipped)

fig=plt.figure(1,figsize=(7,nL))
fig.clf()

dL=1./nL-.02
d=dL-.06
top=(nL-1)/float(nL)+.01

for k in range(nL):

	fig.text(.01,top-(k-.5)*dL,'{0} nt'.format(range(Lmin,Lmax+1)[k]),size=6,rotation=90,ha='center',va='center')

	if np.sum(cov_start[k]) > 0:
		ax=fig.add_axes([.08,top-k*dL,.35,d])
		[ax.vlines(np.arange(-9,99)[i::3],0,cov_start[k,i::3],color='rgb'[i]) for i in range(3)]
		ax.set_xlim([-10,100])
	if k==nL-1:
		ax.set_xlabel('distance from start [nt]')

	if np.sum(cov_stop[k]) > 0:
		ax=fig.add_axes([.51,top-k*dL,.35,d])
		[ax.vlines(np.arange(-99,9)[i::3],0,cov_stop[k,i::3],color='rgb'[i]) for i in range(3)]
		ax.set_xlim([-100,10])
	if k==nL-1:
		ax.set_xlabel('distance from stop [nt]')

	if np.sum(cov_tot[k]) > 0:
		ax=fig.add_axes([.9,top-k*dL,.08,d])
		[ax.bar(i,cov_tot[k,i]/float(np.sum(cov_tot[k])),color='rgb'[i]) for i in range(3)]
		ax.set_ylim([0,1])
		ax.set_xticks([.5,1.5,2.5])
		ax.set_xticklabels([0,1,2])
	if k==nL-1:
		ax.set_xlabel('% frame')

fig.suptitle(options.bam)
fig.savefig(options.out)
