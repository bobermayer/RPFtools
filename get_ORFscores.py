import os
import sys
import numpy as np
import bisect
import pysam
from optparse import OptionParser

parser=OptionParser()
parser.add_option('-b','--bed',dest='bed',help="12-column bed file with ORF definitions, should be sorted")
parser.add_option('-B','--bam',dest='bam',help="comma-separated list of BAM files with mapped reads, should have indices")
parser.add_option('-L','--L',dest='L',default='30',help="read lengths to use [30]")
parser.add_option('-o','--offsets',dest='offsets',default='12',help="offsets of P-sites from 5' end of read [12]")
parser.add_option('-s','--stranded',dest='stranded',default='yes',help="strand information (yes/no/reverse); default: yes")

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

tx_old=''

print >> sys.stderr, 'using bed file',options.bed

sys.stdout.write('# bed file: '+options.bed+'\n# bam files:\n')
for n in range(nB):
	sys.stdout.write('#  {0}: {1} ({2} reads)\n'.format(n+1,options.bam.split(',')[n],nmapped[n]))

sys.stdout.write('# ORF')
for n in range(nB):
	sys.stdout.write('\tRPM_tx_{0}\tRPM_orf_{0}\tORFscore_{0}\tcov_p0_{0}\tfrac_multimapper_{0}'.format(n+1))
if nB > 1:
	sys.stdout.write('\tRPM_tx_tot\tRPM_orf_tot\tORFscore_tot\tcov_p0_tot\tfrac_multimapper_tot\n')
else:
	sys.stdout.write('\n')
				 
with open(options.bed) as inf:

	for line in inf:

		ls=line.split()
		chrom,tstart,tend,name,_,strand,cstart,cend,_,nexons,exon_size,exon_start=ls
		tx='|'.join([chrom,tstart,tend,strand,nexons,exon_size,exon_start])
		tstart,tend,cstart,cend,nexons=map(int,[tstart,tend,cstart,cend,nexons])
		exon_size=map(int,exon_size.strip(',').split(','))
		exon_start=map(int,exon_start.strip(',').split(','))
		txlen=tend-tstart
		fe=bisect.bisect(exon_start,cstart-tstart)-1
		le=bisect.bisect(exon_start,cend-tstart)-1
		rel_start=sum(exon_size[:fe])+(cstart-tstart-exon_start[fe])
		rel_end=sum(exon_size[:le])+(cend-tstart-exon_start[le])

		orflen=rel_end-rel_start

		if tx!=tx_old:

			try:
				cov_all=np.zeros((nB,txlen))
				cov_psite=np.zeros((nB,txlen),dtype=int)
				cov_multimapper=np.zeros((nB,txlen),dtype=int)
			except MemoryError:
				print >> sys.stderr, 'not enough memory for {0}; skipping this transcript'.format(name)
				continue

			for n,bam in enumerate(bam_files):

				if chrom not in bam.references:
					chrom=chrom.strip('chr')

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
					cov_all[n,pos[(pos >= 0) & (pos < txlen)]]+=1./L
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
						cov_psite[n,psite]+=1
						if is_multimapper:
							cov_multimapper[n,psite]+=1

			cov_all_tx=np.concatenate([cov_all[:,est:est+esi] for est,esi in zip(exon_start,exon_size)],axis=1)
			cov_psite_tx=np.concatenate([cov_psite[:,est:est+esi] for est,esi in zip(exon_start,exon_size)],axis=1)
			cov_multimapper_tx=np.concatenate([cov_multimapper[:,est:est+esi] for est,esi in zip(exon_start,exon_size)],axis=1)
			tx_old=tx
			cov=0

		cov_all_orf=cov_all_tx[:,rel_start:rel_end]
		cov_psite_orf=cov_psite_tx[:,rel_start:rel_end]
		cov_multimapper_orf=cov_multimapper_tx[:,rel_start:rel_end]
		if strand=='-':
			cov_all_orf=cov_all_orf[:,::-1]
			cov_psite_orf=cov_psite_orf[:,::-1]
			cov_multimapper_orf=cov_multimapper_orf[:,::-1]

		sys.stdout.write("{0}".format(name))

		for n in range(nB):

			RPM_tx=np.sum(cov_all_tx[n])/nmapped[n]*1.e6
			RPM_orf=np.sum(cov_all_orf[n])/nmapped[n]*1.e6

			if orflen==0 or orflen%3!=0:

				cov_p0=0
				score=0
				frac_multimapper=0

			else:

				# mask most abundant position if it comprises more than 70% of reads
				# ignore start codon and last 9nt near stop

				mask=np.zeros(orflen)
				m=np.argmax(cov_psite_orf[n])
				if cov_psite_orf[n][m] > .7*np.sum(cov_psite_orf[n]):
					mask[m]=cov_psite_orf[n][m]
				f=[np.sum((cov_psite_orf[n]-mask)[3+k:-9:3]) for k in range(3)]
				fbar=np.sum(f)/3.
				if fbar==0:
					cov_p0=0
					score=0
				else:
					cov_p0=np.sum((cov_psite_orf[n]-mask)[::3] > 0)/(orflen/3.)
					score=np.log2(np.sum((f-fbar)**2)/fbar+1.)
					if f[0] < f[1] or f[0] < f[2]:
						score=-score

				if np.sum(cov_psite_orf[n]) > 0:
					frac_multimapper=np.sum(cov_multimapper_orf[n])/float(np.sum(cov_psite_orf[n]))
				else:
					frac_multimapper=0

			sys.stdout.write('\t{0:.4e}\t{1:.4e}\t{2:.4e}\t{3:.4e}\t{4:.2f}'.format(RPM_tx,RPM_orf,score,cov_p0,frac_multimapper))

		if nB > 1:

			RPM_tx=np.sum(cov_all_tx)/np.sum(nmapped)*1.e6
			RPM_orf=np.sum(cov_all_orf)/np.sum(nmapped)*1.e6

			if orflen==0 or orflen%3!=0:

				cov_p0=0
				score=0
				frac_multimapper=0

			else:

				# pool reads from different files
				cov_psite_orf=np.sum(cov_psite_orf,axis=0)
				mask=np.zeros(orflen)
				m=np.argmax(cov_psite_orf)
				if cov_psite_orf[m] > .7*np.sum(cov_psite_orf):
					mask[m]=cov_psite_orf[m]
				f=np.array([np.sum((cov_psite_orf-mask)[3+k:-9:3]) for k in range(3)])
				fbar=np.sum(f)/3.
				if fbar==0:
					cov_p0=0
					score=0
				else:
					cov_p0=np.sum((cov_psite_orf-mask)[::3] > 0)/(orflen/3.)
					score=np.log2(np.sum((f-fbar)**2)/fbar+1.)
					if f[0] < f[1] or f[0] < f[2]:
						score=-score

				if np.sum(cov_psite_orf) > 0:
					frac_multimapper=np.sum(cov_multimapper_orf)/float(np.sum(cov_psite_orf))
				else:
					frac_multimapper=0

			sys.stdout.write('\t{0:.4e}\t{1:.4e}\t{2:.4e}\t{3:.4e}\t{4:.2f}'.format(RPM_tx,RPM_orf,score,cov_p0,frac_multimapper))

		sys.stdout.write("\n")

print >> sys.stderr, 'done'
