import os
import sys
import numpy as np
import bisect
import pysam
from string import maketrans
from twobitreader import TwoBitFile
from optparse import OptionParser

rc_tab=maketrans('ACGTUNacgtun','TGCAANtgcaan')

def RC (s):
	return s.translate(rc_tab)[::-1]

codon_translate={'TTT': 'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

codons=sorted(codon_translate.keys())

parser=OptionParser()
parser.add_option('-b','--bed',dest='bed',help="12-column bed file with ORF definitions")
parser.add_option('-B','--bam',dest='bam',help="comma-separated list of BAM files with mapped reads, should have indices")
parser.add_option('-L','--L',dest='L',default='30',help="read lengths to use [30]")
parser.add_option('-e','--exclude',dest='exclude',default='25,1',help="codons to exclude at beginning and end [25,1]")
parser.add_option('-n','--names',dest='names',default=None,help="header names for bam files (comma-separated)")
parser.add_option('-o','--offsets',dest='offsets',default='15',help="offsets from 5' end of read [15]")
parser.add_option('-G','--genome',dest='genome',help="genome file (2bit format)")
parser.add_option('-s','--stranded',dest='stranded',default='yes',help="strand information (yes/no/reverse) [yes]")

options,args=parser.parse_args()

if options.genome is None:
	raise Exception("you need to provide a 2-bit genome file!")

genome=TwoBitFile(options.genome)

if options.bam is not None:
	print >> sys.stderr, 'using bam files',options.bam
	try:
		bam_files=[pysam.Samfile(bam.strip(),'rb') for bam in options.bam.split(',')]
		nmapped=np.array([bam.mapped for bam in bam_files])
		nB=len(bam_files)
	except:
		raise Exception ("couldn't open bam files")

	if options.names is not None:
		names=dict((n,x.strip()) for n,x in enumerate(options.names.split(',')))
		if len(names)!=nB:
			raise Exception("number of header names doesn't match number of bam files")
	else:
		names=dict(zip(range(nB),range(1,nB+1)))

	LL=[map(int,options.L.split('|')[n].split(',')) for n in range(nB)]
	Lmax=max(map(max,LL))
	Lmin=min(map(min,LL))
	offset=[dict(zip(LL[n],map(int,options.offsets.split('|')[n].split(',')))) for n in range(nB)]

else:
	bam_files=[]
	nB=0

exclude=map(int,options.exclude.split(','))

tx_old=''

print >> sys.stderr, 'using bed file',options.bed

sys.stdout.write('# bed files: '+options.bed)
if nB > 0:
	sys.stdout.write('\n# bam files:\n')
	for n in range(nB):
		sys.stdout.write('#  {0}: {1} ({2} reads)\n'.format(names[n],options.bam.split(',')[n],nmapped[n]))

sys.stdout.write('ORF')
if nB==0:
	sys.stdout.write('\t'+'\t'.join(codons)+'\n')
	
for n in range(nB):
	for c in codons:
		sys.stdout.write('\t{0}_{1}'.format(c,names[n]))
if nB > 1:
	for c in codons:
		sys.stdout.write('\t{0}_pooled'.format(c))
sys.stdout.write('\n')

nskipped=0
				 
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
		if orflen < 3*sum(exclude) or orflen%3!=0:
			nskipped+=1
			continue

		if tx!=tx_old:

			txseq=''.join(genome[chrom][tstart+estart:tstart+estart+esize] for estart,esize in zip(exon_start,exon_size)).upper()

			try:
				cov_site=np.zeros((nB,txlen),dtype=int)
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
					# libraries shouldn't be paired-end but use only one mate here just in case 
					if read.is_paired and read.is_read2:
						continue
					pos=np.array(read.positions)-tstart
					L=len(pos)
					if L not in offset[n]:
						continue
					if strand=='-' and (options.stranded=='yes' and read.is_reverse or\
										options.stranded=='reverse' and not read.is_reverse or\
										options.stranded=='no'):
						site=pos[-offset[n][L]-1]
					elif strand=='+' and (options.stranded=='yes' and not read.is_reverse or\
										  options.stranded=='reverse' and read.is_reverse or\
										  options.stranded=='no'):
						site=pos[offset[n][L]]
					else:
						continue
					if site >= 0 and site < txlen:
						cov_site[n,site]+=1

			cov_site_tx=np.concatenate([cov_site[:,est:est+esi] for est,esi in zip(exon_start,exon_size)],axis=1)
			tx_old=tx

		orf_seq=txseq[rel_start:rel_end]
		if strand=='-':
			orf_seq=RC(orf_seq)

		codons_here=np.array([orf_seq[3*k:3*k+3] for k in range(exclude[0],orflen/3-exclude[1])])
		aa_seq=''.join(codon_translate[c] for c in codons_here)

		if nB > 0:

			cov_site_orf=cov_site_tx[:,rel_start:rel_end]
			if strand=='-':
				cov_site_orf=cov_site_orf[:,::-1]

			# exclude ramp and stop codons as specified by exclude and sum up reads from 3 frames
			if exclude[1] > 0:
				cov_site_orf=np.sum([cov_site_orf[:,3*exclude[0]+k:-3*exclude[1]:3] for k in range(3)],axis=0)
			else:
				cov_site_orf=np.sum([cov_site_orf[:,3*exclude[0]+k::3] for k in range(3)],axis=0)

			if len(codons_here)!=cov_site_orf.shape[1]:
				raise Exception("lengths don't match!")

			if np.sum(cov_site_orf)==0 or '*' in aa_seq.rstrip('*'):
				nskipped+=1
				continue

			sys.stdout.write("{0}".format(name))

			codon_cov=np.array([np.sum(cov_site_orf[:,codons_here==c],axis=1) for c in codons])

			for n in range(nB):
				for k in range(len(codons)):
					sys.stdout.write('\t{0}'.format(codon_cov[k,n]))

			if nB > 1:

				# use pooled reads
				cov_site_orf=np.sum(cov_site_orf,axis=0)
				codon_cov=np.array([np.sum(cov_site_orf[codons_here==c]) for c in codons])

				for k in range(len(codons)):
					sys.stdout.write('\t{0}'.format(codon_cov[k]))
		else:

			sys.stdout.write("{0}".format(name))
			for c in codons:
				sys.stdout.write('\t{0}'.format(np.sum(codons_here==c)))



		sys.stdout.write("\n")

print >> sys.stderr, 'done ({0} skipped)'.format(nskipped)
