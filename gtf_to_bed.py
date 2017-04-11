import os
import sys
import gzip
import numpy as np
from optparse import OptionParser

parser=OptionParser()
parser.add_option('-i','--infile',dest='infile',help="input GTF (default: stdin)")
parser.add_option('-m','--merge',dest='merge',type=int,default=1,help="merge exons separated by n nts (default: n=1)")
parser.add_option('-c','--chrominfo',dest='chrominfo',help="ignore chromosomes not listed here, or transcripts exceeding chromosome dimensions")
parser.add_option('-a','--annotation',dest='annotation',help="also output annotation file with transcript/gene biotypes")

options,args=parser.parse_args()

if options.infile is None:
	inf=sys.stdin
else:
	if options.infile.endswith('.gz'):
		inf=gzip.open(options.infile)
	else:
		inf=open(options.infile)

if options.chrominfo is not None:
	chrominfo={}
	for line in open(options.chrominfo):
		ls=line.split()
		chrominfo[ls[0]]=int(ls[1])
	chroms=chrominfo.keys()

transcripts={}
for line in inf:
	if line.startswith('#'):
		continue
	ls=line.strip().split("\t")
	# this works for Ensembl v82 or Gencode gtf files (start codon is included in CDS so we can ignore it, stop codon is not)
	if ls[2] not in set(['exon','CDS','three_prime_utr','five_prime_utr','stop_codon','UTR']):
		continue
	chrom=ls[0]
	if 'chr' not in chrom:
		chrom='chr'+chrom
	start=int(ls[3])-1
	end=int(ls[4])
	score=0
	strand=ls[6]
	attr=dict((x.split()[0].strip(),x.split()[1].strip().strip('"')) for x in ls[8].strip(';').split(";"))
	if 'transcript_id' not in attr:
		print >> sys.stderr, "no transcript id in line"
		continue
	if options.chrominfo is not None:
		if chrom not in chrominfo:
			continue
		elif end >= chrominfo[chrom]:
			print >> sys.stderr, 'ignoring {0}: tx_end {1} > chromsize {2}'.format(attr['transcript_id'],end,chrominfo[chrom])
			continue
	tx=attr['transcript_id']
	if tx in transcripts:
		vals=transcripts[tx]
		if strand!=vals['strand'] or chrom!=vals['chrom']:
			print >> sys.stderr, "different strand or chrom for", tx
			continue
		vals['nexons']+=1
		vals['exon_start'].append(start)
		vals['exon_end'].append(end)
		vals['score'].append(score)
		if ls[2] in set(['CDS','stop_codon']):
			if 'cstart' in vals:
				vals['cstart']=min(start,vals['cstart'])
			else:
				vals['cstart']=start
			if 'cend' in vals:
				vals['cend']=max(end,vals['cend'])
			else:
				vals['cend']=end
	else:
		transcripts[tx]=dict(nexons=1,exon_start=[start],exon_end=[end],chrom=chrom,strand=strand,score=[score],attr=attr)

sorted_tx=sorted((tx for tx in transcripts), key=lambda tx: (transcripts[tx]['chrom'].strip('chr'),min(transcripts[tx]['exon_start']),min(transcripts[tx]['exon_end'])))

for tx in sorted_tx:
	vals=transcripts[tx]
	exon_start=[]
	exon_end=[]
	# merge adjoining exons from GTF file
	if options.merge is not None:
		for estart,eend in sorted(zip(vals['exon_start'],vals['exon_end']),key=lambda x: x[0]):
			if len(exon_end) > 0 and estart < exon_end[-1]+options.merge:
				exon_end[-1]=eend
			else:
				exon_start.append(estart)
				exon_end.append(eend)
	else:
		exon_start=vals['exon_start']
		exon_end=vals['exon_end']
	tstart=min(exon_start)
	tend=max(exon_end)
	rel_start=[es-tstart for es in exon_start]
	exon_size=[ee-es for ee,es in zip(exon_end,exon_start)]
	score=0
	cstart=vals['cstart'] if 'cstart' in vals else tend
	cend=vals['cend'] if 'cend' in vals else tend
	print '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}'.format(vals['chrom'],tstart,tend,tx,score,vals['strand'],cstart,cend,'0',len(exon_size),','.join(map(str,exon_size))+',',','.join(map(str,rel_start))+',')
	
	
if options.annotation is not None:
	with open(options.annotation,'w') as outf:
		fields=['Ensembl Gene ID','Ensembl Transcript ID','Gene Biotype','Transcript Biotype']
		transcript_attrs=set.union(*[set(tx['attr']) for tx in transcripts.values()])
		if 'gene_biotype' in transcript_attrs and 'transcript_biotype' in transcript_attrs:
			keys=['gene_id','transcript_id','gene_biotype','transcript_biotype']
		elif 'gene_type' in transcript_attrs and 'transcript_type' in transcript_attrs:
			keys=['gene_id','transcript_id','gene_type','transcript_type']
		else:
			raise Exception("can't find gene_type or transcript_type in attributes!")
		outf.write('\t'.join(keys)+'\n')
		for tx in sorted_tx:
			outf.write('\t'.join(transcripts[tx]['attr'][k] for k in keys)+'\n')
		
