import os
import sys
import pandas as pd
import re
import bisect
import gzip
import hashlib
from twobitreader import TwoBitFile
from string import maketrans
from optparse import OptionParser

def parse_bed_file (bed_file):

	with open(bed_file) as inf:

		for line in inf:
			ls=line.split()
			chrom,tstart,tend,tx,_,strand,cstart,cend,_,nexons,exon_size,exon_start=ls
			tstart,tend,cstart,cend,nexons=map(int,[tstart,tend,cstart,cend,nexons])
			if cend < cstart:
				cstart,cend=cend-1,cstart+1
			exon_size=map(int,exon_size.strip(',').split(','))
			exon_start=map(int,exon_start.strip(',').split(','))
			fe=bisect.bisect(exon_start,cstart-tstart)-1
			le=bisect.bisect(exon_start,cend-tstart)-1
			coding=(cstart< cend)
			if coding:
				utr5len=sum(exon_size[:fe])+(cstart-tstart-exon_start[fe])
				utr3len=exon_size[le]-(cend-tstart-exon_start[le])+sum(exon_size[le+1:])
				cdslen=sum(exon_size)-utr5len-utr3len
			else:
				utr5len=0
				utr3len=0
				cdslen=0
			if strand=='-':
				utr5len,utr3len=utr3len,utr5len
			yield (tx,dict(coding=coding,\
						   nexons=nexons,\
						   length=sum(exon_size),\
						   utr5len=utr5len,\
						   utr3len=utr3len,\
						   cdslen=cdslen,\
						   chrom=chrom,\
						   start=tstart,\
						   end=tend,\
						   strand=strand,\
						   exon_size=exon_size,\
						   exon_start=exon_start,\
						   cstart=cstart,\
						   cend=cend,\
						   rel_start=utr5len,\
						   rel_end=utr5len+cdslen))

codon_translate={'TTT': 'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

rc_tab=maketrans('ACGTUNacgtun','TGCAANtgcaan')

def RC (s):
	return s.translate(rc_tab)[::-1]

def translate (s):
	return ''.join(codon_translate[s[3*k:3*k+3]] for k in range(len(s)/3))

def get_sequence (genome, tx_info, rel_start, rel_end):

	tstart=tx_info['start']
	tend=tx_info['end']
	nexons=tx_info['nexons']
	exon_size=tx_info['exon_size']
	exon_start=tx_info['exon_start']
	strand=tx_info['strand']
	chrom=tx_info['chrom']
	txseq=''.join(genome[chrom][tstart+estart:tstart+estart+esize] for estart,esize in zip(exon_start,exon_size)).upper()
	if strand=='-':
		txseq=RC(txseq)
	return txseq[rel_start:rel_end]

def get_rel_coords (tx_info,cstart,cend):

	tstart=tx_info['start']
	tend=tx_info['end']
	nexons=tx_info['nexons']
	exon_size=tx_info['exon_size']
	exon_start=tx_info['exon_start']
	strand=tx_info['strand']
	tlen=tx_info['length']
	if cend < cstart:
		cstart,cend=cend-1,cstart+1
	fe=bisect.bisect(exon_start,cstart-tstart)-1
	le=bisect.bisect(exon_start,cend-tstart)-1
	rel_start=sum(exon_size[:fe])+(cstart-tstart-exon_start[fe])
	rel_end=sum(exon_size[:le])+(cend-tstart-exon_start[le])
	if strand=='-':
		rel_start,rel_end=tlen-rel_end,tlen-rel_start

	return (rel_start,rel_end)

def get_abs_coords (tx_info,rel_start,rel_end):

	tstart=tx_info['start']
	tend=tx_info['end']
	nexons=tx_info['nexons']
	exon_size=tx_info['exon_size']
	exon_start=tx_info['exon_start']
	strand=tx_info['strand']
	chrom=tx_info['chrom']
	tlen=tx_info['length']
	es=[0]+[sum(exon_size[:i+1]) for i in range(nexons)]
	if strand=='-':
		rel_start,rel_end=tlen-rel_end,tlen-rel_start
	fe=min(bisect.bisect(es,rel_start)-1,nexons-1)
	le=min(bisect.bisect_left(es,rel_end)-1,nexons-1)
	cstart=tstart+exon_start[fe]+rel_start-es[fe]+1
	cend=tstart+exon_start[le]+rel_end-es[le]

	return '{0}:{1}-{2}:{3}'.format(chrom,cstart,cend,strand)

def find_ORFs (tx, info, seq):

	stop_codons=['TAG', 'TAA', 'TGA']

	# get positions after stop codons in all 3 frames
	all_stops={0:[0],1:[1],2:[2]}
	for m in re.finditer('|'.join(stop_codons),seq):
		all_stops[m.start()%3].append(m.start()+3)
	# get relative coordinates of stop-stop in the same frame
	stop_stops=[(all_stops[f][n],all_stops[f][n+1]) for f in range(3) for n in range(len(all_stops[f])-1)]
	# get relative coordinates of start-stop in the same frame
	start_stops=[(stop1+seq[stop1:stop2].index('ATG'),stop2) for stop1,stop2 in stop_stops if 'ATG' in seq[stop1:stop2] and seq[stop1:stop2].index('ATG')%3==0]
	for start,stop in start_stops:
		yield (tx+'_'+get_abs_coords(info,start,stop),(start,stop,translate(seq[start:stop])))

parser=OptionParser()
parser.add_option('-b','--bed_file',dest='bed',help="bed-file to use (12-column UCSC style)")
parser.add_option('-G','--genome',dest='genome',help="genome (2bit file)")
parser.add_option('-s','--stats',dest='stats',help="write orf stats to this csv file")
parser.add_option('-o','--outfile',dest='outfile',help="write output to this file (default: stdout)")
parser.add_option('','--minlength',dest='minlength',default=6,help="""minimum ORF length (in nt, including stop) [12]""")

options,args=parser.parse_args()

print >> sys.stderr, 'reading genome from '+options.genome
genome=TwoBitFile(options.genome)
print >> sys.stderr, 'reading bed file '+options.bed

orfs=[]
hash_values=[]
orf_types=[]
orf_length=[]
utr5_length=[]
utr3_length=[]

def match_type(sig):
	if sig=='annotated':
		return 'annotated'
	elif sig=='utr5:cds:0':
		return "N-ext"
	elif 'utr5' in sig:
		return "5'UTR"
	elif 'utr3' in sig:
		return "3'UTR"
	elif 'cds' in sig:
		return "CDS"
	elif 'nc' in sig:
		return "noncoding"
	else:
		return 'other'

if options.outfile is None:
	print >> sys.stderr, 'writing to stdout'
	outf=sys.stdout
else:
	print >> sys.stderr, 'writing to '+options.outfile
	outf=open(options.outfile,'w')

for tx,info in parse_bed_file(options.bed):

	if info['coding']:
		rel_start,rel_end=get_rel_coords(info,info['cstart'],info['cend'])
		seq=get_sequence(genome,info,rel_start,rel_end)
		orf_hash=hashlib.md5(seq).hexdigest()
		orf=tx+'_{0}:{1}-{2}:{3}'.format(info['chrom'],info['cstart']+1,info['cend'],info['strand'])
		outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n'.format(info['chrom'],info['start'],info['end'],orf+':'+orf_hash+'_annotated',0,info['strand'],info['cstart'],info['cend'],0,info['nexons'],','.join(map(str,info['exon_size']))+',',','.join(map(str,info['exon_start']))+','))
		orfs.append(orf)
		hash_values.append(orf_hash)
		orf_types.append('annotated')
		orf_length.append(rel_end-rel_start)
		utr5_length.append(rel_start)
		utr3_length.append(info['length']-rel_end)

	for orf,(rel_start,rel_end,seq) in find_ORFs(tx, info, get_sequence(genome,info,0,info['length'])):
		if rel_end-rel_start <= options.minlength:
			continue
		rel_frame=(info['rel_start']-rel_start)%3

		if info['coding']:
			if rel_start < info['rel_start']:
				if rel_end < info['rel_start']:
					orf_type='utr5:utr5:{0}'.format(rel_frame)
				else:
					orf_type='utr5:cds:{0}'.format(rel_frame)
			elif rel_start < info['rel_end']:
				if rel_end <= info['rel_end']:
					orf_type='cds:cds:{0}'.format(rel_frame)
				else:
					orf_type='cds:utr3:{0}'.format(rel_frame)
			else:
				orf_type='utr3:utr3:{0}'.format(rel_frame)
		else:
			orf_type='nc:nc'

		coords=orf.split('_')[1]
		abs_start,abs_end=map(int,coords.split(':')[1].split('-'))
		if abs_start!=info['cstart'] and abs_end!=info['cend']:
			orf_hash=hashlib.md5(seq).hexdigest()
			outf.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n'.format(info['chrom'],info['start'],info['end'],orf+':'+orf_hash+'_'+orf_type,0,info['strand'],abs_start-1,abs_end,0,info['nexons'],','.join(map(str,info['exon_size']))+',',','.join(map(str,info['exon_start']))+','))
			orfs.append(orf)
			hash_values.append(orf_hash)
			orf_types.append(match_type(orf_type))
			orf_length.append(rel_end-rel_start)
			utr5_length.append(rel_start)
			utr3_length.append(info['length']-rel_end)

if options.stats is not None:
	print >> sys.stderr, 'saving stats as '+options.stats
	df=pd.DataFrame(data={'orf_type':orf_types,\
						  'orf_length':orf_length,\
						  'utr5_length':utr5_length,\
						  'utr3_length':utr3_length,\
						  'orf_hash': hash_values},index=orfs)
	df.to_csv(options.stats)

