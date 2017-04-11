import sys
import bisect
import pandas as pd
from optparse import OptionParser

def get_CDSlen (ls):

	chrom,tstart,tend,_,strand,cstart,cend,_,nexons,exon_size,exon_start=ls[range(11)]
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
		return sum(exon_size)-utr5len-utr3len
	else:
		return 0

parser=OptionParser()
parser.add_option('-b','--bed',dest='bed',help="""12-column bed file with ORF definitions""")
parser.add_option('-a','--annotation',dest='annotation',help="""annotation file with transcript-gene mappings (from gtf_to_bed.py)""")
parser.add_option('-o','--out',dest='out',help="""output file [stdout]""")

options,args=parser.parse_args()

annotation=pd.read_csv(options.annotation,header=None,index_col=3,sep='\t')
tx_stats=pd.read_csv(options.bed,header=0,index_col=1,sep='\t')
annotation=annotation.join(tx_stats)
annotation['CDSlen']=annotation.apply(get_CDSlen,axis=1)
annotation=annotation.groupby('gene_id').aggregate(lambda df: df.loc[df['CDSlen'].argmax()]).reset_index()
annotation[3]=annotation['gene_id']
annotation.to_csv(sys.stdout if options.out is None else options.out,sep='\t',header=None,index=None,columns=range(12))
