import sys
import bisect
import pandas as pd
from string import maketrans
from optparse import OptionParser
from twobitreader import TwoBitFile

rc_tab=maketrans('ACGTUNacgtun','TGCAANtgcaan')

def RC (s):
    return s.translate(rc_tab)[::-1]

def get_CDSlen (ls, only_complete=False, genome=None):

    chrom,tstart,tend,_,strand,cstart,cend,_,nexons,exon_size,exon_start=ls[range(11)]
    tstart,tend,cstart,cend,nexons=map(int,[tstart,tend,cstart,cend,nexons])
    if cend < cstart:
        cstart,cend=cend-1,cstart+1
    exon_size=map(int,exon_size.strip(',').split(','))
    exon_start=map(int,exon_start.strip(',').split(','))
    fe=bisect.bisect(exon_start,cstart-tstart)-1
    le=bisect.bisect(exon_start,cend-tstart)-1
    rel_start=sum(exon_size[:fe])+(cstart-tstart-exon_start[fe])
    rel_end=sum(exon_size[:le])+(cend-tstart-exon_start[le])
    coding=(cstart < cend)
    if not coding:
        return None
    if only_complete and chrom in genome:
        seq=''.join(genome[chrom][tstart+estart:tstart+estart+esize].upper() for estart,esize in zip(exon_start,exon_size))[rel_start:rel_end]
        if strand=='-':
            seq=RC(seq)
        if len(seq)%3==0 and seq[:3]=='ATG' and seq[-3:] in ['TAA','TGA','TAG']:
            return rel_end-rel_start
        else:
            return None
    else:
        return rel_end-rel_start

parser=OptionParser()
parser.add_option('-b','--bed',dest='bed',help="""12-column bed file with ORF definitions""")
parser.add_option('-a','--annotation',dest='annotation',help="""annotation file with transcript-gene mappings (from gtf_to_bed.py)""")
parser.add_option('-o','--out',dest='out',help="""output file [stdout]""")
parser.add_option('-G','--genome',dest='genome',help="""genome file (2bit format)""")
parser.add_option('','--only_complete',dest='only_complete',action='store_true',default=False,help="use only complete ORFs (ATG-stop) ")

options,args=parser.parse_args()

# read in bed file and combine with transcript-gene mappings from options.annotation
annotation=pd.read_csv(options.bed,header=None,index_col=3,sep='\t').join(pd.read_csv(options.annotation,header=0,index_col=1,sep='\t'))

# compute CDS length
if options.only_complete:
    from twobitreader import TwoBitFile
    genome=TwoBitFile(options.genome)
    annotation['CDSlen']=annotation.apply(lambda df: get_CDSlen(df, only_complete=True, genome=genome),axis=1)
else:
    annotation['CDSlen']=annotation.apply(get_CDSlen,axis=1)

# group by gene ID and choose entry with longest CDS
annotation=annotation.dropna().groupby('gene_id').aggregate(lambda df: df.loc[df['CDSlen'].idxmax()]).reset_index()
# set "name" column in bed to "gene_id"
annotation[3]=annotation['gene_id']

# write in bed format
annotation.to_csv(sys.stdout if options.out is None else options.out,sep='\t',header=None,index=None,columns=range(12))
