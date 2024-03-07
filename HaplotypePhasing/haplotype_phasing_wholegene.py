import re
import argparse
from datetime import datetime
import numpy as np
import numpy.matlib
import sys
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import scipy
from scipy import stats
import pysam
import statistics as stats
from collections import Counter
from scipy.stats import fisher_exact
import subprocess

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'small',
  'axes.labelsize': 'medium',
  'axes.titlesize':'medium',
  'xtick.labelsize':'small',
  'ytick.labelsize':'medium'}
pylab.rcParams.update(params)
import networkx as nx

def log(info):
  now = datetime.now()
  dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
  print("["+dt_string+"] "+info)

def loadSVs(file):
  SVs = dict()
  with open(file, 'r') as fin:
    line = fin.readline().strip(" \r\n")
    while line:
      if line[0]=="#":
        line = fin.readline().strip(" \r\n")
        continue
      eles = line.split("\t")
      if eles[6] != "PASS": # Filter out these not "PASS"
        line = fin.readline().strip(" \r\n")
        continue
      chr = eles[0]
      start = int(eles[1])-1
      alleles = [eles[3]]
      isSNP=True
      for alle in eles[4].split(","):
        if len(alle) != 1:
          isSNP = False
        #   break
        # else:
        alleles.append(alle)
      if len(alleles) == 1: # or (not isSNP):
        line = fin.readline().strip(" \r\n")
        continue
      assert(len(eles[9])>3)
      if eles[9][1] == "/":
        phased=False
        eles = eles[9].split(":")[0].split("/")
      elif eles[9][1] == "|":
        phased=True
        eles = eles[9].split(":")[0].split("|")
      else:
        line = fin.readline().strip(" \r\n")
        continue
      assert (len(eles)==2 and eles[0].isdigit() and int(eles[0])>=0 and int(eles[0])<len(alleles) and eles[1].isdigit()) and int(eles[1])>=0 and int(eles[1])<len(alleles), line
      alles = {"phased":phased, "isHeter":(alleles[int(eles[0])]!=alleles[int(eles[1])]), "alleles":[alleles[int(eles[0])], alleles[int(eles[1])]]}
      SVs[chr+":"+str(start)] = alles
      line = fin.readline().strip(" \r\n")
  return SVs

#assemblies: [["chrx:pos:allele"], ]
def compareWithGIAB(SVs, assemblies):
  res=[0, [], False, False, False, False, [], []]
  #[mean(MaxConsistentProportions),
  # MaxConsistentProportions(separated by ":"),
  # True|False(variant site and allele in unphased SVs)]
  # True|False(variant site in but allele not in phased SVs)
  # True|False(variant site in but allele not in unphased SVs)
  # True|False(variant site not in SVs)
  # number of variants in assemblies(separated by ":")
  # max. number of variants consistent with either haplotypes in the annotation
  all_types = []
  for assembly in assemblies:
    types = []
    # 0 : haplotye0;
    # 1 : haplotye1;
    # 2 : variant site is het. site and allele in unphased SVs;
    # 3 : variant site is het. site but allele not in phased SVs;
    # 4 : variant site is het. site but allele not in unphased SVs;
    # 5 : variant site not in SVs;
    # 6 : variant site is homo. site and allele in (but we called heterozygous variants.)
    # 7 : variant site is homo. site, but allele not in
    for allele in assembly:
      tmp = allele
      eles = allele.split(":")
      ID = eles[0]+":"+eles[1]
      allele = eles[2]
      if ID not in SVs:
        # print('5:', tmp)
        # print(ID)
        # print(assemblies)
        types.append(5)
        continue
      svInRef = SVs[ID]
      if svInRef["isHeter"] == False:
        if allele in svInRef["alleles"]:
          types.append(6)
        else:
          types.append(7)
        continue
      if svInRef["phased"]:
        if allele not in svInRef["alleles"]:
          idx = -1
        else:
          idx = svInRef["alleles"].index(allele)
        if idx == -1:
          # print('3:', tmp)
          types.append(3)
        else:
          types.append(idx)
      else:
        if allele not in svInRef["alleles"]:
          idx = -1
        else:
          idx = svInRef["alleles"].index(allele)
        if idx == -1:
          # print('4:', tmp)
          types.append(4)
        else:
          # print('2:', tmp)
          types.append(2)
    print("\t".join([assembly[x]+"("+str(types[x])+")" for x in range(0, len(types))]))
    #process types
    all_types = all_types+types
    vsAnnotated = len(types) - types.count(5)
    score = 0
    h0 = types.count(0) + types.count(6)
    h1 = types.count(1) + types.count(6)
    if h0 + h1 > 0:
      score = max(h0, h1)/float(h0+h1)
    res[1].append(score)
    res[6].append(len(assembly))
    res[7].append(max(h0, h1))
  # res[0] = stats.mean(res[1]) if len(res[1])>0 else 0
  res[0] = max(res[1]) if len(res[1])>0 else 0
  for x in range(2, 6):
    if x in all_types:
      res[x] = True
  for x in [1, 6, 7]:
    res[x] = ":".join([str(y) for y in res[x]])
  return res

def getArguments():
  parser = argparse.ArgumentParser(description="Reference-based haplotype phasing.")
  parser.add_argument('-a', '--align', required=True, help="alignment file in bam format.")
  parser.add_argument('-r', '--ref', help="Reference genome file", required=True)
  parser.add_argument('-g', '--gRNA', help="gRNAs designed for cancer genes", required=True)
  parser.add_argument('--gene', default="", help="for a specific gene")
  parser.add_argument('-v', '--vcf', help="Variant annotation in VCF format.", required=True)
  # parser.add_argument('-i', '--inputs', required=True, help='Reads aligned to linkers. Files ordered by linker and separated by ";"')
  # parser.add_argument('-l', '--linkers', required=True, help='Linker sequences. Files ordered linker and separated by ";"')
  parser.add_argument('-rn', '--min_read_num', type=int, default=10, help="Min. number of reads supporting the top two alleles.")
  parser.add_argument('-sr', '--max_SNP_ratio', type=float, default=2, help="Max. allowed ratio for the two alleles.")
  parser.add_argument('-rr', '--min_read_ratio', type=float, default=0.5, help="To call SNPs at bases covered by at least <-rr> reads.")
  parser.add_argument('-mq', '--min_base_quality', type=int, default=10, help="Min. base quality used to call variants")
  parser.add_argument('--noRepeat', action='store_true', help='filter out these falling within repetitive seqs.')
  parser.add_argument('-rl', '--repeat_len', type=int, default=5, help="check for repeat of <repeat_len>.")
  parser.add_argument('--repeat_ratio', type=float, default=0.8, help="ration of the repeating unit.")
  parser.add_argument('--min_cov', type=int, default=5, help="min. cov. of bases to build consensus.")
  parser.add_argument('--min_ratio', type=float, default=0.5, help="min. ratio of most frequency alleles to be included in consensus.")
  parser.add_argument('--mafft', required=True, help="abs path of mafft program for MSA.")
  parser.add_argument('--min_SNPs_inHap', type=int, default=5, help="filter out haplotypes with less than <min_SNPs_inHap> SNPs.")
  parser.add_argument('-o', '--output', default="output", help='output directory (default: output)')
  parser.add_argument('--version', default="v0", help='output version (helpful if multiple runs with different parameters.)')
  args = parser.parse_args()
  return args

args = getArguments()

log("Start")
log("Load SVs by GIAB")
GIAB_SVs = loadSVs(args.vcf)
log("Load reference genome.")
chr2seq=dict()
with open(args.ref, 'r') as fin:
  line = fin.readline().strip(" \r\n")
  while line:
    name = line[1:]
    seq = ""
    line = fin.readline().strip(" \r\n")
    while line and line[0] != '>':
      seq = seq + line
      line = fin.readline().strip(" \r\n")
    chr2seq[name]=seq.upper()

log("Processing reads...")
nick_gRNAs = []
cut_gRNAs = []
genes = dict()
longGenes = set()
gene2gRNAs = dict()
gene2linkers = dict()
# linkers = []
with open(args.gRNA, 'r') as fin:
  data = list(csv.reader(fin, delimiter="\t"))
  for record in data[1:]:
    if "Long" in record[0]:
      nick_gRNAs.append(record)
      longGenes.add(record[-5])
    else:
      cut_gRNAs.append(record)
    genes[record[-5]]=[record[-3], int(record[-2]), int(record[-1])]
    if record[-5] not in gene2gRNAs:
      gene2gRNAs[record[-5]] = []
    gene2gRNAs[record[-5]].append(record)
  # for x in range(0, len(nick_gRNAs), 2):
  #   linkers.append([nick_gRNAs[x][-3], int(nick_gRNAs[x][-6]), int(nick_gRNAs[x+1][-6])])
for gene,gRNAs in gene2gRNAs.items():
  if gene not in longGenes:
    continue
  gene2linkers[gene]=[]
  for x in range(0, len(gRNAs), 2):
    gene2linkers[gene].append([gRNAs[x][-3], int(gRNAs[x][-6]), int(gRNAs[x+1][-6])])
    gene2linkers[gene].sort()
# linkers.sort()

log("Load alignment")
alignFin = pysam.AlignmentFile(args.align, 'rb')

if False:
# alignedRegions=[]
# for read in alignFin.fetch():
#   if (not read.is_unmapped) and (not read.is_duplicate) and (not read.is_secondary) and (not read.is_supplementary):
#     # alignedRead_num = alignedRead_num + 1
#     # alignedPercentages.append(read.query_alignment_length*100/float(read.query_length))
#     #if read.query_alignment_length/float(read.query_length) >= 0.95:
#     if abs(read.query_alignment_length - read.query_length) <= 50:
#       # alignedReadGood_num = alignedReadGood_num + 1
#       alignedRegions.append([read.reference_name, read.reference_start, read.reference_end])#, read.cigartuples])#0-based, right open coordinate
# alignedRegions.sort()
# print("\t".join(["Gene", "Chr", "Start", "End",
#   "alignedBaseInTotal", "DepthOfCoverage",
#   "#linkers", "#variants", "#ins", "#del", "#delFromRepeats", "#reads", "Ave.No.OfAllelesInReads",
#   "hap_lens",
#   "Max.Consi.Score", "Consi.Scores"]))
  x=1

BASES = ["A", "C", "G", "T"]
base2idx = {'A':0, 'C':1, 'G':2, 'T':3}
matchcases = ['A', 'C', 'G', 'T', "In", "De"]

for gene,Info in genes.items():
  # if gene not in longGenes:# or
  if args.gene and gene != args.gene:#"PIK3R1":#"BRCA2":#"RNF43":#
    continue
  # print(gene+str(Info))
  chr = Info[0]
  start = Info[1]
  end = Info[2]

  # compile a list of SNPs inside gene body region
  geneLen = end - start
  # chr.encode('utf-8')
  gene_covs = alignFin.count_coverage(bytes(chr, 'utf-8'), start=start, stop=end, quality_threshold=args.min_base_quality)

  SNPs = []
  for loc in range(0, end-start):
      cov_total = 0
      alleleFreq = []
      for x in range(0, 4): #ACGT
          cov_total = cov_total + gene_covs[x][loc]
          alleleFreq.append([gene_covs[x][loc], BASES[x]])
      alleleFreq.sort(reverse=True)
      # choose top 2 SNPs
      top2_cov = alleleFreq[0][0] + alleleFreq[1][0]

      if alleleFreq[1][0] > 0 and \
        top2_cov >= args.min_read_num and \
        top2_cov/float(cov_total) >= args.min_read_ratio and \
        alleleFreq[0][0]/alleleFreq[1][0] <= args.max_SNP_ratio:
          SNPs.append([loc+start, alleleFreq[0][1], alleleFreq[1][1]])

  # Filter out variants nearby repetitive regions
  log("Filter out variants nearby repetitive regions.")
  if args.noRepeat:
    SNPs_filter = []
    for ref_loc, var1, var2 in SNPs:
      inRepeat = False
      # Check one-base repeats
      # for seq_tmp in [chr2seq[chr][max(0, ref_loc-args.repeat_len):ref_loc], chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len+1, len(chr2seq[chr]))]]:
      #   rep2count = list(Counter(seq_tmp).items())
      #   rep2count.sort(key = lambda x : x[1], reverse=True)
      #   rep_ratio = rep2count[0][1]/float(len(seq_tmp))
      #   if rep_ratio >= args.repeat_ratio:
      #     inRepeat = True
      #     break
      # if inRepeat:
      #   print([chr2seq[chr][max(0, ref_loc-args.repeat_len):ref_loc], chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len+1, len(chr2seq[chr]))]])
      #   continue
      if (len(set([x for x in chr2seq[chr][max(0, ref_loc-args.repeat_len):ref_loc]])) == 1) or \
        (len(set([x for x in chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len+1, len(chr2seq[chr]))]])) == 1):
        print([chr2seq[chr][max(0, ref_loc-args.repeat_len):ref_loc], chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len+1, len(chr2seq[chr]))]])
        continue
      # Check two-bases repeats
      for seq_tmp in [chr2seq[chr][max(0, ref_loc-args.repeat_len*2):ref_loc], chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len*2+1, len(chr2seq[chr]))]]:
        patterns = set()
        for x in range(0, len(seq_tmp)-1, 2):
          patterns.add(seq_tmp[x:(x+2)])
        if len(patterns) == 1:
          inRepeat = True
          break
        # rep2count = dict()
        # for x in range(0, len(seq_tmp)-1, 2):
        #   tmp = seq_tmp[x:(x+2)]
        #   if tmp not in rep2count:
        #     rep2count[tmp] = 0
        #   rep2count[tmp] += 1
        # rep2count = [[key, value] for key, value in rep2count.items()]
        # rep2count.sort(key = lambda x : x[1], reverse=True)
        # rep_ratio = rep2count[0][1]*2/float(len(seq_tmp))
        # if rep_ratio >= args.repeat_ratio:
        #   inRepeat = True
        #   break
      if inRepeat:
        print([chr2seq[chr][max(0, ref_loc-args.repeat_len*2):ref_loc], chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len*2+1, len(chr2seq[chr]))]])
        continue
      SNPs_filter.append([ref_loc, var1, var2])
    print("[Info] filter out", len(SNPs)-len(SNPs_filter), "SNPs nearby repetitive regions.", len(SNPs_filter), "SNPs remains after filtering.")
    SNPs = SNPs_filter

  if False:
  # # cov_total >= cov_total_min and (SNP_covs[1][0]+SNP_covs[0][0])/cov_total >=cov_top2_ratio and SNP_covs[1][0]/SNP_covs[0][0] >=top2_ratio_min:
  # linkers = gene2linkers[gene]
  # linkers.sort()
  # # linkers = linkers[2:]
  # # print("Linkers" + str(linkers))
  # #compile a list of alleles inside linkers
  # log("Extract heterozygous alleles inside linkers")
  # linker2vs = {tuple(x):[] for x in linkers}
  # for linker in linkers:
  #   # print(linker)
  #   linkerLen = linker[2] - linker[1]
  #   # altAlles = dict()
  #   covOfLinker = [{x:0 for x in matchcases} for y in range(0, linkerLen)] # A C G T Insert Del
  #   for read in alignFin.fetch(chr, linker[1], linker[2]):
  #     if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
  #       continue
  #     #clippings should have already been removed from read_seq and ref_loci
  #     read_seq = read.query_sequence#query_alignment_sequence BUGs in pysam
  #     ref_loci = read.get_reference_positions(full_length=True)
  #     read_qual = read.query_qualities
  #     skip_front = 0
  #     while skip_front < len(read_seq) and (ref_loci[skip_front] == None or ref_loci[skip_front] < linker[1]):
  #       skip_front = skip_front + 1
  #     skip_end = len(read_seq)-1
  #     while skip_end >= 0 and (ref_loci[skip_end] == None or ref_loci[skip_end] >= linker[2]):
  #       skip_end = skip_end - 1
  #     read_seq = read_seq[skip_front:skip_end+1]
  #     ref_loci = ref_loci[skip_front:skip_end+1]
  #     read_qual = read_qual[skip_front:skip_end+1]
  #     assert len(read_seq) == len(ref_loci), "\n".join([str(len(ref_loci)), str(len(read.get_reference_positions(full_length=True))), str(read.query_alignment_length), read.to_string()])#"\t".join([read.query_name, str(len(read_seq)), str(len(ref_loci))])
  #     read_idx = 0
  #     # while read_idx < len(read_seq) and (ref_loci[read_idx] is not None) and ref_loci[read_idx] < linker[1]:
  #     #   read_idx = read_idx + 1
  #     while read_idx < len(read_seq):
  #       if ref_loci[read_idx] is None:
  #         read_idx = read_idx + 1
  #         continue
  #       ref_loc = ref_loci[read_idx]
  #       # if ref_loc not in altAlles:
  #       #   altAlles[ref_loc] = dict()
  #       # alles = set()
  #       # add bases as SNPs
  #       # if read_qual[read_idx] >= seq_qual_min:
  #       #   alles.add(read_seq[read_idx])
  #       stop = read_idx + 1
  #       # if read_seq[read_idx] == chr2seq[chr][ref_loc]:
  #         #mark the end of the read
  #         # if read_idx == len(read_seq)-1:
  #         #   alles.add(read_seq[read_idx]+"$")
  #       # insertion
  #       while stop < len(read_seq) and (ref_loc == ref_loci[stop] or ref_loci[stop] is None):
  #         stop = stop + 1
  #       if stats.mean(read_qual[read_idx:stop]) < args.min_base_quality:
  #         read_idx = stop
  #         continue
  #       # insertion
  #       if stop > read_idx+1:
  #         covOfLinker[ref_loc-linker[1]]["In"] +=1 # insertion
  #
  #       # allele = read_seq[read_idx:stop].upper()
  #       # mark insertion as "-"
  #       # if len(allele) > 1:
  #       #   allele = allele[0]+"+"
  #
  #       # deletion
  #       elif stop == read_idx+1 and stop < len(read_seq) and ref_loci[stop] > ref_loci[read_idx]+1:
  #         # allele = allele + "*"*(ref_loci[stop]-ref_loci[read_idx]-1)
  #         # allele = allele[0]+"-"
  #         for x in range(ref_loc+1, ref_loci[stop]):
  #           covOfLinker[x-linker[1]]["De"] +=1 # deletion
  #       else:
  #         # base match or mismatch
  #         covOfLinker[ref_loc-linker[1]][read_seq[read_idx].upper()] +=1
  #       # if stop > read_idx + 1 and stats.mean(read_qual[read_idx:stop]) >= seq_qual_min:
  #       #   alles.add(read_seq[read_idx:stop])
  #       #deletion
  #       # if read_idx+1 < len(read_seq) and (ref_loci[read_idx+1] is not None) and ref_loci[read_idx+1] > ref_loci[read_idx]+1 and read_qual[read_idx] >= seq_qual_min:
  #       #   alles.add(read_seq[read_idx]+"".join(["*" for x in range(0, ref_loci[read_idx+1]-ref_loci[read_idx]-1)]))
  #       # for alle in alles:
  #       # if allele not in altAlles[ref_loc]:
  #       #   altAlles[ref_loc][allele] = 0
  #       # altAlles[ref_loc][allele] = altAlles[ref_loc][allele] + 1
  #       read_idx = stop
  #   # compile a list of heterozygous variants
  #   # if linker[1] == 32388404:
  #   #   for ref_loc in range(32388540, 32388545):
  #   #     print(ref_loc, covOfLinker[ref_loc-linker[1]])
  #   #   sys.exit(0)
  #   ref_loc = linker[1]
  #   while ref_loc < linker[2]:
  #     alles = covOfLinker[ref_loc-linker[1]]
  #     read_num = 0
  #     for alle,count in alles.items():
  #       read_num = read_num + count
  #     alleleFreq = [[alle, count] for alle,count in alles.items()]
  #     alleleFreq.sort(key=lambda x:x[1], reverse=True)
  #     if alleleFreq[1][1] == 0:
  #       ref_loc += 1
  #       continue
  #     top2_cov = alleleFreq[0][1] + alleleFreq[1][1]
  #     top2_ratio = float(alleleFreq[0][1])/alleleFreq[1][1]
  #     if ref_loc == 32388541:
  #       print(top2_cov, top2_ratio)
  #     if top2_cov >= args.min_read_num and top2_cov/float(read_num) >= args.min_read_ratio and alleleFreq[0][1]/alleleFreq[1][1] <= args.max_SNP_ratio: #and alleleFreq[1][1]/alleleFreq[0][1] <= args.max_SNP_ratio:
  #       # # filter out repetitive regions
  #       # seq_before = chr2seq[chr][max(0, ref_loc-args.repeat_len):ref_loc]
  #       # bases_before = list(Counter(seq_before).items())
  #       # bases_before.sort(key = lambda x : x[1], reverse=True)
  #       # before_top_ratio = bases_before[0][1]/float(len(seq_before))
  #       # # before_other_freq = sum([x[1] for x in bases_before[1:]])
  #       #
  #       # seq_after = chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len+1, len(chr2seq[chr]))]
  #       # bases_after = list(Counter(seq_after).items())
  #       # bases_after.sort(key = lambda x : x[1], reverse=True)
  #       # after_top_ratio = bases_after[0][1]/float(len(seq_after))
  #       # # after_other_freq = sum([x[1] for x in bases_after[1:]])
  #       # # base_before = {consensus_seq[x] for x in range(max(0, loc-repetitive_seq_len), loc)}
  #       # # base_after = {consensus_seq[x] for x in range(loc+1, min(loc+repetitive_seq_len+1, consensus_len))}
  #       # #consider only the first base inside each allele
  #       # # if args.noRepeat and ((len(base_before)==1 and base_before.issubset({alleleFreq[0][0][0], alleleFreq[1][0][0]})) or (len(base_after)==1 and base_after.issubset({alleleFreq[0][0][0], alleleFreq[1][0][0]}))):
  #       # if args.noRepeat and (len(alleleFreq[0][0]) > 1 or len(alleleFreq[1][0]) > 1) and (before_top_ratio >= args.repeat_ratio or after_top_ratio >= args.repeat_ratio):
  #       #   continue
  #
  #       # when indel found, extend to multiple bases
  #       cases = {alleleFreq[0][0], alleleFreq[1][0]}
  #       # if ref_loc == 32388541:
  #       #   print(cases)
  #       # find deletion center
  #       if ("De" in cases) and ("In" not in cases):
  #         next_loc = ref_loc + 1
  #         # if linker[1] == 68271342:
  #         #   print(ref_loc)
  #         while next_loc < linker[2]:
  #           alles = covOfLinker[next_loc-linker[1]]
  #           read_num_tmp = 0
  #           for alle,count in alles.items():
  #             read_num_tmp = read_num_tmp + count
  #           alleleFreq_tmp = [[alle, count] for alle,count in alles.items()]
  #           alleleFreq_tmp.sort(key=lambda x:x[1], reverse=True)
  #           if alleleFreq_tmp[1][1] == 0:
  #             break
  #           top2_cov_tmp = alleleFreq_tmp[0][1] + alleleFreq_tmp[1][1]
  #           top2_ratio_tmp = float(alleleFreq_tmp[0][1])/alleleFreq_tmp[1][1]
  #           if top2_cov_tmp >= args.min_read_num  and \
  #             top2_ratio_tmp <= args.min_read_ratio and \
  #             (((alleleFreq_tmp[0][0] in cases) and len(alleleFreq_tmp[0][0])>1) or \
  #              ((alleleFreq_tmp[1][0] in cases) and len(alleleFreq_tmp[1][0])>1)):
  #             if top2_ratio_tmp < top2_ratio and \
  #               top2_cov_tmp/float(read_num_tmp) >= args.min_read_ratio:
  #               # top2_cov_tmp >= args.min_read_num and \
  #               # if linker[1] == 68271342:
  #               #   print("Replace", ref_loc, "by", next_loc)
  #               ref_loc = next_loc
  #               top2_ratio = top2_ratio_tmp
  #               alleleFreq = alleleFreq_tmp
  #             next_loc += 1
  #           else:
  #             break
  #         linker2vs[tuple(linker)].append([ref_loc, alleleFreq[0][0], alleleFreq[1][0]])
  #         ref_loc = next_loc
  #       elif ("De" not in cases) and ("In" in cases):
  #         # if len(alleleFreq[0][0]) == 1:
  #         #   linker2vs[tuple(linker)].append([ref_loc, alleleFreq[0][0], alleleFreq[0][0]+"+"])
  #         # else:
  #         #   linker2vs[tuple(linker)].append([ref_loc, alleleFreq[1][0], alleleFreq[1][0]+"+"])
  #         linker2vs[tuple(linker)].append([ref_loc, alleleFreq[0][0], alleleFreq[1][0]])
  #         ref_loc += 1
  #       elif ("De" not in cases) and ("In" not in cases):
  #         linker2vs[tuple(linker)].append([ref_loc, alleleFreq[0][0], alleleFreq[1][0]])
  #         ref_loc += 1
  #       else:
  #         ref_loc += 1
  #     else:
  #       ref_loc += 1
    # for ref_loc, alles in altAlles.items():
    #   if len(alles) < 2:
    #     continue
    #   read_num = 0
    #   for alle,count in alles.items():
    #     # if len(alle) == 1:
    #     read_num = read_num + count
    #   alleleFreq = [[alle, count] for alle,count in alles.items()]
    #   alleleFreq.sort(key=lambda x:x[1], reverse=True)
    #
    #   top2_cov = alleleFreq[0][1] + alleleFreq[1][1]
    #   if top2_cov >= args.min_read_num and top2_cov/float(read_num) >= args.min_read_ratio and alleleFreq[0][1]/alleleFreq[1][1] <= args.max_SNP_ratio and alleleFreq[1][1]/alleleFreq[0][1] <= args.max_SNP_ratio:
    #     seq_before = chr2seq[chr][max(0, ref_loc-args.repeat_len):ref_loc]
    #     bases_before = list(Counter(seq_before).items())
    #     bases_before.sort(key = lambda x : x[1], reverse=True)
    #     before_top_ratio = bases_before[0][1]/float(len(seq_before))
    #     # before_other_freq = sum([x[1] for x in bases_before[1:]])
    #
    #     seq_after = chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len+1, len(chr2seq[chr]))]
    #     bases_after = list(Counter(seq_after).items())
    #     bases_after.sort(key = lambda x : x[1], reverse=True)
    #     after_top_ratio = bases_after[0][1]/float(len(seq_after))
    #     # after_other_freq = sum([x[1] for x in bases_after[1:]])
    #     # base_before = {consensus_seq[x] for x in range(max(0, loc-repetitive_seq_len), loc)}
    #     # base_after = {consensus_seq[x] for x in range(loc+1, min(loc+repetitive_seq_len+1, consensus_len))}
    #     #consider only the first base inside each allele
    #     # if args.noRepeat and ((len(base_before)==1 and base_before.issubset({alleleFreq[0][0][0], alleleFreq[1][0][0]})) or (len(base_after)==1 and base_after.issubset({alleleFreq[0][0][0], alleleFreq[1][0][0]}))):
    #     if args.noRepeat and (len(alleleFreq[0][0]) > 1 or len(alleleFreq[1][0]) > 1) and (before_top_ratio >= args.repeat_ratio or after_top_ratio >= args.repeat_ratio):
    #       continue
    #     linker2vs[tuple(linker)].append([ref_loc, alleleFreq[0][0], alleleFreq[1][0]])
        # SNPs.append([loc, {alleleFreq[0][0]:set(), alleleFreq[1][0]:set()}])
      # if len(base2cov[0][1]) == 1 and len(base2cov[1][1]) == 1:
      #   if base2cov[0][0] >= allele_min_cov and base2cov[1][0] >= allele_min_cov and base2cov[0][0]+base2cov[1][0] >= top2Cov_min*allbasecov and (base2cov[1][0]>0 and base2cov[0][0]/base2cov[1][0] <= top2_dif_max):
      #     # print(str(ref_loc)+" : "+str(base2cov))
      #     linker2vs[tuple(linker)].append([ref_loc, base2cov[0][1], base2cov[1][1]])
      # elif (len(base2cov[0][1]) == 1 and base2cov[0][1] == chr2seq[chr][ref_loc]) or (len(base2cov[1][1]) == 1 and base2cov[1][1] == chr2seq[chr][ref_loc]):
      #   indel_cov = 0
      #   for alle,count in alles.items():
      #     if len(alle) > 1:
      #       indel_cov = indel_cov + count
      #   read_end = chr2seq[chr][ref_loc]+"$"
      #   read_end_count = 0
      #   if read_end in alles:
      #      read_end_count = alles[read_end]
      #   if len(base2cov[0][1]) == 1:
      #     base2cov[0][0] = base2cov[0][0] - indel_cov - read_end_count
      #   if len(base2cov[1][1]) == 1:
      #     base2cov[1][0] = base2cov[1][0] - indel_cov - read_end_count
      #   #base2cov =
      #   base2cov[0:2].sort(reverse=True)
      #   #base2cov.sort(reverse=True)
      #   if base2cov[0][0] >= allele_min_cov and base2cov[1][0] >= allele_min_cov and base2cov[0][0]+base2cov[1][0] >= top2Cov_min*allbasecov and (base2cov[1][0]>0 and base2cov[0][0]/base2cov[1][0] <= top2_dif_max):
      #     # print(str(ref_loc)+" : "+str(base2cov))
      #     linker2vs[tuple(linker)].append([ref_loc, base2cov[0][1], base2cov[1][1]])

    # readNum = alignFin.count(linker[0], linker[1], linker[2])#number of reads aligned to the linker region
    # coverage = alignFin.count_coverage(linker[0], linker[1], linker[2])#, quality_threshold=0)
    # coverage_noFilter = alignFin.count_coverage(linker[0], linker[1], linker[2], quality_threshold=0)
    # for base in range(0, linkerLen):
    #   base2cov = [[coverage[x][base], BASES[x]] for x in range(0, 4)]
    #   base2cov_noFilter = [[coverage_noFilter[x][base], BASES[x]] for x in range(0, 4)]
    #   allbasecov = sum([coverage[x][base] for x in range(0, 4)])
    #   base2cov.sort(reverse = True)
    #   base2cov_noFilter.sort(reverse = True)
    #   if base2cov[0][0]+base2cov[1][0] >= top2Cov_min*allbasecov and (base2cov[1][0]>0 and base2cov[0][0]/base2cov[1][0] <= top2_dif_max and base2cov_noFilter[0][0]/base2cov_noFilter[1][0] <= top2_dif_max):
    #     linker2vs[tuple(linker)].append([base+linker[1], base2cov[0][1], base2cov[1][1]])
  # VSs = set()#variant sites
  # VS_NumInTotal = 0
  # # VS_del_fromRepeats_NumInTotal = 0
  # VS_ins_NumInTotal = 0
  # VS_del_NumInTotal = 0
  # for linker,vss in linker2vs.items():
  #   # for each linker, select the one SV with the higest confidence (SNP > InDels).
  #   vss_SNPs = []
  #   if len(vss) > 1:
  #     for vs in vss:
  #       if len(vs[1])==1 and len(vs[2])==1:
  #         vss_SNPs.append(vs)
  #   if len(vss_SNPs) > 0:
  #     vss = vss_SNPs
  #   for vs in vss:
  #     VS_NumInTotal = VS_NumInTotal + 1
  #     #filter out deletions detected near the repetitve loci
  #     # FromRepeats = False
  #     # if "*" in vs[1] or "*" in vs[2]:
  #     #   # FromRepeats=True
  #     #   # continue #To skip indels
  #     #   VS_del_NumInTotal = VS_del_NumInTotal + 1
  #     #   rep_len = 3
  #     #   nearbySeq = chr2seq[chr][min(vs[0]+1, len(chr2seq[chr])):min(vs[0]+rep_len+1, len(chr2seq[chr]))]
  #     #   nearbySeq = nearbySeq.upper()
  #     #   chars = {nearbySeq[x] for x in range(0, len(nearbySeq))}
  #     #   if len(chars) == 1:
  #     #     FromRepeats = True
  #     #     VS_del_fromRepeats_NumInTotal = VS_del_fromRepeats_NumInTotal + 1
  #     #     continue
  #     # elif len(vs[1]) > 1 or len(vs[2]) > 1:
  #     #   VS_ins_NumInTotal = VS_ins_NumInTotal + 1
  #     # if not FromRepeats:
  #     VSs.add(tuple(vs))
  # VSs = list(VSs)
  # VSs.sort()
  # print("Variant sites: "+str(VSs))
  # for linker,svs in linker2vs.items():
  #   print(linker, svs)
  # #   print(svs)
  # # print(linker2vs)
  # sys.exit(0)
    x=1

  VSs = SNPs
  VSs.sort()

  log("Count co-ocurrence of heterozygous variants inside reads.")
  # reads = []
  # for read in alignFin.fetch(chr, start, end):
  #   if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
  #     continue
  #   #use read only if it overlaps with at least one linkers
  #   for linker in linkers:
  #     if len(linker2vs[tuple(linker)]) == 0:
  #       continue
  #     if read.reference_start < linker[2] and linker[1] < read.reference_end:
  #       reads.append(read)#[read.reference_name, read.reference_start, read.reference_end])
  #       break
  # alignedBaseInTotal = sum([min(read.reference_end, end)-max(read.reference_start, start) for read in reads])
  # log(str(alignedBaseInTotal)+" bp aligned to gene in total, and with an average coverage of depth: "+str(alignedBaseInTotal/(end-start)))
  # log()
  readWithalles = []#[read_name, [[pos,seq], [pos,seq]]]
  for read in alignFin.fetch(chr, start, end): #reads:
    # seq = read.query_alignment_sequence#query_sequence
    #print(read.cigarstring)
    #print(read.get_reference_positions(full_length=True))
    # print(len(read.query_alignment_sequence))
    # print("Length of loci:")
    # print(len(read.get_reference_positions(full_length=True)))
    # sys.exit(0)
    # ref_pos = read.get_reference_positions()#full_length=True)
    read_seq = read.query_sequence#query_alignment_sequence BUGs in pysam
    ref_loci = read.get_reference_positions(full_length=True)
    read_qual = read.query_qualities
    skip_front = 0
    while skip_front < len(read_seq) and ref_loci[skip_front] == None:
      skip_front = skip_front + 1
    skip_end = len(read_seq)-1
    while skip_end >= 0 and ref_loci[skip_end] == None:
      skip_end = skip_end - 1
    read_seq = read_seq[skip_front:skip_end+1]
    ref_loci = ref_loci[skip_front:skip_end+1]
    read_qual = read_qual[skip_front:skip_end+1]
    if ref_loci[0] > ref_loci[-1]:
      print(ref_loci)
      sys.exit(0)
    alles = []
    for vs in VSs:
      if vs[0] in ref_loci:
        pos = ref_loci.index(vs[0])
        allele = read_seq[pos].upper()
        if allele in vs[1:]:# and read_qual[pos] >= args.min_base_quality:
          alles.append([vs[0], allele])
      if False:
      # if vs[0] != 89968575:
      #   continue
      # if vs[0] in ref_loci: # if vs[0] can be found, then it is not an deletion
      #   pos = ref_loci.index(vs[0])
      #   # if vs[0] == 68271585:
      #   #   print(read.query_name, ref_loci[pos-1:pos+5], read_seq[pos-1:pos+5])
      #   #   sys.exit(0)
      #   allele = read_seq[pos].upper()
      #   if (allele not in vs[1:]) or (stats.mean(read_qual[pos:pos+1]) < args.min_base_quality):
      #     continue
      #   if "In" in vs[1:]:
      #     stop = pos + 1
      #     while stop < len(read_seq) and ((ref_loci[stop] is None) or ref_loci[stop]==ref_loci[pos]):
      #       stop = stop + 1
      #     if (stop > pos + 1) and stats.mean(read_qual[pos:stop]) >= args.min_base_quality: # insertion
      #       allele="In"
      #       alles.append([vs[0], allele])
      #   else:
      #     alles.append([vs[0], allele])
      #   # stop = pos + 1
      #   #
      #   # # need to trace back for insertion
      #   # pos_tmp = pos - 1
      #   # while pos_tmp > 0 and ((ref_loci[pos_tmp] is None) or ref_loci[pos_tmp]==ref_loci[pos]):
      #   #   pos_tmp = pos_tmp - 1
      #   # while ref_loci[pos_tmp]!=ref_loci[pos]:
      #   #   pos_tmp = pos_tmp + 1
      #
      #   # while stop < len(read_seq) and ((ref_loci[stop] is None) or ref_loci[stop]==ref_loci[pos]):
      #   #   stop = stop + 1
      #   # if stop > pos_tmp + 1: # insertion
      #   #   allele="In"
      #   # # allele = read_seq[pos:stop]
      #   # # if len(allele) > 1:
      #   # #   allele = allele[0]+"+"
      #   # #deletions
      #   # if pos+1 == stop and stop < len(read_seq) and ref_loci[stop] > ref_loci[pos]+1:
      #   #   # allele = allele + "*"*(ref_loci[stop]-ref_loci[pos]-1)
      #   #   allele = "De"#allele[0]+"-"
      #   # # if stop > read_idx + 1 and stats.mean(read_qual[read_idx:stop]) >= seq_qual_min:
      #   # #   alles.add(read_seq[read_idx:stop])
      #   # if (allele in vs[1:]) and stats.mean(read_qual[pos:stop]) >= args.min_base_quality:
      #   #   alles.append([vs[0], allele])
      # else:
      #   if ref_loci[0] > vs[0] or ref_loci[-1] < vs[0]:
      #     continue
      #   tmp_array = np.array(ref_loci)
      #   tmp_array[tmp_array==None] = 1000000000
      #   difference_array = np.absolute(tmp_array-vs[0])
      #   pos = difference_array.argmin()
      #   stop = pos + 1
      #   while pos >= 0 and ((ref_loci[pos] is None) or (ref_loci[pos] > vs[0])):
      #     pos = pos - 1
      #   if (ref_loci[pos+1] > ref_loci[pos]+1) and (ref_loci[pos+1] > vs[0]):
      #     allele = "De"
      #   # else:
      #   #   allele = "In"
      #     if (allele in vs[1:]) and stats.mean(read_qual[pos:pos+1]) >= args.min_base_quality:
      #       alles.append([vs[0], allele])
      #       # if vs[0] == 68271585:
            #   print(read.query_name, ref_loci[pos-1:pos+5], read_seq[pos-1:pos+5])
            #   sys.exit(0)
        # pos = stop

        # if read_seq[pos].upper() in vs[1:]:
        #   #could be indels
        #   stop = pos + 1
        #   if stop < len(read_seq) and (ref_loci[stop] is None):
        #     while stop < len(read_seq) and (ref_loci[stop] is None):
        #       stop = stop + 1
        #     alle = read_seq[pos:stop].upper()
        #   elif stop < len(read_seq) and ref_loci[stop]-ref_loci[pos] > 1:
        #     alle = read_seq[pos].upper()+"".join(["*" for x in range(0, ref_loci[stop]-ref_loci[pos]-1)])
        #   else:
        #     alle = read_seq[pos].upper()
        #   if alle in vs[1:] and stats.mean(read_qual[pos:stop]) >= seq_qual_min:
        #     alles.append([vs[0], alle])
        x=1
    # if len(alles) > 1:
    readWithalles.append([read, alles])
  log(str(len(readWithalles))+" reads with at least one allele.")
  alignedBaseInTotal = sum([min(read.reference_end, end)-max(read.reference_start, start) for read, tmp in readWithalles])
  log(str(alignedBaseInTotal)+" bp aligned to gene in total, and with an average coverage of depth: "+str(alignedBaseInTotal/(end-start)))

  read_num = len(readWithalles)
  # sys.exit(0)
  pos2SNPs = dict()
  for record in readWithalles:
    if len(record[1]) <= 1:
      continue
    for item in record[1]:
      pos = item[0]
      allele = item[1]
      if pos not in pos2SNPs.keys():
        pos2SNPs[pos] = set()
      pos2SNPs[pos].add(allele)
  # for pos, snps in pos2SNPs.items():
  #   if len(snps) <= 1:
  #     del pos2SNPs[pos]
  # print(pos2SNPs[68271585])
  adjacency = dict()
  node2degree = dict()
  for record in readWithalles:
    for x in range(0, len(record[1])):
      nodeX = ":".join([str(t) for t in record[1][x]])
      if nodeX not in node2degree:
        node2degree[nodeX] = 0
      node2degree[nodeX] = node2degree[nodeX] + 1
      for y in range(x+1, len(record[1])):
        nodeY = ":".join([str(t) for t in record[1][y]])
        if nodeX not in adjacency:
          adjacency[nodeX] = dict()
        if nodeY not in adjacency[nodeX]:
          adjacency[nodeX][nodeY] = 0
        adjacency[nodeX][nodeY] = adjacency[nodeX][nodeY] + 1
  # print(node2degree["68271585:De"])
  # build haplotye graph
  # node2degree = dict()
  # for nodeX,nexts in adjacency.items():
  #   for nodeY, w in nexts.items():
  #     if nodeX not in node2degree:
  #       node2degree[nodeX] = 0
  #     if nodeY not in node2degree:
  #       node2degree[nodeY] = 0
  #     node2degree[nodeX] = node2degree[nodeX] + w
  #     node2degree[nodeY] = node2degree[nodeY] + w

  if len(pos2SNPs) == 0:
    print("\t".join([str(x) for x in [gene, chr, start, end,
      alignedBaseInTotal, alignedBaseInTotal/(end-start),
      len(linkers), VS_NumInTotal, VS_ins_NumInTotal, VS_del_NumInTotal, VS_del_fromRepeats_NumInTotal, read_num, (stats.mean([len(x[1]) for x in readWithalles]) if len(readWithalles)>0 else 0)
      ]]))
    continue

  log("draw haplotype graph")
  hapgraph={"Nodes":[], "Edges":[]}
  width_unit = 2
  node_dis_hor = 5
  node_dis_ver = 2
  nodes = []
  poses = dict()
  labels = {}
  X = 0
  Y = 0
  vs_site = 0
  min_weight = 5
  posInAscending = {89949812}
  pos_labels= []
  for pos in sorted(pos2SNPs.keys()):
    # plt.text(X, -node_dis_ver/4, f"{pos:,}", ha="center", va="center")
    SNPs = pos2SNPs[pos]
    SNPs = list(SNPs)
    if pos in posInAscending:
      SNPs.sort()
    else:
      SNPs.sort(reverse=True)
    # Extra round of node weight checking
    check_pass = False
    if len(SNPs) == 2:
      node1 = str(pos)+":"+SNPs[0]
      node2 = str(pos)+":"+SNPs[1]
      if (node1 in node2degree) and (node2 in node2degree):
        if node2degree[node1]/node2degree[node2]<=args.max_SNP_ratio and node2degree[node2]/node2degree[node1]<=args.max_SNP_ratio:
          check_pass = True
    if not check_pass:
      continue
    pos_labels.append([X, -node_dis_ver/4, f"{pos:,}"])
    for SNP in SNPs:
      node = str(pos)+":"+SNP
      nodes.append(node)
      poses[node] = (X, Y)
      if node in node2degree:
        labels[node] = SNP+"("+str(node2degree[node])+")"
      Y = Y + node_dis_ver
    X = X + node_dis_hor
    vs_site = vs_site + 1
    if vs_site % 2 == 1:
      Y = node_dis_ver/3
    else:
      Y = 0
    hapgraph["Nodes"].append([pos, SNPs[0], SNPs[1]])
  fig = plt.figure(figsize=((len(pos2SNPs))*width_unit, width_unit))
  G=nx.Graph()
  G.add_nodes_from(nodes)
  for pos_label in pos_labels:
    plt.text(pos_label[0], pos_label[1], pos_label[2], ha="center", va="center")
  for node, nexts in adjacency.items():
    for next, w in nexts.items():
      if (node not in nodes) or (next not in nodes):
        continue
      if w >= min_weight:
        G.add_edge(node, next, weight=w, color='b')
      else:
        G.add_edge(node, next, weight=w, color='r')
      hapgraph["Edges"].append([node, next, w])
  edges = G.edges()
  colors = [G[u][v]['color'] for u,v in edges]
  nx.draw(G, pos=poses, labels = labels, with_labels=True, node_color="white", node_size=1300, edge_color=colors)#"whitesmoke"
  nx.draw_networkx_edge_labels(G, pos=poses, edge_labels=nx.get_edge_attributes(G,'weight'), label_pos=0.6)
  plt.savefig(args.output+gene+"_raw_coloredByWeight"+str(min_weight)+"_"+args.version+".pdf")
  plt.close()

  fig = plt.figure(figsize=((len(pos2SNPs))*width_unit, width_unit))
  G=nx.Graph()
  G.add_nodes_from(nodes)
  for pos_label in pos_labels:
    plt.text(pos_label[0], pos_label[1], pos_label[2], ha="center", va="center")
  for node, nexts in adjacency.items():
    for next, w in nexts.items():
      if (node not in nodes) or (next not in nodes):
        continue
      if w >= min_weight:
        G.add_edge(node, next, weight=w, color='grey')
      # else:
      #   G.add_edge(node, next, weight=w, color='r')
  edges = G.edges()
  colors = [G[u][v]['color'] for u,v in edges]
  nx.draw(G, pos=poses, labels = labels, with_labels=True, node_color="white", node_size=1300, edge_color=colors)#"whitesmoke"
  nx.draw_networkx_edge_labels(G, pos=poses, edge_labels=nx.get_edge_attributes(G,'weight'), label_pos=0.6)
  plt.savefig(args.output+gene+"_raw_filterWeight"+str(min_weight)+"_"+args.version+".pdf")
  plt.close()
  # sys.exit(0)

  log("Resolve haplotypes")
  # Resolve haplotyes
  nodes_enoughWeight = set()
  for nodex, nexts in adjacency.items():
    for nodey, w in nexts.items():
      if w >= min_weight:
        nodes_enoughWeight.add(nodex)
        nodes_enoughWeight.add(nodey)
  nodes = list(nodes_enoughWeight)
  nodes_withHetAllele = dict()
  for node in nodes:
    eles = node.split(":")
    ID = eles[0]
    if ID not in nodes_withHetAllele:
      nodes_withHetAllele[ID] = set()
    nodes_withHetAllele[ID].add(eles[1])
  nodes_ = set()
  for pos,bases in nodes_withHetAllele.items():
    if len(bases) == 2:
      for base in bases:
        nodes_.add(pos+":"+base)
  nodes = list(nodes_)
  nodes.sort()

  # correlate position to node
  loc2nodes = dict()
  for node in nodes:
    eles = node.split(":")
    loc = eles[0]
    if loc not in loc2nodes:
      loc2nodes[loc] = set()
    loc2nodes[loc].add(node)

  # resolve haplotypes by resolving neighboring SNP pairs
  labels = {node:set() for node in nodes}
  assert (nodes[0].split(":")[0] == nodes[1].split(":")[0]), "first two nodes are not at same loc."
  labels[nodes[0]].add(0)
  labels[nodes[1]].add(1)
  hap_id = 2
  node_idx = 0
  while node_idx < len(nodes)-1:
    if nodes[node_idx].split(":")[0] != nodes[node_idx+1].split(":")[0]:
      print(nodes[node_idx], nodes[node_idx+1])
      node_idx += 1
      continue
    # assert (nodes[node_idx].split(":")[0] == nodes[node_idx+1].split(":")[0]), "Two nodes are not at same loc."+str(nodes[node_idx])+str(nodes[node_idx+1])
    # if no label, assign one
    if len(labels[nodes[node_idx]]) == 0:
      labels[nodes[node_idx]].add(hap_id)
      hap_id += 1
    if len(labels[nodes[node_idx+1]]) == 0:
      labels[nodes[node_idx+1]].add(hap_id)
      hap_id += 1
    # propagrate the haplotype ID
    x1 = nodes[node_idx]
    x2 = nodes[node_idx+1]
    if (x1 not in adjacency) or (x2 not in adjacency):
      # handle single linkage case
      if x1 in adjacency:
        nextnodes = dict()
        for node in adjacency[x1].keys():
          tmp = node.split(":")[0]
          if tmp not in nextnodes:
            nextnodes[tmp] = 0
          nextnodes[tmp] += 1
        for node, count in adjacency[x1].items():
          tmp = node.split(":")[0]
          if nextnodes[tmp] == 1 and count >= args.min_read_num/2:
            labels[node].update(labels[x1])
            # if node in ["68262796:G", "68262796:A"]:
            #   print(x1, x2, node)
            #   print(adjacency[x1])
            loc = node.split(":")[0]
            if len(loc2nodes[loc])==2:
              tmp = list(loc2nodes[loc])
              theOtherNode = tmp[0] if tmp[1]==node else tmp[1]
              labels[theOtherNode].update(labels[x2])
      elif x2 in adjacency:
        nextnodes = dict()
        for node in adjacency[x2].keys():
          tmp = node.split(":")[0]
          if tmp not in nextnodes:
            nextnodes[tmp] = 0
          nextnodes[tmp] += 1
        for node, count in adjacency[x2].items():
          tmp = node.split(":")[0]
          if nextnodes[tmp] == 1 and count >= args.min_read_num/2:
            labels[node].update(labels[x2])
            # if node in ["68262796:G", "68262796:A"]:
            #   print(x1, x2, node)
            #   print(adjacency[x2])
            loc = node.split(":")[0]
            if len(loc2nodes[loc])==2:
              tmp = list(loc2nodes[loc])
              theOtherNode = tmp[0] if tmp[1]==node else tmp[1]
              labels[theOtherNode].update(labels[x1])
      node_idx += 2
      continue
    # find common nodes
    nextnodes_common = set(adjacency[x1].keys())
    nextnodes_common.update(set(adjacency[x2].keys()))
    # for node in adjacency[x1].keys():
    #   if node in adjacency[x2].keys():
        # nextnodes_common.add(node)
    nextnodes_common = list(nextnodes_common)
    nextnodes_common.sort()
    # find common pairs
    nextnodes_pairs = []
    x = 0
    while x < len(nextnodes_common)-1:
      if nextnodes_common[x].split(":")[0] == nextnodes_common[x+1].split(":")[0]:
        nextnodes_pairs.append([nextnodes_common[x], nextnodes_common[x+1]])
        x += 2
      else:
        x += 1
    # examine each pair
    for y1, y2 in nextnodes_pairs:
      table = [[adjacency[x1][y1] if y1 in adjacency[x1] else 0,
                adjacency[x1][y2] if y2 in adjacency[x1] else 0],
               [adjacency[x2][y1] if y1 in adjacency[x2] else 0,
                adjacency[x2][y2] if y2 in adjacency[x2] else 0]]
      if not ((table[0][0]>0 and table[1][1]) or (table[0][1]>0 and table[1][0])):
        continue
      if table[0][0]+table[0][1]+table[1][0]+table[1][1] < args.min_read_num:
        continue
      res = fisher_exact(table, alternative='two-sided')
      # print(table)
      # print(res)
      if res[1] < 0.05:
        # propagate idx
        if table[0][0] == table[0][1] and table[1][0] == table[1][1]:
          # print("Weried case", table)
          continue
        elif table[0][0] > table[0][1] and table[1][0] < table[1][1]:
          x1_2 = y1
          x2_2 = y2
        elif table[0][0] < table[0][1] and table[1][0] > table[1][1]:
          x1_2 = y2
          x2_2 = y1
        else:
          continue
          # print("Weried case", table)
          # sys.exit(0)
        # if x1_2 in ["68262820:A", "68262820:G"]:
        #   print(x1, x2, x1_2, x2_2)
        #   print(table)
        labels[x1_2].update(labels[x1])
        labels[x2_2].update(labels[x2])
    node_idx += 2

  if False:
    # # propagate the haplotype idx
    # labels = [set() for x in nodes]
    # node_idx = 0
    # hap_idx = 0
    # while node_idx < len(nodes):
    #   if len(labels[node_idx]) == 0:
    #     labels[node_idx].add(hap_idx)
    #     nexts = {nodes[node_idx]}
    #     while len(nexts) != 0:
    #       nexts_ = set()
    #       for n1 in nexts:
    #         if n1 not in adjacency.keys():
    #           continue
    #         for n2, w in adjacency[n1].items():
    #           if w >= min_weight:
    #             idx2 = nodes.index(n2)
    #             labels[idx2].add(hap_idx)
    #             nexts_.add(n2)
    #       nexts = nexts_
    #     hap_idx = hap_idx + 1
    #   node_idx = node_idx + 1
    x =1

  # # print the label assignment and adjacency of haplotype graph
  # for node in nodes:
  #   if node in adjacency:
  #     print(node, labels[node], adjacency[node])
  #   else:
  #     print(node, labels[node])

  uniquelabels = set()
  uniquelabelsets = set()
  for node,label in labels.items():
    uniquelabels.update(label)
    uniquelabelsets.add(tuple(label))
  print("All possible labels:", uniquelabels)
  print("All possible label sets:", uniquelabelsets)

  nodelabel2hapID = dict()
  # equallabels = []
  hap_idx = 0
  uniquelabelsets_remained = list(uniquelabelsets)
  while len(uniquelabelsets_remained) > 0 :
    labelset = uniquelabelsets_remained[0]
    haslabel=set()
    for label in labelset:
      if label in nodelabel2hapID:
        haslabel.add(nodelabel2hapID[label])
    if len(haslabel) == 0:
      hapID = hap_idx
      hap_idx += 1
    elif len(haslabel) == 1:
      hapID = list(haslabel)[0]
    else:
      print("Werid case")
      sys.exit(0)
    for label in labelset:
      if label not in nodelabel2hapID:
        nodelabel2hapID[label] = hapID
    tmp_remained = []
    for tmp in uniquelabelsets_remained[1:]:
      hasID = False
      for label in tmp:
        if label in nodelabel2hapID:
          hapID = nodelabel2hapID[label]
          hasID = True
          break
      if not hasID:
        tmp_remained.append(tmp)
        continue
      for label in tmp:
        if label not in nodelabel2hapID:
          nodelabel2hapID[label] = hapID
    uniquelabelsets_remained = tmp_remained
  print("Node label to hapID", nodelabel2hapID)

  node2hapID = dict()
  for node, labelset in labels.items():
    hapIDs = set([nodelabel2hapID[label] for label in labelset])
    if len(hapIDs) > 1:
      print(hapIDs)
      sys.exit(0)
    node2hapID[node] = list(hapIDs)[0]

  hapID2nodes = dict()
  for node,hapID in node2hapID.items():
    if hapID not in hapID2nodes:
      hapID2nodes[hapID] = []
    hapID2nodes[hapID].append(node)

  for hapID, nodes_ in hapID2nodes.items():
    print("HapID", hapID, nodes_)

  # An extra round of computing of adjusted p-value
  hapPairs = set()
  hap2other = dict()
  for node, hapID in node2hapID.items():
    if hapID not in hap2other:
      loc = node.split(":")[0]
      locsNodes = list(loc2nodes[loc])
      h0 = node2hapID[locsNodes[0]]
      h1 = node2hapID[locsNodes[1]]
      hap2other[h0] = h1
      hap2other[h1] = h0
      hapPairs.add((h0, h1))
  hapID2nodes_new = dict()
  for h0, h1 in hapPairs:
    if len(hapID2nodes[h0])==1:
      continue
    hapID2nodes_new[h0] = []
    hapID2nodes_new[h1] = []
    h0_nodes = hapID2nodes[h0]
    h0_nodes.sort()
    h1_nodes = hapID2nodes[h1]
    h1_nodes.sort()
    assert len(h0_nodes)==len(h1_nodes), "No way!"
    test_num = [0]*len(h0_nodes)
    testPassed = [[] for x in range(0, len(h0_nodes))]
    for x in range(0, len(h0_nodes)-1):
      x1 = h0_nodes[x]
      x2 = h1_nodes[x]
      for y in range(x+1, len(h0_nodes)):
        y1 = h0_nodes[y]
        y2 = h1_nodes[y]
        table = [[adjacency[x1][y1] if (x1 in adjacency and y1 in adjacency[x1]) else 0,
                  adjacency[x1][y2] if (x1 in adjacency and y2 in adjacency[x1]) else 0],
                 [adjacency[x2][y1] if (x2 in adjacency and y1 in adjacency[x2]) else 0,
                  adjacency[x2][y2] if (x2 in adjacency and y2 in adjacency[x2]) else 0]]
        if table[0][0]+table[0][1]+table[1][0]+table[1][1] < args.min_read_num:
          continue
        res = fisher_exact(table, alternative='two-sided')
        test_num[x] += 1
        test_num[y] += 1
        if res[1] < 0.05 and x1 != y1 and x2 != y2 and table[0][0] > table[0][1] and table[1][0] < table[1][1]:
            testPassed[x].append(res[1])
            testPassed[y].append(res[1])
    passed = [False] * len(h0_nodes)
    test_num_total = sum(test_num)
    for x in range(0, len(h0_nodes)):
      pass_num = 0
      if test_num[x] > 0 and len(testPassed[x]) > 0:
        testAdjusted = [tmp*test_num[x] for tmp in testPassed[x]]
        # testAdjusted = [tmp*test_num_total for tmp in testPassed[x]]
        # testAdjusted = [tmp*(len(h0_nodes)-1) for tmp in testPassed[x]]
        # testAdjusted = [tmp for tmp in testPassed[x]]
        pass_num = 0
        for tmp in testAdjusted:
          if tmp < 0.05:
            pass_num += 1
        if pass_num >= test_num[x]/3:
          passed[x] = True
      if not passed[x]:
        print(h0_nodes[x], h1_nodes[x], "failed.", pass_num, test_num[x])
      else:
        hapID2nodes_new[h0].append(h0_nodes[x])
        hapID2nodes_new[h1].append(h1_nodes[x])
    print(sum(passed), "/", len(h0_nodes), "nodes of hap", (h0, h1), " passed adjusted Fisher's Exact test.")

  # Save hapGraph to file
  new_node_sets = set()
  for hapID, nodes_ in hapID2nodes_new.items():
    new_node_sets.update(nodes_)
  with open(args.output+gene+"_"+args.version+".hapgraph", 'w') as fout:
    fout.write("# Gene"+ gene+"\n")
    fout.write("[Nodes]\n")
    for node in hapgraph["Nodes"]:
      if str(node[0])+":"+node[1] in new_node_sets:
        fout.write(chr+":"+str(node[0])+","+node[1]+","+node[2]+"\n")
    fout.write("[Edges]\n")
    for edge in hapgraph["Edges"]:
      if (edge[0] in new_node_sets) and (edge[1] in new_node_sets):
        fout.write(chr+":"+edge[0]+","+chr+":"+edge[1]+","+str(edge[2])+"\n")

  # haplotypes = [[chr+":"+node for node in nodes_] for hapID,nodes_ in hapID2nodes.items()]
  haplotypes = [[chr+":"+node for node in nodes_] for hapID,nodes_ in hapID2nodes_new.items()]
  # for label in labels:
  #   if len(label) == 1:
  #     hap_ids.add(list(label)[0])
  # for hap_id in hap_ids:
  #   haplotype = []
  #   for x in range(0, len(nodes)):
  #     if len(labels[x]) == 1 and (hap_id in labels[x]):
  #       haplotype.append(chr+":"+nodes[x])
  #   if len(haplotype) > 1:
  #     haplotypes.append(haplotype)

  # G.add_node(1,pos=(1,1))
  # pos=nx.get_node_attributes(G, 'pos')
  # nx.draw(G, pos=pos, labels={1:"C"}, with_labels=True, node_color="whitesmoke")
  # plt.savefig("network.pdf")
  # sys.exit(0)
  hap_lens = []
  for hap in haplotypes:
    loci = [int(x.split(":")[1]) for x in hap]
    hap_lens.append(max(loci)-min(loci))
  hap_lens.sort(reverse=True)

  with open(args.output+"/"+gene+"_SNPsInWholeGene_"+args.version+".haps", 'w') as fout:
    for hap in haplotypes:
      fout.write(",".join(hap)+"\n")

  print("\t".join(["Gene", "Chr", "Start", "End", "Len"
    "alignedBaseInTotal", "DepthOfCoverage",
    # "#linkers", "#variants", "#ins", "#del", "#delFromRepeats", "#reads", "Ave.No.OfAllelesInReads",
    "hap_lens(#SNPs)", "hap_lens(bp)",
    "Max.Consi.Score", "Consi.Scores",
    "Class#1", "Class#2", "Class#3", "Class#4",
    "Hap.Lens", "Hap.Consi.Lens"
    ]))
  print("\t".join([str(x) for x in [gene, chr, start, end, end-start,
    alignedBaseInTotal, alignedBaseInTotal/(end-start)]+
    # len(linkers), VS_NumInTotal, VS_ins_NumInTotal, VS_del_NumInTotal, VS_del_fromRepeats_NumInTotal, read_num, stats.mean([len(x[1]) for x in readWithalles])]+
    [":".join([str(x) for x in [len(hap) in haplotypes]])] +
    [":".join([str(x) for x in hap_lens])] +
    compareWithGIAB(GIAB_SVs, haplotypes)]))
  haplotypes_SNPs = haplotypes

  # skip the following
  sys.exit(0)

  #*****  Whole-gene haplotype sequence assembly  *****
  #step 1: Cluster reads into haplotypes
  hapwithreads = {hapID:set() for hapID in hapID2nodes.keys()}
  hap_readratio_min = 2
  for read, alles in readWithalles:
    hap2count = dict()
    for alle in alles:
      node = ":".join([str(x) for x in alle])
      if node not in node2hapID:
        continue
      if node2hapID[node] not in hap2count:
        hap2count[node2hapID[node]] = 0
      hap2count[node2hapID[node]] += 1
    count2hap = [(x[1], x[0]) for x in hap2count.items()]
    count2hap.sort(reverse=True)
    if len(count2hap) == 0:
      # for hapID in hapID2nodes.keys():
      #   hapwithreads[hapID].add(read)
      continue
    else: # TODO: A read can be assigned to multiple haplotypes, depending on the propotion of agreement with haplotypes
      hapwithreads[count2hap[0][1]].add(read)
      # if count2hap[0][1] > 1:
      #   print(read.query_name, count2hap)
      # if len(count2hap) > 1 and count2hap[1][0]>0 and count2hap[0][0]/count2hap[1][0] < hap_readratio_min:
      #   if count2hap[1][1] not in hapwithreads:
      #     hapwithreads[count2hap[1][1]] = set()
      #   hapwithreads[count2hap[1][1]].add(read)
  print("HapID #supportingReads")
  for hap, reads in hapwithreads.items():
    print(hap, len(reads))
  #step 2: assemble each haplotypes
  hapID2seq = dict()
  hapID2cons = dict()
  for hapID, reads in hapwithreads.items():
    if len(reads) < args.min_cov or len(hapID2nodes[hapID]) < args.min_SNPs_inHap:
      continue
    freqTable = [{'A':0, 'T':0, 'G':0, 'C':0, 'In':[], 'De':0} for x in range(start, end)]
    #step 2.1: build table of all candidicate occurrences
    for read in reads:
      read_seq = read.query_sequence.upper()#query_alignment_sequence BUGs in pysam
      ref_loci = read.get_reference_positions(full_length=True)
      read_qual = read.query_qualities
      idx = 0
      while ref_loci[idx] == None or ref_loci[idx] < start:
        idx += 1
        continue
      while idx < len(ref_loci):
        loc = ref_loci[idx]
        if loc >= end:
          break
        freqTable[loc-start][read_seq[idx]] += 1
        next = idx + 1
        while next < len(ref_loci) and (ref_loci[next] == None or ref_loci[next] == loc): # insertion
          next += 1
        if next >= len(ref_loci):
          break
        if next - idx > 1: #confirm insertion
           freqTable[loc-start]['In'].append(read_seq[(idx+1):next])
        if ref_loci[next] <= end and ref_loci[next] > ref_loci[idx]+1: #confirm deletion
          for x in range(ref_loci[idx]+1, ref_loci[next]):
            freqTable[x-start]['De'] += 1
        idx = next
    # print(hapID, freqTable[89965487-start:89965495-start])
    #step 2.2 Generate consensus
    top_vs_in_max = 1.5
    cons = ['' for x in range(start, end)] #seq of each loc
    for idx in range(0, end-start):
      freqs = [(len(y) if x=="In" else y, x) for x,y in freqTable[idx].items()]
      freqs.sort(key=lambda x:x[0],reverse=True)
      # print(freqs)
      if freqs[0][0]==0:#werid, no coverage, but not classified as deletion
        print("check out this werid case", chr, start+idx, idx, freqs)
        cons[idx] = 'N'
      elif freqs[0][0] < args.min_cov:
        cons[idx] = 'N'
        continue
      elif freqs[0][1] == 'In' or (freqs[1][1] == 'In' and freqs[1][0] >= args.min_cov and freqs[0][0]/freqs[1][0] <= top_vs_in_max):
        if False: # use MSA to compile the consensus of insertion
          # print("insert to be resolved", chr, start+idx, idx, freqs)
          # print(freqTable[idx]["In"])
          tmp_filename=args.output+"/tmp.fa"
          with open(tmp_filename, 'w') as fout:
            x = 0
            for seq in freqTable[idx]["In"]:
              fout.write(">"+str(x)+"\n")
              fout.write(seq+"\n")
              x += 0
          # mafft="/public/hshi/tools/miniconda2/envs/condaPython3/bin/mafft"
          # results = subprocess.check_output([mafft, "--auto", tmp_filename])
          results = subprocess.run([args.mafft, "--auto", tmp_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE)#, text=True)
          # results = result.stdout.split("\r")
          results = [line for line in results.stdout.decode('utf-8').split("\n") if line.strip()]
          MSAs = []
          x = 0
          while x < len(results):
            if results[x] and results[x][0] != '>':
              x += 1
              continue
            x += 1
            msa = ""
            while x < len(results) and results[x] and results[x][0] != '>':
              msa = msa + results[x]
              x += 1
            MSAs.append(msa.upper())
          msa_len = len(MSAs[0])
          DNAbases = {'A', 'T', 'G', 'C'}
          base2cov = [{base:0 for base in DNAbases} for x in range(0, msa_len)]
          for msa in MSAs:
            for x in range(0, msa_len):
              if msa[x] in DNAbases:
                base2cov[x][msa[x]] = base2cov[x][msa[x]] + 1
          # print(base2cov)
          # base2check = []
          baseWithGoodCoverage = []
          baseWithGoodCoverage_num = 0
          for covs in base2cov:
            covs = [[base, freq] for base,freq in covs.items()]
            covs.sort(key=lambda x:x[1], reverse=True)
            if covs[0][1] >= args.min_cov and covs[0][1]/len(MSAs) >= args.min_ratio: #or (len(covs) >= 2 and (covs[0][1]+covs[1][1]) >= min_cov):
              baseWithGoodCoverage.append(covs[0][0])
              # base2check.append(True)
              baseWithGoodCoverage_num = baseWithGoodCoverage_num + 1
            else:
              baseWithGoodCoverage.append(" ")
              # base2check.append(False)
          # print(base2cov)
          # print(str(baseWithGoodCoverage_num)+"/"+str(msa_len)+" bases with coverage >= "+str(args.min_cov))
          cons_tmp = ""
          for base in baseWithGoodCoverage:
            if base == " ":
              continue
            cons_tmp = cons_tmp + base
        elif True: # select the most frequent insertion sequences
          seq2count = dict()
          for seq_tmp in freqTable[idx]["In"]:
            if seq_tmp not in seq2count:
              seq2count[seq_tmp] = 0
            seq2count[seq_tmp] += 1
          count2seq = [[count, seq_tmp] for seq_tmp, count in seq2count.items()]
          count2seq.sort(reverse=True)
          cons_tmp = count2seq[0][1]
        extra=''
        if freqs[0][1] in BASES:
          extra=freqs[0][1]
        cons[idx]=extra + cons_tmp
      elif freqs[0][1] == "De": # stringent criteria for deletion
        if freqs[1][0] >= args.min_cov and (freqs[1][1] not in ["De", "In"]) and freqs[0][0]/freqs[1][0] <= top_vs_in_max:
          cons[idx] = freqs[1][1]
        else:
          cons[idx] = ""
      else:
        # print("Hi")
        cons[idx] = freqs[0][1]
      # print(cons[idx])
      # sys.exit(0)
    hapID2cons[hapID] = cons
    hapseq = "".join(cons)
    hapID2seq[hapID] = hapseq
    # print(cons[0:10])
  with open(args.output+"/"+gene+"_hapseqs"+"_"+args.version+".fa", 'w') as fout:
    for hapID, seq in hapID2seq.items():
      fout.write(">"+str(hapID)+"\n")
      fout.write(seq+"\n")
  # compile a list of new variants based on the whole-gene consensus
  hetVars = dict()
  hetVar_total = 0
  hetVar_inGIAB = 0
  hetVar_inGIAB_consAllele = 0
  hetVar_morethan2cand = 0
  for idx in range(0, end-start):
    vars = list(set([cons[idx] for hapID, cons in hapID2cons.items()]))
    if len(vars) == 2:# heterozygous variants
      if ('N' in vars) or ('N' in vars[0]) or ('N' in vars[1]):
        continue
      var2hapIDs = dict()
      for hapID, cons in hapID2cons.items():
        if cons[idx] not in var2hapIDs:
          var2hapIDs[cons[idx]] = set()
        var2hapIDs[cons[idx]].add(hapID)
      hetVars[start+idx] = var2hapIDs
    elif len(vars) > 2:
      hetVar_morethan2cand += 1
  hetVars = [[pos, vars] for pos,vars in hetVars.items()]
  hetVars.sort(key=lambda x: x[0])
  hetVars_merged = [] #merge continuous deletions
  idx = 0
  while idx < len(hetVars):
    if "" in hetVars[idx][1]:
      pos = hetVars[idx][0]
      vars = hetVars[idx][1]
      del_hapID = hetVars[idx][1][""]
      seq_hap= "".join(list(hetVars[idx][1].keys()))
      seq_hapID = hetVars[idx][1][seq_hap]
      next = idx+1
      while next < len(hetVars) and ("" in hetVars[next][1]):
        if hetVars[next-1][0]+1 == hetVars[next][0]:
          if hetVars[next][1][""] == del_hapID and hetVars[next][1]["".join(list(hetVars[next][1].keys()))] == seq_hapID:
            seq_hap = seq_hap + "".join(list(hetVars[next][1].keys()))
            next += 1
          else:
            break
        else:
          break
      hetVars_merged.append([pos, {"":del_hapID, seq_hap:seq_hapID}])
      idx = next
    else:
      hetVars_merged.append(hetVars[idx])
      idx += 1
  # Filter out variants nearby repetitive regions
  log("Filter out variants nearby repetitive regions.")
  if args.noRepeat:
    hetVars_merged_filter = []
    for ref_loc, vars in hetVars_merged:
      inRepeat = False
      if False not in [len(x)==1 for x in list(vars.keys())]:# SNPs
        # Check one-base repeats
        if (len(set([x for x in chr2seq[chr][max(0, ref_loc-args.repeat_len):ref_loc]])) == 1) or \
          (len(set([x for x in chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len+1, len(chr2seq[chr]))]])) == 1):
          print([chr2seq[chr][max(0, ref_loc-args.repeat_len):ref_loc], chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len+1, len(chr2seq[chr]))]])
          continue
        # Check two-bases repeats
        for seq_tmp in [chr2seq[chr][max(0, ref_loc-args.repeat_len*2):ref_loc], chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len*2+1, len(chr2seq[chr]))]]:
          patterns = set()
          for x in range(0, len(seq_tmp)-1, 2):
            patterns.add(seq_tmp[x:(x+2)])
          if len(patterns) == 1:
            inRepeat = True
            break
        if inRepeat:
          print([chr2seq[chr][max(0, ref_loc-args.repeat_len*2):ref_loc], chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len*2+1, len(chr2seq[chr]))]])
          continue
      else:# indels
        # Check one-base repeats
        for seq_tmp in [chr2seq[chr][max(0, ref_loc-args.repeat_len):ref_loc], chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len+1, len(chr2seq[chr]))]]:
          rep2count = list(Counter(seq_tmp).items())
          rep2count.sort(key = lambda x : x[1], reverse=True)
          rep_ratio = rep2count[0][1]/float(len(seq_tmp))
          if rep_ratio >= args.repeat_ratio:
            print([chr2seq[chr][max(0, ref_loc-args.repeat_len):ref_loc], chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len+1, len(chr2seq[chr]))]])
            inRepeat = True
            break
        if inRepeat:
          continue
        # Check two-bases repeats
        for seq_tmp in [chr2seq[chr][max(0, ref_loc-args.repeat_len*2):ref_loc], chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len*2+1, len(chr2seq[chr]))]]:
          rep2count = dict()
          for x in range(0, len(seq_tmp)-1):
            tmp = seq_tmp[x:(x+2)]
            if tmp not in rep2count:
              rep2count[tmp] = 0
            rep2count[tmp] += 1
          rep2count = [[key, value] for key, value in rep2count.items()]
          rep2count.sort(key = lambda x : x[1], reverse=True)
          rep_ratio = rep2count[0][1]*2/float(len(seq_tmp))
          if rep_ratio >= args.repeat_ratio:
            print([chr2seq[chr][max(0, ref_loc-args.repeat_len*2):ref_loc], chr2seq[chr][(ref_loc+1):min(ref_loc+args.repeat_len*2+1, len(chr2seq[chr]))]])
            inRepeat = True
            break
      if inRepeat:
        continue
      hetVars_merged_filter.append([ref_loc, vars])
    print("[Info] filter out", len(hetVars_merged)-len(hetVars_merged_filter), "variants nearby repetitive regions.", len(hetVars_merged_filter), "variances remains after filtering.")
    hetVars_merged = hetVars_merged_filter

  # Count the number of variants
  hetVars_counts = {"SNP":[], "Ins":[], "Del":[]}
  for pos, vars in hetVars_merged:
    alleles = list(vars.keys())
    if len(alleles[0]) == 1 and len(alleles[1]) == 1:
      loc_found=False
      alleles_found=False
      if chr+":"+str(pos) in GIAB_SVs:
        loc_found = True
        if (alleles[0] in GIAB_SVs[chr+":"+str(pos)]["alleles"]) and (alleles[1] in GIAB_SVs[chr+":"+str(pos)]["alleles"]):
          alleles_found=True
      hetVars_counts["SNP"].append([pos, vars, loc_found, alleles_found])
    elif len(alleles[0]) >= 1 and len(alleles[1]) >= 1:
      loc_found=False
      alleles_found=False
      if chr+":"+str(pos) in GIAB_SVs:
        loc_found = True
        if (alleles[0] in GIAB_SVs[chr+":"+str(pos)]["alleles"]) and (alleles[1] in GIAB_SVs[chr+":"+str(pos)]["alleles"]):
          alleles_found=True
      hetVars_counts["Ins"].append([pos, vars, loc_found, alleles_found])
    else:
      loc_found=False
      alleles_found=False
      if chr+":"+str(pos-1) in GIAB_SVs:
        loc_found = True
        if (hapID2cons[list(vars[alleles[0]])[0]][pos-1-start]+alleles[0] in GIAB_SVs[chr+":"+str(pos-1)]["alleles"]) and (hapID2cons[list(vars[alleles[1]])[0]][pos-1-start]+alleles[1] in GIAB_SVs[chr+":"+str(pos-1)]["alleles"]):
          alleles_found=True
      # if loc_found and (not alleles_found):
      #   print(pos, vars)
      hetVars_counts["Del"].append([pos, vars, loc_found, alleles_found])
  for type, vars in hetVars_counts.items():
    loc_found = 0
    alleles_found = 0
    for var in vars:
      if var[2]:
        loc_found += 1
      if var[3]:
        alleles_found += 1
    print(type, len(vars), loc_found, alleles_found)
  # print(hetVars_counts["Ins"])
  # Resolve haplotypes
  hapIDs = set()
  for loc,vars in hetVars_merged:
    for hapID in vars.values():
      hapIDs.add(tuple(hapID))
  haplotypes = dict()
  for loc,vars in hetVars_merged:
    isIndel = False
    for seq in vars.keys():
      if seq == "":
        isIndel = True
        break
    for seq,hapID in vars.items():
      hapID = tuple(hapID)
      if hapID not in haplotypes:
        haplotypes[hapID] = []
      if isIndel:
        seq = hapID2cons[hapID[0]][loc-1-start]+seq
        haplotypes[hapID].append(chr+":"+str(loc-1)+":"+seq)
      else:
        haplotypes[hapID].append(chr+":"+str(loc)+":"+seq)
  hap2nodes_stage2 = {key:value for key,value in haplotypes.items()}
  haplotypes = list(haplotypes.values())

  with open(args.output+"/"+gene+"_VarsInWholeGene_"+args.version+".haps", 'w') as fout:
    for hap in haplotypes:
      fout.write(",".join(hap)+"\n")

  # Cross-check with the SNPs only haplotypes
  nodes_tmp = set()
  for hap in haplotypes:
    nodes_tmp.update(hap)
  node_common = 0
  for hap in haplotypes_SNPs:
    for node in hap:
      if node in nodes_tmp:
        node_common += 1
      else:
        print(node)
  print(node_common, "SNPs also detected after whole genome assembly.")
  print(compareWithGIAB(GIAB_SVs, haplotypes))
  # print(hap2nodes_stage2)
  # print(hapID2nodes)
  # Merge the stage 1 and stage 2 variant calls
  hap2nodes_final = dict()
  for hapID, nodes in hapID2nodes.items():
    if len(nodes) == 1:
      continue
    pos2alleles = dict()
    for node in nodes:
      eles = node.split(":")
      pos2alleles[chr+":"+eles[0]] = eles[1]
    inconsistent_alleles = 0
    # for key,value in hap2nodes_stage2.items():
    #   print(key,value)
    # print(type(hap2nodes_stage2))
    # print(tuple(hapID), hap2nodes_stage2.keys())
    # if (hapID,) not in hap2nodes_stage2:
    #   print(tuple(hapID))
    for node in hap2nodes_stage2[(hapID,)]:
      eles = node.split(":")
      pos_tmp = eles[0]+":"+eles[1]
      if pos_tmp in pos2alleles:
        if eles[2] != pos2alleles[pos_tmp]:
          pos2alleles[pos_tmp] = eles[2]
          inconsistent_alleles += 1
      else:
        pos2alleles[pos_tmp] = eles[2]
    print("For hap.", hapID, inconsistent_alleles, "alleles are inconsistent between first and second stage.")
    hap2nodes_final[hapID] = [pos_tmp+":"+allele for pos_tmp, allele in pos2alleles.items()]
  hap_final = hap2nodes_final.values()
  # print("After merging, we have ", len(hap_final[0]), "")
  for hapID, nodes in hap2nodes_final.items():
    print(hapID, nodes)
  print(compareWithGIAB(GIAB_SVs, hap_final))
