# # # # # # # # # # # # # # # # # # # # # #
# Author: Jun He
# Date: 07/04/2016
#
# This is a Python script to annotate variants type according to the position and allel. It supports two models of annotation, one is tree-based while the other is array-based.
# The tree-based model is efficient for random position annotation and small scale annotation. And the array-based model is efficient for sorted large scale annotation.
# Make sure your vcf file is sorted when you are using array model. If you don't want to sort your vcf file, try tree model.
# 
# This script can be easily extended to support an efficient online annotation tool, if you are interested.
#
# This script is written with the concept of Object Oriented Design.
#
# Usage:
#   e.g. input for annotation:
#   
#   python annotate_variant.py gene_array.pkl mutect.sorted.vcf output.vcf
#    gene_array.pkl / gene_tree.pkl     :  An array or tree model file that is built by this script. If your vcf file is unsorted, use tree model.
#    mutect.sorted.vcf                  :  A vcf file to be annotated.
#    output.vcf                         :  The annotated vcf file
#
#   e.g. input for building and saving a new model:
#
#   python annotate_variant.py -build_tree/-build_array dogGene.gtf genome.fa mutect.vcf gene_name.tsv output.vcf
#
#    -build_tree/-build_array  :  -build_tree, command to build a tree-based model
#                                 -build_array, command to build an array-based model
#
#    dogGene.gtf    :       A GTF file describing transcripts.  Transcripts must be grouped by Gene ID, i.e., all transcripts from the same Gene ID must be in a consecutive
#                           block of lines. Gene IDs and transcript IDs must be of the same type as the gene and transcript identifiers in columns 1 and 2 of "gene_name_file"
#    
#    genome_file    :       A genome sequence file (in multipart FASTA format); sequence names must be string-identical to the chromosome identifiers in the VCF file.
#    
#    
#    gene_name_file :       A tab-delimited text file with the Ensembl transcript id in the 1st column, Ensembl gene id in the 2nd column and gene name on the 3rd 
#                           column. No quotes. You can obtain this file from Ensembl BioMart.
#    
#    output_file    :       The output file for storing the built model.
#
#
# References:   The Variant Call Format (VCF) Version 4.2 Specification, 26 January 2015
#
# This script was written in Python 2.7
#
# # # # # # # # # # # # # # # # # # # # # #

import sys
import cPickle as pickle

variant_interpreter = {
    'intron'    : 'intron_variant|MODIFIER',
    'intergenic': 'intergenic_region|MODIFIER',
    '5UTR'      : '5_prime_UTR_variant|MODIFIER',
    '3UTR'      : '3_prime_UTR_variant|MODIFIER',
    'non_coding': 'non_coding_exon_variant|MODIFIER',
    'synonymous': 'synonymous_variant|LOW',
}

cmpl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}

# Global function to translate the DNA sequence into its opposite direction.
def complement(s):
    s = list(s)
    for i in xrange(len(s)):
	s[i] = cmpl[s[i]]
    return ''.join(s)

amino_codon = {    
    'ATT':     'I',
    'ATC':     'I',
    'ATA':     'I',
    'CTT':     'L',
    'CTC':     'L',
    'CTA':     'L',
    'CTG':     'L',
    'TTA':     'L',
    'TTG':     'L',
    'GTT':     'V',
    'GTC':     'V',
    'GTA':     'V',
    'GTG':     'V',
    'TTT':     'F',
    'TTC':     'F',
    'ATG':     'M',
    'TGT':     'C',
    'TGC':     'C',
    'GCT':     'A',
    'GCC':     'A',
    'GCA':     'A',
    'GCG':     'A',
    'GGT':     'G',
    'GGC':     'G',
    'GGA':     'G',
    'GGG':     'G',
    'CCT':     'P',
    'CCC':     'P',
    'CCA':     'P',
    'CCG':     'P',
    'ACT':     'T',
    'ACC':     'T',
    'ACA':     'T',
    'ACG':     'T',
    'TCT':     'S',
    'TCC':     'S',
    'TCA':     'S',
    'TCG':     'S',
    'AGT':     'S',
    'AGC':     'S',
    'TAT':     'Y',
    'TAC':     'Y',
    'TGG':     'W',
    'CAA':     'Q',
    'CAG':     'Q',
    'AAT':     'N',
    'AAC':     'N',
    'CAT':     'H',
    'CAC':     'H',
    'GAA':     'E',
    'GAG':     'E',
    'GAT':     'D',
    'GAC':     'D',
    'AAA':     'K',
    'AAG':     'K',
    'CGT':     'R',
    'CGC':     'R',
    'CGA':     'R',
    'CGG':     'R',
    'AGA':     'R',
    'AGG':     'R',
    'TAA':     '*',
    'TAG':     '*',
    'TGA':     '*',    
    }

class interval(object):
    def __init__(self, lo=sys.maxint, hi=0):
	self.lo = lo
	self.hi = hi
    
    def set(self, lo, hi):
	self.lo = lo
	self.hi = hi
    
    def overlap(self, pos):
	return self.lo <= pos and self.hi >= pos

# exon object stores the range of an exon
class exon(object):
    def __init__(self, lo, hi):
	self.range = interval(lo, hi)
	self.offset = 0
    
    def lo(self):
	return self.range.lo
    
    def hi(self):
	return self.range.hi
    
    def overlap(self, pos):
	return self.range.overlap(pos)

# gene object stores all information about a gene
class gene(object):
    def __init__(self):
	self.range = interval()
        self.exons = []
	self.seq = ''
        self.start = interval()
        self.stop = interval()
        self.direction = '+'
	self.transcript_id = ''
	self.ensembl_id = ''
	self.gene_name = ''
    
    # API for non-intergenic annotation. Make sure the position is inside the gene.
    def get_variant_type_at(self, pos, nuc):
	if not self.overlap(pos):
	    return variant_interpreter['intergenic']
	if self.start.lo > self.start.hi:
	    return variant_interpreter['non_coding']
	
	for i, exon in enumerate(self.exons):
	    if exon.overlap(pos):
		if self.is_5_UTR(pos):
		    return variant_interpreter['5UTR']
		if self.is_3_UTR(pos):
		    return variant_interpreter['3UTR']	
		
		return self.get_amino_variant_at(exon, pos, nuc)
	
	return variant_interpreter['intron']
    
    def lo(self):
	return self.range.lo
    def hi(self):
	return self.range.hi
    
    def set_range(self):
	self.range.set(self.exons[0].lo(), self.exons[-1].hi())
    
    def overlap(self, pos):
	return self.range.overlap(pos)
    
    # Make sure pos overlaps an exon before using this function
    def is_5_UTR(self, pos):
	return (self.direction == '+' and pos < self.start.lo) or (self.direction == '-' and pos > self.start.hi)
    
    # Make sure pos overlaps an exon before using this function
    def is_3_UTR(self, pos):
	return (self.direction == '+' and pos >= self.stop.lo) or (self.direction == '-' and pos <= self.stop.hi)
	
    def get_amino_variant_at(self, exon, pos, nuc):
	if self.direction == '+':
	    offset = exon.offset + pos - max(exon.lo(), self.start.lo)
	else:
	    offset = exon.offset + min(exon.hi(), self.start.hi) - pos
	    nuc = complement(nuc)
	    
	origin_amino = self.seq[offset/3*3: offset/3*3+3]
	t = list(origin_amino)
	t[offset%3] = nuc
	variant_amino = ''.join(t)
	
	if 'N' in origin_amino:
	    return 'non_sequencing_position'
	
	origin_amino = amino_codon[origin_amino]
	variant_amino = amino_codon[variant_amino]
	
	if origin_amino == variant_amino:
	    return variant_interpreter.get('synonymous')
	
	return'p.' + origin_amino + str(offset / 3 + 1) + variant_amino + '|MODERATE'
    

# Tree Node for gene_tree
class gene_tree_node(object):
    def __init__(self, gene):
	self.gene = gene
	self.left = None
	self.right = None
	self.mini = gene.lo()
	self.maxi = gene.hi()
    
    def overlap(self, pos):
	return self.gene.overlap(pos)
    
    def contain(self, pos):
	return self.mini <= pos and self.maxi >= pos

# An interval tree that stores gene objects in a chromosome. One tree per chromosome.
# Efficient for random or small scale annotation
class gene_tree(object):
    def __init__(self, ar):
	self.root = None
	self.build_balanced_bst(ar)
	
	
    def build_balanced_bst(self, gene_seq):
	self.root = self.balanced_insert(gene_seq, 0, len(gene_seq))
    
    def balanced_insert(self, gene_seq, l, r):
	if l>=r:
	    return None
	
	mid = (l+r)/2
	root = gene_tree_node(gene_seq[mid])
	
	a1, b1, a2, b2 = 0, sys.maxint, 0, sys.maxint
	
	if l < mid:
	    root.left = self.balanced_insert(gene_seq, l, mid)
	    a1, b1 = root.left.maxi, root.left.mini
	if mid+1 < r:
	    root.right = self.balanced_insert(gene_seq, mid+1, r)
	    a2, b2 = root.right.maxi, root.right.mini

	root.maxi = max(root.maxi, a1, a2)
	root.mini = min(root.mini, b1, b2)
	return root
    
    def get_genes_at(self, pos):
	genes = []
	self.query_all_at(self.root, pos, genes)
	return genes
    
    def query_all_at(self, root, pos, genes):
	if root is None:
	    return
	if not root.contain(pos):
	    return
	
	if root.overlap(pos):
	    genes.append(root.gene)
	
	self.query_all_at(root.left, pos, genes)
	self.query_all_at(root.right, pos, genes)
	
# An sorted array that stores gene objects in a chromosome. One array per chromosome.
# Efficient for sorted large scale annotation
class gene_array(object):
    def __init__(self, ar):
	self.genes = ar
	self.i = 0
    
    def reset(self):
	self.i = 0
    
    def get_genes_at(self, pos):
	ret = []
	while self.i < len(self.genes) and self.genes[self.i].lo() <= pos:
	    if self.genes[self.i].overlap(pos):
		break
	    self.i += 1
	
	j = self.i
	while j<len(self.genes) and self.genes[j].lo() <= pos:
	    if self.genes[j].overlap(pos):
		ret.append(self.genes[j])
	    j += 1
	
	return ret
	
# Read DNA sequence for each gene and exon from a genome file.
class genome(object):
    def __init__(self, filename):
        self.f = open(filename, 'r')
        self.chr = None
        self.next_chr = self.f.readline()[1:-1]
        self.chr_seq = ''
    
    def __del__(self):
        self.f.close()
    
    def read_chr(self, chr):
        if self.chr == chr:
            return
        
        if self.next_chr != chr:
            for line in self.f:
                if line[0] == '>' and line[1:] == chr:
                    self.next_chr = line[1:-1]
                    break
        
        buffer = []
        for line in self.f:
            
            if line[0] == '>':
                
                self.next_chr = line[1:-1]
                break
            
            buffer.append(line[:-1])
        self.chr = chr
        self.chr_seq = ''.join(buffer)

    def get_gene_seq(self, chr, itv, direction):
        self.read_chr(chr)
        
        if direction == '+':
	    seq = self.chr_seq[itv.lo-1:itv.hi]
	else:
	    seq = self.chr_seq[itv.lo-1:itv.hi][::-1]
	    seq = complement(seq)
	return seq
    

# core annnotation module
# superclass for tree_annotation and array_annotation
class annotation(object):
    def __init__(self):
	self.chromo = None
	
    def build_chromosomes(self, gtf, fa, gn=None):
	raise NotImplementedError("Subclass must implement abstract method")
    
    def build_genes(self, gtf, fa, gn=None):
	self.chromo = {}
	DNA = genome(fa)
	gname = None
	if gn is not None:
	    gname = self.gn_read(gn)	
	chromosomes = self.gtf_read(gtf, DNA, gname)
		

	for chr, genes in chromosomes.iteritems():
	    yield chr, sorted(genes.itervalues(), key = lambda k: k.exons[0].lo())
	    
    
    def save(self, filename):
	f = open(filename, 'wb')
	pickle.dump(self.chromo, f)
	f.close()
	
    def load(self, filename):
	f = open(filename, 'rb')
	self.chromo = pickle.load(f)
	f.close()
	
    
    def gene_insert(self, dict, chr, gene_id, d):
	if (chr not in dict):
	    dict[chr] = {}
    
	if (gene_id not in dict[chr]):
	    dict[chr][gene_id] = gene()
	
	dict[chr][gene_id].direction = d
	dict[chr][gene_id].transcript_id = gene_id
	    
    def update_offset(self, chr, gene, DNA):
	if gene == '':
	    return
	
	itv = interval()
	gene.set_range()
	
	if gene.direction == '+':
	    for i in xrange(0, len(gene.exons)):
		if i>0:
		    l = max(0, gene.exons[i-1].hi() + 1 - max(gene.start.lo, gene.exons[i-1].lo()))
		    gene.exons[i].offset = gene.exons[i-1].offset + l
		    
		itv.set(max(gene.start.lo, gene.exons[i].lo()) , gene.exons[i].hi())
		gene.seq += DNA.get_gene_seq(chr, itv, gene.direction)		
	else:
	    for i in xrange(len(gene.exons)-1, -1, -1):
		if i < len(gene.exons)-1:
		    l = max(0, min(gene.start.hi, gene.exons[i+1].hi()) + 1 - gene.exons[i+1].lo())
		    gene.exons[i].offset = gene.exons[i+1].offset + l
		
		itv.set(gene.exons[i].lo(), min(gene.start.hi, gene.exons[i].hi()) )
		gene.seq += DNA.get_gene_seq(chr, itv, gene.direction)
    
    def update_gene_name(self, gene, gname=None):
	if gname is not None:
	    p = gname.get(gene.transcript_id, ['',''])
	    gene.ensembl_id, gene.gene_name = p[0], p[1]
	
    def gtf_read(self, gtf, DNA, gname=None):
	f = open(gtf)
	
	ret = {}
	last_gene = ''
	last_chr = ''
	
	for line in f:
	    p = line[:-1].split('\t');
	    if p[2] == 'CDS':
		continue
	    
	    chr = p[0][3:]
	    if chr == 'M':
		chr = 'MT'
	    gene_id = p[8][9:27]
	    if last_gene != gene_id:
		if last_gene != '':
		    self.update_offset(last_chr, ret[last_chr][last_gene], DNA)
		    self.update_gene_name(ret[last_chr][last_gene], gname)
		last_gene = gene_id
		last_chr = chr
		self.gene_insert(ret, chr, gene_id, p[6])
		
	    
	    lo = long(p[3])
	    hi = long(p[4])

	    if 'exon' == p[2]:
		ret[chr][gene_id].exons.append(exon(lo, hi))
	    
	    if 'start_codon' == p[2]:
		ret[chr][gene_id].start.set(lo, hi)
    
	    if 'stop_codon' == p[2]:
		ret[chr][gene_id].stop.set(lo, hi)
	
	self.update_offset(chr, ret[chr][gene_id], DNA)
	self.update_gene_name(ret[chr][gene_id], gname)
	f.close()
	return ret    
    
    def gn_read(self, gname):
	f = open(gname)
	ret = {}
	for line in f:
	    p = line[:-1].split('\t')
	    ret[p[0]] = [p[1], p[2]]
	f.close()
	return ret
    
    def get_genes_at(self, chr, pos):
	return self.chromo[chr].get_genes_at(pos)
    
    def build_annotation(self, nuc, t, gene):
	if gene is None:
	    return nuc + '|' + t + '|||'
	return nuc + '|' + t + '|' + gene.transcript_id + '|' + gene.ensembl_id + '|' + gene.gene_name  
    
    def annotate(self, chr, pos, nuc):
	genes = self.get_genes_at(chr, pos)
	if len(genes) == 0:
	    return [ self.build_annotation(nuc, variant_interpreter['intergenic'], None) ]
	
	ret = []
	for gene in genes:
	    ret.append( self.build_annotation(nuc, gene.get_variant_type_at(pos, nuc),  gene) )
	
	return ret	
	    
    def parse(self, string):
	chr, pos, nuc = string.split()
	
	return self.annotate(chr, long(pos), nuc)
    
    def annotate_vcf(self, vcf, output):
	f = open(vcf)
	fo = open(output, 'w')
	
	insert_info = False
	count = 0
	
	for line in f:
	    if line[0] == '#': 
		if insert_info and line[2:6] != 'INFO':
		    fo.write('##INFO=<ID=ANN,Number=1,Type=String,Description="variant annotation">\n')
		    fo.write('##INFO=<ID=VSF,Number=1,Type=String,Description="variant source file">\n')
		    insert_info = False
		elif line[2:6] == 'INFO':
		    insert_info = True
		fo.write(line)
		continue
	    
	    p = line.split()
	    
	    chr = p[0]
	    pos = long(p[1])
	    nuc = p[4]	    
	    
	    annotation = self.annotate(chr, pos, nuc)
	    info = 'ANN=' + ','.join(annotation)
	    p[7] = p[7] + ';' + info + ';' + 'VSF=' + vcf
	    
	    fo.write('\t'.join(p) + '\n')
	    
	    count = count + 1
	    if count & 15 == 0:
		print ('\rAnnotated: ' + str(count)),
		
	print '\nCompleted! Annotated '+ str(count) + ' variants.'
	
	f.close()
	fo.close()
# Search genes on balanced binary search tree structure.
# Support both random annotation and sorted annotation
# Efficient for random or small scale annotation
class tree_annotation(annotation):
    
    # API for build a new gene balanced binary search tree model
    def build_chromosomes(self, gtf, fa, gn=None):
	for chr, genes in self.build_genes(gtf, fa, gn):
	    self.chromo[chr] = gene_tree(genes)
    
    # API for save the current model to a file
    def save(self, filename):
	super(array_annotation, self).save(filename)  
    
    # API for load a model from a file
    def load(self, filename):
	super(array_annotation, self).load(filename)  
    
    # API for VCF file annotation
    def annotate_vcf(self, vcf, output):
	super(array_annotation, self).annotate_vcf(vcf, output)  
    
    # API for one random annotation
    # the query string in the form: chromosome pos allele
    def parse(self, string):
	super(array_annotation, self).parse(string)  


# Search genes on sorted array structure.
# Support only sorted VCF file annotation
# Efficient for sorted large scale annotation
class array_annotation(annotation):
    
    # API for build a new gene array model
    def build_chromosomes(self, gtf, fa, gn=None):
	for chr, genes in self.build_genes(gtf, fa, gn):
	    self.chromo[chr] = gene_array(genes)    

    # API for save the current model to a file
    def save(self, filename):
	super(array_annotation, self).save(filename)  
    
    # API for load a model from a file
    def load(self, filename):
	super(array_annotation, self).load(filename)  
    
    # API for sorted VCF file annotation
    def annotate_vcf(self, vcf, output):
	self.reset()
	super(array_annotation, self).annotate_vcf(vcf, output)   
	
    def reset(self):
	for gene_array in self.chromo.itervalues():
	    gene_array.reset()

usage = ['Usage:',
  'e.g. input for annotation:',
  '',
  '   python annotate_variant.py gene_array.pkl mutect.sorted.vcf output.vcf',
  '',
  '    gene_array.pkl / gene_tree.pkl     :  An array or tree model file that is built by this script. If your vcf file is unsorted, use tree model.',
  '    mutect.sorted.vcf                  :  A vcf file to be annotated.',
  '    output.vcf                         :  The annotated vcf file',
  '',
  'e.g. input for building and saving a new model:',
  '',
  '   python annotate_variant.py -build_tree/-build_array dogGene.gtf genome.fa mutect.vcf gene_name.tsv output.vcf',
  '',
  '    -build_tree/-build_array  :  -build_tree, command to build a tree-based model',
  '                                 -build_array, command to build an array-based model',
  '',
  '    dogGene.gtf    :       A GTF file describing transcripts.  Transcripts must be grouped by Gene ID, i.e., all transcripts from the same Gene ID must be in a consecutive',
  '                           block of lines. Gene IDs and transcript IDs must be of the same type as the gene and transcript identifiers in columns 1 and 2 of "gene_name_file"',
  '    ',
  '    genome_file    :       A genome sequence file (in multipart FASTA format); sequence names must be string-identical to the chromosome identifiers in the VCF file.',
  '    ',
  '    gene_name_file :       A tab-delimited text file with the Ensembl transcript id in the 1st column, Ensembl gene id in the 2nd column and gene name on the 3rd ',
  '                           column. No quotes. You can obtain this file from Ensembl BioMart.',
  '    ',
  '    output_file    :       The output file for storing the built model.',
  ]

if __name__=='__main__':
    
    if len(sys.argv) != 4 and len(sys.argv) != 6:
	print '\n'.join(usage)
	exit()
    

    if '-build_tree' in sys.argv[1]:
	tree_model = tree_annotation()
	tree_model.build_chromosomes(sys.argv[2], sys.argv[3], sys.argv[4])
	tree_model.save(sys.argv[5])
    elif '-build_array' in sys.argv[1]:
	array_model = array_annotation()
	array_model.build_chromosomes(sys.argv[2], sys.argv[3], sys.argv[4])
	array_model.save(sys.argv[5])	
    else:
	model = annotation()
	model.load(sys.argv[1])
	model.annotate_vcf(sys.argv[2], sys.argv[3])

	


