# # # # # # # # # # # # # # # # # # # # # #
# Author: Jun He
#
# Usage:
#    gtf_file       :       A GTF file sorted first by chromosome then by gene_id then by coordinate.
#    genome_file    :       A genome sequence file.
#    vcf_file       :       The VCF file to be annotated. Must be sorted first by chromosome then by coordinate.
#    gene_name_file :       A file with transcript id on the 1st column, gene id on the 2nd column and gene name on the 3rd 
#                           column.
#    output_file    :       A designated output file.
#
# Note: When the number of input arguments is not correct, an usage message will be printed on the screen.
#
# # # # # # # # # # # # # # # # # # # # # #
import sys
from subprocess import Popen, PIPE

variant_interpreter = {
        1: 'intron_variant|MODIFIER',
        2: 'intergenic_region|MODIFIER',
        3: '5_prime_UTR_variant|MODIFIER',
        4: '3_prime_UTR_variant|MODIFIER',
        5: 'non_coding_exon_variant|MODIFIER',
        6: 'synonymous_variant|LOW',
    }

gene_name = {}

cmpl = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

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

class gene(object):
    def __init__(self, range=None, start=None, stop=None, direction=None):
        self.range = []
        self.start = start
        self.stop = stop
        self.direction = direction

class genome(object):
    def __init__(self, filename):
        self.f = open(filename, 'r')
        self.chr = None
        self.next_chr = self.f.readline()[1:-1]
        self.chr_seq = ''
        self.gene_seq = ''
    
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
    
    def test(self, variant_pos, s, e, d):
        if d == '+':
            return self.chr_seq[variant_pos + s-1: variant_pos + e-1]
        
        if s == 0:
            seq = self.chr_seq[variant_pos -e: variant_pos]
        elif e ==0:
            seq = self.chr_seq[variant_pos -1: variant_pos - s]
        else:
            seq = self.chr_seq[variant_pos + s-1: variant_pos + e-1]

        return seq[::-1]

    def get_gene_seq(self, chr, gene_object, variant_pos):
        self.read_chr(chr)
        self.gene_seq = ''
        
        offset = 0
        quit = False
        if gene_object.direction == '+':
            for r in gene_object.range:
                if gene_object.start[0] > r[1]:
                    continue
                if quit:
                    break
                
                if gene_object.start[0] > r[0]:
                    s = gene_object.start[0]
                else:
                    s = r[0]
                
                if gene_object.stop[0] <= r[1]:
                    e = gene_object.stop[0]
                    quit = True
                else:
                    e = r[1] + 1
                
                if variant_pos < e:
                    if variant_pos >= s:
                        offset = offset + variant_pos - s
                else:
                    offset = offset + e - s
                
                #print r, chr, variant_pos, offset, gene_object.start, gene_object.stop, gene_object.direction
                self.gene_seq = self.gene_seq + self.chr_seq[s-1:e-1]
                #print r, gene_object.start, gene_object.stop, variant_pos, offset, len(self.gene_seq), s, e
        else:
            for r in gene_object.range:
                if gene_object.stop[1] > r[1]:
                    continue
                if quit:
                    break
                
                if gene_object.stop[1] > r[0]:
                    s = gene_object.stop[1] + 1
                else:
                    s = r[0]
                
                if gene_object.start[1] <= r[1]:
                    e = gene_object.start[1] + 1
                    quit = True
                else:
                    e = r[1] + 1
                
                #print chr, variant_pos, s, e, self.chr_seq[s-1:e-1]
                if variant_pos >= s:
                    if variant_pos < e:
                        offset = offset + e - variant_pos - 1
                else:
                    offset = offset + e - s
                
                self.gene_seq = self.gene_seq + self.chr_seq[s-1:e-1]
                #print r, gene_object.start, gene_object.stop, variant_pos, offset, len(self.gene_seq), s, e, self.chr_seq[s-1:e-1]
            
            #print self.gene_seq[offset-2 : offset+2]
            self.gene_seq = self.gene_seq[::-1]
            #print self.gene_seq[offset-2 : offset+2]
        
        return [self.gene_seq, offset]
    
    
def complement(s):
    s = list(s)
    for i in xrange(len(s)):
        s[i] = cmpl[s[i]]
    return ''.join(s)

def gene_insert(dict, chr, gene_id, d):
    chr = chr[3:]
    
    if (chr not in dict):
        dict[chr] = {}

    if (gene_id not in dict[chr]):
        dict[chr][gene_id] = gene()
    
    dict[chr][gene_id].direction = d

def gn_read(gname):
    f = open(gname)
    ret = {}
    for line in f:
        p = line[:-1].split('\t')
        ret[p[0]] = [p[1], p[2]]
    f.close()
    return ret

def gtf_read(gtf):
    f = open(gtf)
    
    ret = {}
    
    cur_gene = ''
    for line in f:
        p = line[:-1].split('\t');
        if p[2] == 'CDS':
            continue
        
        gene_id = p[8][9:27]
        if cur_gene != gene_id:
            cur_gene = gene_id
            gene_insert(ret, p[0], gene_id, p[6])
        
        if 'exon' == p[2]:
            ret[p[0][3:]][gene_id].range.append((long(p[3]), long(p[4])))
        
        if 'start_codon' == p[2]:
            ret[p[0][3:]][gene_id].start = (long(p[3]), long(p[4]))

        if 'stop_codon' == p[2]:
            ret[p[0][3:]][gene_id].stop = (long(p[3]), long(p[4]))
    
    
    f.close()
    return ret
    
def find_amino(chr, seq_object, gene_object, variant_pos, origin_nuc, variant_nuc):
    [seq, offset] = seq_object.get_gene_seq(chr, gene_object, variant_pos)
    
    if seq[offset] != origin_nuc:
        print 'error:', chr, variant_pos
        sys.stderr.write('error:' + chr + ', ' + str(variant_pos))
        return 'error_original_seq_not_match'
    
    origin_seq = seq[offset / 3 * 3 : offset / 3 * 3 + 3]
    if gene_object.direction == '-':
        origin_seq = complement(origin_seq)
    origin_amino = amino_codon[origin_seq]
    
    t = list(origin_seq)
    if gene_object.direction == '-':
        variant_nuc = complement(variant_nuc)    
    t[offset % 3] = variant_nuc
    variant_seq = ''.join(t)
    variant_amino = amino_codon[variant_seq]
    
    if origin_amino == variant_amino:
        return variant_interpreter.get(6)
    
    return'p.' + origin_amino + str(offset / 3 + 1) + variant_amino + '|MODERATE'

def append_type(variant_type, t1, trancript_id):
    gid = gene_name.get(trancript_id, ['',''])
    variant_type.append(t1 + '|' + trancript_id + '|' + gid[0] + '|' + gid[-1]);

def detect_variant_type(seq_object, origin_nuc, variant_nuc, gene_data, chr, sorted_key, i, variant_pos):
    length = len(sorted_key)
    variant_type = []
    
    for i in xrange(i, length):
        if gene_data[chr][sorted_key[i]].range[0][0] > variant_pos:
            append_type(variant_type, variant_interpreter.get(2), '')
            break

        if gene_data[chr][sorted_key[i]].range[-1][1] >= variant_pos:
            break
        
    for j in xrange(i, length):
        if gene_data[chr][sorted_key[j]].range[0][0] > variant_pos:
            break
        if gene_data[chr][sorted_key[j]].range[-1][1] < variant_pos:
            continue
        
        is_intron = True
        for exon in gene_data[chr][sorted_key[j]].range:
            if variant_pos >= exon[0] and variant_pos <= exon[1]:
                is_intron = False
                if gene_data[chr][sorted_key[j]].start == None :
                    append_type(variant_type, variant_interpreter.get(5), sorted_key[j])
                    break
                if gene_data[chr][sorted_key[j]].direction == '+':
                    if variant_pos < gene_data[chr][sorted_key[j]].start[0]:
                        append_type(variant_type, variant_interpreter.get(3), sorted_key[j])
                        break
                    elif variant_pos >= gene_data[chr][sorted_key[j]].stop[0]:
                        append_type(variant_type, variant_interpreter.get(4), sorted_key[j])
                        break
                else:
                    if variant_pos >= gene_data[chr][sorted_key[j]].start[1]:
                        append_type(variant_type, variant_interpreter.get(3), sorted_key[j])
                        break
                    elif variant_pos <= gene_data[chr][sorted_key[j]].stop[1]:
                        append_type(variant_type, variant_interpreter.get(4), sorted_key[j])
                        break                   
                append_type(variant_type, find_amino(chr, seq_object, gene_data[chr][sorted_key[j]], variant_pos, origin_nuc, variant_nuc), sorted_key[j])

        if is_intron:
            append_type(variant_type, variant_interpreter.get(1), sorted_key[j])
    
    return [variant_type, i] 
    

def annotate_variant(gtf, fa, vcf, gname, fileout):
    gene_data = gtf_read(gtf)
    global gene_name
    gene_name = gn_read(gname)
    seq_object = genome(fa)
    fv = open(vcf)
    fo = open(fileout, 'w')
    
    output = Popen('grep -v "^#" ' + vcf + ' | wc -l', stdout=PIPE, shell=True)
    Total = int(output.stdout.read().split()[0])
    print 'Total: ', Total

    count = 0
    cur_chr = ''
    info = 0
    for line in fv:
        if line[0] == '#': 
            if info == 1 and line[2:6] != 'INFO':
                fo.write('##INFO=<ID=ANN,Number=1,Type=String,Description="variant annotation">\n')
                info = 2
            elif info == 0 and line[2:6] == 'INFO':
                info = 1
            fo.write(line)
            continue
        
        p = line.split()
        if cur_chr != p[0]:
            cur_chr = p[0]
            
            sorted_list = sorted(gene_data[cur_chr].items(), key=lambda k:k[1].range[0][0])
            sorted_key = [key[0] for key in sorted_list]
            i = 0
        
        variant_pos = long(p[1])
        origin_nuc  = p[3]
        variant_nuc = p[4]
        
        # get variant type
        [variant_type, i] = detect_variant_type(seq_object, origin_nuc, variant_nuc, gene_data, cur_chr, sorted_key, i, variant_pos)
        
        info = 'ANN=' + ','.join(variant_type)
        
        p[7] = p[7] + ';' + info
        fo.write('\t'.join(p) + '\n')
        
        
        count = count + 1
        if count %10 ==0:
            print ('\rCurrent: ' + str(count) + '\t' + str(round(float(count) / Total * 100, 2)) + '%'),
    
    print ('\rCurrent: ' + str(count) + '\t' + str(round(float(count) / Total * 100, 2)) + '%')
    
    fv.close()
    fo.close()

def append_indel_type(indel_type, trancript_id):
    indel_type.append(gene_name.get(trancript_id, ['',''])[-1]);

def detect_indel_type(gene_data, chr, indel_pos):
    indel_type = []
    
    for trans_id, gene in gene_data[chr].iteritems():
        if indel_pos >= gene.range[0][0] and indel_pos <= gene.range[-1][1]:
            indel_type.append(gene_name.get(trans_id, ['',''])[-1])
    
    return indel_type
    
def annotate_indel(gtf, indel, gname, fileout):
    gene_data = gtf_read(gtf)
    global gene_name
    gene_name = gn_read(gname)
    
    fv = open(indel)
    fo = open(fileout, 'w')
    
    cur_chr = ''
    for line in fv:
        if line[0] == '#':
            continue
        p = line[:-1].split()
        
        chrom = p[0].split('_')[0]
        coor = p[0].split('_')[1]
        
        if cur_chr != chrom:
            cur_chr = chrom
            
            sorted_list = sorted(gene_data[cur_chr].items(), key=lambda k:k[1].range[0][0])
            sorted_key = [key[0] for key in sorted_list]
            i = 0
        
        indel_type = detect_indel_type(gene_data, cur_chr, int(coor))
        
        p.append(','.join(indel_type))
        
        p.append(str(len(p[2].split(',')[0]) - len(p[1])))
        
        fo.write('\t'.join(p) + '\n')
    
    fv.close()
    fo.close()

        
if __name__=='__main__':
    if len(sys.argv) != 7 and len(sys.argv) != 6:
        print 'invalid arguments.'
        print 'Usage:'
        print '\tvariant_annotate.py -indel/-vcf gtf_file genome_file vcf_file gene_name_file output_file\n'
        print '\t-vcf         :       annotate vcf file'
        print '\t-indel       :       annotate indel file'
        print '\tgtf_file       :       A GTF file sorted first by chromosome then by gene_id then by coordinate.'
        print '\tgenome_file    :       A genome sequence file.'
        print '\tvcf_file       :       The VCF file to be annotated. Must be sorted first by chromosome then by coordinate.'
        print '\tgene_name_file :       A file with transcript id on the 1st column, gene id on the 2nd column and gene name on the 3rd column.'
        print '\toutput_file    :       A designated output file.'
        exit(1)
    
    if 'indel' in sys.argv[1]:
        annotate_indel(sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5] )
    if 'vcf' in sys.argv[1]:
        annotate_variant(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],sys.argv[6] )
    
