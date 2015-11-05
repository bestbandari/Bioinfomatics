# # # # # # # # # # # # # # # # # # # # # #
# Author: Jun He
#
# Usage:
#    input_file         :       The input file to be sorted.
#    chromosome_column  :       The column of chromosome in the input file.
#    gene_column        :       The column of gene_id in the input file.
#    coordinate_column  :       The column of coordinate in the input file.
#    gene_name_file     :       A mutect output containing gene name. (Optional)
#    output_file        :       A designated output file.
#
# Note: When the number of input arguments is not correct, an usage message will be printed on the screen.
#
# # # # # # # # # # # # # # # # # # # # # #
import sys


def mysort(filein, chr, gene, pos, fileout):
    fi = open(filein)
    
    data = []
    for line in fi:
        p = line[:-1].split('\t')
        data.append(p)

    fi.close()
    
    data = sorted(data, key=lambda i: (100 if 'X' in i[chr] else 
                                       101 if 'Y' in i[chr] else
                                       102 if 'M' in i[chr] else
                                       int(i[chr][3:]) if i[chr][0] == 'c' else 
                                       int(i[chr]), i[gene], int(i[pos])))

    fo = open(fileout, 'w')
    for row in data:
        line = '\t'.join(row) + '\n'
        fo.write(line)

    fo.close()


def usage():
    print 'usage:'
    print '\tsortby_chr_pos.py input_file chromosome_column gene_column position_column output_file'
    print 'e.g.'
    print '\t python sortby_chr_pos.py dogGenes.gtf 0 8 3 dogGenes.sorted.gtf'

if __name__=='__main__':
    if len(sys.argv) != 6:
        usage()
        exit(1)

    mysort(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
