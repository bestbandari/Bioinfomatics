# # # # # # # # # # # # # # # # # # # # # #
# Author: Jun He
#
# Usage:
#    input_file         :       The input file to be sorted.
#    chromosome_column  :       The column of chromosome in the input file.
#    coordinate_column  :       The column of coordinate in the input file.
#    gene_name_file     :       A mutect output containing gene name. (Optional)
#    output_file        :       A designated output file.
#
# Note: When the number of input arguments is not correct, an usage message will be printed on the screen.
#
# # # # # # # # # # # # # # # # # # # # # #
import sys


def mysort(filein, chr, pos, fileout):
    fi = open(filein)
    fo = open(fileout, 'w')
    
    data = []
    for line in fi:
        if line[0] == '#':
            fo.write(line)
        else:
            p = line[:-1].split('\t')
            data.append(p)

    fi.close()
    
    data = sorted(data, key=lambda i: (100 if 'X' in i[chr] else 
                                       101 if 'Y' in i[chr] else
                                       102 if 'M' in i[chr] else
                                       int(i[chr][3:]) if i[chr][0] == 'c' else 
                                       int(i[chr]), int(i[pos])))

    
    for row in data:
        line = '\t'.join(row) + '\n'
        fo.write(line)

    fo.close()


def usage():
    print 'usage:'
    print '\tsortby_chr_pos.py input_file chromosome_column coordinate_column output_file'
    print 'e.g.'
    print '\t python sortby_chr_pos.py mutect.vcf 0 1 mutect.sorted.gtf'    


if __name__=='__main__':
    if len(sys.argv) != 5:
        usage()
        exit(1)

    mysort(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
