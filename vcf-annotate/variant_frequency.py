import sys
from subprocess import Popen, PIPE

def read_gene_mutation_frequency(in_file):
    f = open(in_file)
    gene = {}
    line = f.readline()
    
    for line in f:
        p = line.split()
        p.append(float(p[1]) / float(p[2]))
        p.append([])
        p.append(0)
        
        gene[p[0]] = p

    f.close()
    
    return gene
    
    
def count_gene_mutation(vcf, gene):
    f = open(vcf)
    
    for line in f:
        if line[0] == '#':
            continue
        
        p = line[:-1].split('\t')
        info = p[7].split(';')
        for field in info:
            if field[:2] == 'TL':
                q = field[3:]
            if field[:3] == 'ANN':
                if 'intergenic_region' not in field:
                    ann = field.split(',')
                    nl =[]
                    for a in ann:
                        name= a.split('|')[-1]
                        nl.append(name)
                    ns = set(nl)
                    for name in ns:
                        if name in gene:
                            gene[name][4].append(q)
                            gene[name][5] = gene[name][5] + 1
    
    f.close()


def write_to_file(gene,f):
    f = open(f,'w')
    
    for line in gene:
        line[3] = str(line[3])
        line[4] = ','.join(line[4])
        line[5] = str(line[5])
        
        f.write('\t'.join(line) + '\n')
    
    f.close()
    
if __name__=='__main__':
    if len(sys.argv) != 4:
        exit(1)
        
    gene = read_gene_mutation_frequency(sys.argv[1])
    
    count_gene_mutation(sys.argv[2], gene)
    
    gene = sorted(gene.itervalues(), key=lambda k:k[5], reverse=True)
    write_to_file(gene,sys.argv[3])