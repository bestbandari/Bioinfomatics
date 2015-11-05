import sys

def read_gene_mutation(fi):
    f = open(fi)
    gene = {}
    
    for line in f:
        p = line[:-1].split('\t')
        p[5] = int(p[5])
        gene[p[0]] = p

    f.close()
    
    return gene

def write_to_file(gene,f):
    f = open(f,'w')
    
    for line in gene:
        line[5] = str(line[5])
        
        f.write('\t'.join(line) + '\n')
    
    f.close()
    

def add_indels(f_indel, f_gene, f_out):
    gene = read_gene_mutation(f_gene)
    f = open(f_indel)
    
    for line in f:
        p = line[:-1].split('\t')
        name = p[5].split(',')[0]
        
        if name in gene:
            gene[name][5] = gene[name][5] + 1
    
    top = sorted(gene.itervalues(), key=lambda k:k[5], reverse=True)
    write_to_file(top, f_out)
    
    
    
if __name__=='__main__':
    if len(sys.argv) != 4:
        print 'Usage:'
        print 'python add_indel_frequency.py indels_in_cancer_not_in_blood.annotated.indel top_dog.tsv top_dog.add_indel.tsv'
        exit(1)
        
    add_indels(sys.argv[1], sys.argv[2], sys.argv[3])
        