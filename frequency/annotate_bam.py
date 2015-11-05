import sys,re

def annotate(pile,output):
    result = []
    
    f = open(pile)
    r = re.compile('[a-zA-Z]')
    for line in f:
        p = line[:-1].split('\t')
        p[4] = ''.join(r.findall(p[4]))
        result.append(p)
    f.close()
    
    result = sorted(result, key=lambda k:len(k[4]), reverse=True)
    
    f = open(output, 'w')
    for item in result:
        f.write('\t'.join(item) + '\n')
    
    f.close()
        


if __name__=='__main__':
    if len(sys.argv) != 3:
        print 'Invalid arguments'
        exit(1)
    
    annotate(sys.argv[1], sys.argv[2])
    print 'Finished successfully'