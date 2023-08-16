ref=['k__','p__','c__','o__','f__','g__','s__']

fo=open('ncbi_lineages_2023-04-10_reformat.csv','w')
with open('ncbi_lineages_2023-04-10_simp.csv','r') as fi:
     for line in fi:
         if 'tax_id' in line:
            fo.write(line)
         else:
            hh=line.split(',')[0]
            mem=line.strip('\n').split(',')[1:]
            wl=''
            rec=''
            for i,j in zip(ref,mem):
                if j=='':
                   fin='unassigned%s'%rec
                else:
                   rec=j
                   fin=j
                wl+='%s%s,'%(i,fin)
            fo.write('%s\t%s\n'%(hh,wl.strip(',')))
fo.close()
