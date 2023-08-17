data={}
with open('host_phage_pair_info.txt','r') as fi:
     for line in fi:
         a=line.split('\t')[0]
         b=line.split('\t')[1]
         c=line.strip('\n').split('\t')[2]
         d='%s|%s'%(a,b)
         if d in data:
            if c in data[d]:
               data[d][c]+=1
            else:
               data[d][c]=1
         else:
            data[d]={}
            data[d][c]=1

head=['CRISPR','tRNA','prophage']
fo=open('host_phage_pair_summary.00.txt','w')
wl='host\tphage\t'
for h in head:
    wl+='%s\t'%h
fo.write('%s\n'%wl.strip('\t'))
for i in data:
    wl='%s\t'%i.replace('|','\t')
    for j in head:
        if j in data[i]:
           wl+='%s\t'%str(data[i][j])
        else:
           wl+='0\t'
    fo.write('%s\n'%wl.strip('\t'))

fo.close()
