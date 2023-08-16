import argparse
import pandas as pd

parse = argparse.ArgumentParser(description="prepare for inext")
parse.add_argument("--input", help="input file", required=True)
parse.add_argument("--source", help="the source of data", required=True)
parse.add_argument("--out", help="output file", required=True)
args = parse.parse_args()

df=pd.read_csv(args.input)
RF,genome,source=df['RiPP family'],df['contigs/genomes'],df['source']


data=pd.DataFrame(columns=[args.source])
for i,j,k in zip(RF,genome,source):
    if k in data:
       if j in data[k]:
          if i not in data[k][j]:
             data[k][j][i]=1
       else:
          data[k][j]={}
          data[k][j][i]=1


_data={}
head={}
for ss in data:
    _data[ss]={}
    head[ss]=[]
    for h in data[ss]:
        head[ss].append(h)
        for c in data[ss][h]:
            if c in _data[ss]:
               _data[ss][c][h]=data[ss][h][c]
            else:
               _data[ss][c]={}
               _data[ss][c][h]=data[ss][h][c]


fo2=open(args.out,'w')
wl='\t'
for j,ele in enumerate(list(set(head[args.source]))):
    wl+='%s\t'%j
fo2.write('%s\n'%wl.strip('\t'))

for i in range(c):
    wl=''
    if i in _data[args.source]:
       for ele in list(set(head[args.source])):
           if ele in _data[args.source][i]:
              wl+='1\t'
           else:
              wl+='0\t'
    else:
       for ele in list(set(head[args.source])):
           wl+='0\t'
    fo2.write('%s\n'%wl.strip('\t'))
fo2.close()
