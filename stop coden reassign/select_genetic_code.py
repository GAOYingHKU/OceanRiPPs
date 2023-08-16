from Bio import SeqIO
import numpy as np
import argparse


parse = argparse.ArgumentParser(description="select genetic code")
parse.add_argument("--input_nostand", help="input file", required=True)
parse.add_argument("--input_stand", help="input file", required=True)
parse.add_argument("--out", help="output file", required=True)
args = parse.parse_args()





elv={}
with open(args.input_stand,'r') as fi:
     for line in fi:
         elv[line.split()[0]]=line.strip('\n').split()[-1]

ni={}
with open(args.input_nostand,'r') as fi:
     for line in fi:
         ni[line.split()[0]]=line.strip('\n').split()[-1]


fo=open(args.out,'w')
for ele in elv:
    delta=float(ni[ele])-float(elv[ele])
    if float(ni[ele])<150 or float(elv[ele])<150:continue
    try:
       if (delta/float(ni[ele])+delta/float(elv[ele]))/2>=0.2:
          fo.write('%s\t%s\t%s\n'%(ele,elv[ele],ni[ele]))
    except ZeroDivisionError:
       print(ele,elv[ele],ni[ele])

fo.close()
