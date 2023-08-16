import argparse

parse = argparse.ArgumentParser(description="extract coding score")
parse.add_argument("--input", help="input file", required=True)
parse.add_argument("--out", help="output file", required=True)
args = parse.parse_args()


score={}
with open(args.input,'r') as fi:
     for line in fi:
         if line.split()[0]=='DEFINITION':
            idx=line.split('seqhdr="')[1].split('";')[0]
            score[idx]=0
         else:
            if ';score=' in line:
               score[idx]+=float(line.split(';score=')[1].split(';')[0])

fo=open(args.out,'w')
for idx in score:
    fo.write('%s\t%s\n'%(idx,score[idx]))
fo.close()
