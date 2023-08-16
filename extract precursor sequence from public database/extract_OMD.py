import os
from Bio import SeqIO
import argparse


parse = argparse.ArgumentParser(description="extract precursor sequences of BGCs from OMD database")
parse.add_argument("--infolder", help="input folder", required=True)
parse.add_argument("--out", help="output file", required=True)
args = parse.parse_args()



ripp=['bottromycin','cyanobactin','cyclic-lactone-autoinducer','epipeptide','fungal-RiPP','glycocin','LAP','lantipeptide','lassopeptide','linaridin','lipolanthine','microviridin','proteusin','ranthipeptide','RaS-RiPP','redox-cofactor','RiPP-like','sactipeptide','spliceotide','thioamitides','thiopeptide','bacteriocin','head_to_tail','lanthidin','TfuA-related','microcin']
ripp_exclude=['DUF692','PF13471','Asn_synthase','PF00881','YcaO','PF04055','TIGR03975','TIGR01193','cyanobactin_synth','TfuA','PoyD','cypI','Lant_dehydr_C','Lant_dehydr_N','skfc','TIGR03882','TIGR01716','PF06968','PF07366','Flavoprotein']

sig=False
fo=open(args.out,'w')
for fn in os.listdir(args.infolder):
    fa=SeqIO.parse(indir+fn,'genbank')
    try:
       for gb_seq in fa:
           for ele in gb_seq.features:
               if ele.type=='cand_cluster':
                  product=ele.qualifiers['product'][0]
                  if product in ripp:
                     sig=True
                  else:
                     sig=False
               if ele.type=='CDS':
                  if sig:
                     try:
                        gf=ele.qualifiers['gene_functions'][0]
                        if gf=='':continue
                        if product in gf:
                           if ': ' in gf and gf.split(': ')[-1] in ripp_exclude:continue
                           seq=ele.qualifiers['translation'][0]
                           if seq=='':continue
                           fo.write('>%s___%s\n%s\n'%(fn,gf,seq))
                     except KeyError:
                        continue
    except UnicodeDecodeError:
       print(fn)
fo.close()
