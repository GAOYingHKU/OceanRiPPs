reference: https://gitlab.com/mcfrith/last
1. construct lastal database
   lastdb -P10 -c -m1111110  ripp_db/ripp_db_nt ripp_nt.fasta
2. mapping
   lastal -P 50 -Q 1 -e 120 ripp_db/ripp_db_nt xx.fastq >  map/xx.maf
3. reformat
   maf-convert blasttab map/xx.maf > reformat/xx.txt

