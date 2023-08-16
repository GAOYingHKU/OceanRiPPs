1. run prodigal using the standard code (11), opal code (25), amber code (15) and ochre code (90)
   prodigal source code with an added ochre code (90) are available on JGI public ftp site https://portal.nersc.gov/dna/microbial/prokpubs/recoding/
   i.e. prodigal -i in.fasta -a out_90.fas -d out_90.fna -o out_90.gbk -p meta -g 90
2. extract coding score of different genetic codes
   python extract_coding_score.py --input out_90.gbk --out out_90.txt
3. 
