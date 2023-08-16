1. search against NT database using Megablast
   blastn -task megablast -query xx -db /data/database/blastdb/nt -out megablast_out.tbl -outfmt "6 std qlen slen staxid ssciname sskingdom" 
2.create reformated ncbi lineages
  2a. pip install -U ncbitax2lin; get the "ncbi_lineages_2023-04-10.csv.gz" according to the guidance (https://github.com/zyxue/ncbitax2lin)
  2b. cut -d ',' -f1-8 ncbi_lineages_2023-04-10.csv > ncbi_lineages_2023-04-10_simp.csv
  2c. python reformat_ncbi_lineage.py (output: ncbi_lineages_2023-04-10_reformat.csv)
3.voting-based taxonomy classification
  
