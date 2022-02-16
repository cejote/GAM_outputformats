!/bin/bash

for res in 1000000  500000 250000 200000 100000 50000 40000 30000 25000 20000; do
  mkdir ${res}
  for chr1 in {1..22} X Y; do
    for chr2 in {1..22} X Y; do
      echo $res $chr1 $chr2
      python3 ./generate_matrices_220124.py /data/pombo/Sasha/gamtools/H1_hESC_processed/Segregation_tables/Curated/Merged_PMR20_OW60/H1hESC.hg38.merged.R1.R2.segregation_at$res.table NPMI chr$chr1 chr$chr2 ./${res}/ &
      #nice -n 15 /fast/home/cthieme/.guix-profile/bin/python3 ./generate_matrices_220124.py /data/pombo/Sasha/gamtools/H1_hESC_processed/Segregation_tables/Curated/Merged_PMR20_OW60/H1hESC.hg38.merged.R1.R2.segregation_at$res.table NPMI chr$chr1 chr$chr2 ./${res}/ &
    done
    #wait
  done
  #wait
done
#wait
