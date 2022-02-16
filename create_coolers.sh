#!/bin/bash

for res in 1000000 500000 250000 200000 100000 50000 40000 30000 25000 20000 ; do
  echo $res
  BINFILE=bins.$res.bed
  CHROMSIZES_FILE=../hg38.chrom.sizes
  COOLFILE="./NPMI.H1hESC.hg38.merged.R1.R2.segregation_at${res}.cool"

  #cooler makebins $CHROMSIZES_FILE $res > $BINFILE
  zcat ./${res}/*long.gz* | cooler load -f bg2 --count-as-float --input-copy-status unique --assembly hg38 ${CHROMSIZES_FILE}:${res} - $COOLFILE &
done

