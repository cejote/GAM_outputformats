library(tidyverse)
library(future)

plan(multiprocess)


options(scipen=999999999)
options(stringsAsFactors=F)


readmatrix=function(fn, upper=T, maxdist=MAXDIST){
  print(fn)
  #2021 updated clean version using dplyr
  inmatrix=read_tsv(fn, col_names = T, na = c("", "NA", "NaN"), col_types = "c") %>% #, col_types = cols(...1="c", .default = "d")) %>%
    pivot_longer(cols=-1, names_to="bin2", values_to="value") %>% 
    separate(1, into=c("chrom_x","start_x", "end_x")) %>%
    separate(bin2, into=c("chrom_y","start_y", "end_y")) %>%
    #remove diagonal & return upper part of matrix
    mutate(start_x=as.numeric(start_x), start_y=as.numeric(start_y)) %>%
    dplyr::filter(start_y>start_x, (start_y-start_x)<maxdist)
  return(inmatrix)
}

getlist=function(fn){
  contacts=tibble()
  for(chm in rev(chrlist)){
    print(sprintf(fn, chm))
    #contacts = contacts %>% bind_rows(readmatrix(sprintf(fn, chm)))
    
    f=future(
      write_tsv(readmatrix(sprintf(fn, chm)),
      str_replace(sprintf(fn, chm), ".txt.gz", ".long.gz"), col_names=F)
    )
    
  }
  #return(contacts)
}





chrlist = paste0("chr", c(1:22, "X", "Y"))
res=50000

#maxdist
MAXDIST=Inf
MAXDIST.bin=Inf #MAXDIST/res

# should make use of map function to speed-up reading...
setwd("/data/pombo/christoph/datasets/GAM/coseg.H1hESC.hg38.merged.R1.R2/")
getlist("H1hESC.hg38.merged.R1.R2.segregation_at50000.coseg.%s.txt.gz")
