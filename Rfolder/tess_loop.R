## ultimate calculation of div rates.

library(ape)
library(phangorn)
library(dplyr)
library(TESS)

## Birth-death processes with constant rates


input_dir <- "/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees"

test_input_dir <- "/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees/Dataset1/"

out_dir <- "/Users/eren/Documents/GitHub/Chapter3/outdir"

model <- "Dataset1"

for(tree in list.files(test_input_dir)){
      print(paste(tree))
  tree_object <- read.tree(paste(test_input_dir,tree, sep = "")) ## will change based on final directory structure
  
  #variables:
  #tree time
  
  times <- ass.numeric(branching.times(tree_object))
  
  }