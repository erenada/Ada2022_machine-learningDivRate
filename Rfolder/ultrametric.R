library(ape)
library(geiger)



#location
base_dir <- "/Users/eren/Documents/GitHub/Chapter3/"

input_dir <- "/Users/eren/Documents/GitHub/Chapter3/AlexTrees"

output_dir <- "/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees"

tree_data_sets <- c(list.files(output_dir)) # char or list??

treefile <- "inference_fong_best600.txt.treefile"

dataset1_outgroup <- 5
dataset2_outgroup 
dataset3_outgroup
dataset4_outgroup 

  
  
## outgroup loop  BUT NEEDS TO CHANGE!! 

#for(dataset in list.dirs(paste(input_dir),recursive = F,full.names = F)){
#  for(treefile in list.files(paste(input_dir,dataset,sep = "/"))){
#    treeobject <- read.tree(paste(input_dir,dataset,treefile,sep = "/"))
#    treeobject <- root(treeobject, outgroup="5", resolve.root=T)
#    chronotree <- chronoMPL(treeobject)
#    print(is.ultrametric(chronotree))
#    write.tree(chronotree, file = paste(output_dir,"Dataset1",treefile,sep = "/"))
#}
#}





test_tree <- read.tree("../AlexOriginalTrees/Dataset1/s_tree.trees")

is.ultrametric(test_tree)

test_tree1 <- geiger::rescale(test_tree, model = "depth", 1)


min(test_tree$edge.length)


test_tree<- root(test_tree, outgroup = "5", resolve.root = T)

test_tree <- chronoMPL(test_tree)

plot(test_tree)






