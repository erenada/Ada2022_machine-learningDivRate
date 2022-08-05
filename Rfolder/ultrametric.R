library(ape)



#location
base_dir <- "/Users/eren/Documents/GitHub/Chapter3/"

input_dir <- "/Users/eren/Documents/GitHub/Chapter3/AlexTrees"

output_dir <- "/Users/eren/Documents/GitHub/Chapter3/UltrametricTrees"

tree_data_sets <- c(list.files(output_dir)) # char or list??


for(treefile in list.files(paste(input_dir,"Dataset1",sep = "/"))){
    treeobject <- read.tree(paste(input_dir,"Dataset1",treefile,sep = "/"))
    treeobject <- root(treeobject, outgroup="5", resolve.root=T)
    chronotree <- chronoMPL(treeobject)
    print(is.ultrametric(chronotree))
    write.tree(chronotree, file = paste(output_dir,"Dataset1",treefile,sep = "/"))
}


test_tree <- read.tree("../UltrametricTrees/Dataset1/inference_fong_best600.txt.treefile")

test_tree<- root(test_tree, outgroup = "5", resolve.root = T)

test_tree <- chronoMPL(test_tree)

plot(test_tree)




