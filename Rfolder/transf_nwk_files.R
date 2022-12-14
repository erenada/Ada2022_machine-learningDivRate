## converting tree files to Newick for Reltime use

library(ape)

library(phytools)

library(phangorn)

#Dataset 1

s_tree_Dataset1 <- root(read.tree("./AlexOriginalTrees/Dataset1/s_tree.trees"), "5", resolve.root = T)
inf_test_Dataset1 <- root(read.tree("./AlexOriginalTrees/Dataset1/inferenceTest.treefile"), "5", resolve.root = T)
inf_best800_Dataset1 <- root(read.tree("./AlexOriginalTrees/Dataset1/inference_fong_best800.txt.treefile"), "5", resolve.root = T)
inf_best600_Dataset1 <- root(read.tree("./AlexOriginalTrees/Dataset1/inference_fong_best600.txt.treefile"), "5", resolve.root = T)
inf_worst600_Dataset1 <- root(read.tree("./AlexOriginalTrees/Dataset1/inference_fong_worst600.txt.treefile"), "5", resolve.root = T)
rand0_Dataset1 <- root(read.tree("./AlexOriginalTrees/Dataset1/inference_fong_random0.txt.treefile"), "5", resolve.root = T)
rand1_Dataset1 <- root(read.tree("./AlexOriginalTrees/Dataset1/inference_fong_random1.txt.treefile"), "5", resolve.root = T)
rand2_Dataset1 <- root(read.tree("./AlexOriginalTrees/Dataset1/inference_fong_random2.txt.treefile"), "5", resolve.root = T)
rand3_Dataset1 <- root(read.tree("./AlexOriginalTrees/Dataset1/inference_fong_random3.txt.treefile"), "5", resolve.root = T)


#drop branch support
inf_test_Dataset1$node.label <- NULL
inf_best800_Dataset1$node.label <- NULL
inf_best600_Dataset1$node.label <- NULL
inf_worst600_Dataset1$node.label <- NULL
rand0_Dataset1$node.label <- NULL
rand1_Dataset1$node.label <- NULL
rand2_Dataset1$node.label <- NULL
rand3_Dataset1$node.label <- NULL

##write trees

write.tree(s_tree_Dataset1, file = "./RelTime/input/Dataset1/s_tree.nwk")
write.tree(inf_test_Dataset1, file = "./RelTime/input/Dataset1/inferenceTest.nwk")
write.tree(inf_best800_Dataset1, file = "./RelTime/input/Dataset1/inference_fong_best800.nwk")
write.tree(inf_best600_Dataset1, file = "./RelTime/input/Dataset1/inference_fong_best600.nwk")
write.tree(inf_worst600_Dataset1, file = "./RelTime/input/Dataset1/inference_fong_worst600.nwk")
write.tree(rand0_Dataset1, file = "./RelTime/input/Dataset1/inference_fong_random0.nwk")
write.tree(rand1_Dataset1, file = "./RelTime/input/Dataset1/inference_fong_random1.nwk")
write.tree(rand2_Dataset1, file = "./RelTime/input/Dataset1/inference_fong_random2.nwk")
write.tree(rand3_Dataset1, file = "./RelTime/input/Dataset1/inference_fong_random3.nwk")


#Dataset 2

s_tree_Dataset2 <- read.tree("./AlexOriginalTrees/Dataset2/s_tree.trees")
inf_tree_Dataset2 <-  root(read.tree("./AlexOriginalTrees/Dataset2/inferenceTest.treefile"), "59", resolve.root = T)
inf_best800_Dataset2 <- root(read.tree("./AlexOriginalTrees/Dataset2/inference_wickett_best800.txt.treefile"), "59", resolve.root = T)
inf_best600_Dataset2 <- root(read.tree("./AlexOriginalTrees/Dataset2/inference_wickett_best600.txt.treefile"), "59", resolve.root = T)
inf_worst600_Dataset2 <- root(read.tree("./AlexOriginalTrees/Dataset2/inference_wickett_worst600.txt.treefile"), "59", resolve.root = T)
rand0_Dataset2 <- root(read.tree("./AlexOriginalTrees/Dataset2/inference_wickett_random0.txt.treefile"), "59", resolve.root = T)
rand1_Dataset2 <- root(read.tree("./AlexOriginalTrees/Dataset2/inference_wickett_random1.txt.treefile"), "59", resolve.root = T)
rand2_Dataset2 <- root(read.tree("./AlexOriginalTrees/Dataset2/inference_wickett_random2.txt.treefile"), "59", resolve.root = T)
rand3_Dataset2 <- root(read.tree("./AlexOriginalTrees/Dataset2/inference_wickett_random3.txt.treefile"), "59", resolve.root = T)


#drop branch support
inf_tree_Dataset2$node.label <- NULL
inf_best800_Dataset2$node.label <- NULL
inf_best600_Dataset2$node.label <- NULL
inf_worst600_Dataset2$node.label <- NULL
rand0_Dataset2$node.label <- NULL
rand1_Dataset2$node.label <- NULL
rand2_Dataset2$node.label <- NULL
rand3_Dataset2$node.label <- NULL


# write trees

write.tree(s_tree_Dataset2, file = "./RelTime/input/Dataset2/s_tree.nwk")
write.tree(inf_tree_Dataset2, file = "./RelTime/input/Dataset2/inferenceTest.nwk")
write.tree(inf_best800_Dataset2, file = "./RelTime/input/Dataset2/inference_wickett_best800.nwk")
write.tree(inf_best600_Dataset2, file = "./RelTime/input/Dataset2/inference_wickett_best600.nwk")
write.tree(inf_worst600_Dataset2, file = "./RelTime/input/Dataset2/inference_wickett_worst600.nwk")
write.tree(rand0_Dataset2, file = "./RelTime/input/Dataset2/inference_wickett_random0.nwk")
write.tree(rand1_Dataset2, file = "./RelTime/input/Dataset2/inference_wickett_random1.nwk")
write.tree(rand2_Dataset2, file = "./RelTime/input/Dataset2/inference_wickett_random2.nwk")
write.tree(rand3_Dataset2, file = "./RelTime/input/Dataset2/inference_wickett_random3.nwk")


#Dataset 3
s_tree_Dataset3 <- root(read.tree("./AlexOriginalTrees/Dataset3/s_tree.trees"), "7", resolve.root = T)
inf_tree_Dataset3 <- root(read.tree("./AlexOriginalTrees/Dataset3/inferenceTest.treefile"), "7", resolve.root = T)
inf_best800_Dataset3 <- root(read.tree("./AlexOriginalTrees/Dataset3/inference_liu_best800.txt.treefile"), "7", resolve.root = T)
inf_best600_Dataset3 <- root(read.tree("./AlexOriginalTrees/Dataset3/inference_liu_best600.txt.treefile"), "7", resolve.root = T)
inf_worst600_Dataset3 <- root(read.tree("./AlexOriginalTrees/Dataset3/inference_liu_worst600.txt.treefile"), "7", resolve.root = T)
rand0_Dataset3 <- root(read.tree("./AlexOriginalTrees/Dataset3/inference_liu_random0.txt.treefile"), "7", resolve.root = T)
rand1_Dataset3 <- root(read.tree("./AlexOriginalTrees/Dataset3/inference_liu_random1.txt.treefile"), "7", resolve.root = T)
rand2_Dataset3 <- root(read.tree("./AlexOriginalTrees/Dataset3/inference_liu_random2.txt.treefile"), "7", resolve.root = T)
rand3_Dataset3 <- root(read.tree("./AlexOriginalTrees/Dataset3/inference_liu_random3.txt.treefile"), "7", resolve.root = T)


#drop branch support
inf_tree_Dataset3$node.label <- NULL
inf_best800_Dataset3$node.label <- NULL
inf_best600_Dataset3$node.label <- NULL
inf_worst600_Dataset3$node.label <- NULL
rand0_Dataset3$node.label <- NULL
rand1_Dataset3$node.label <- NULL
rand2_Dataset3$node.label <- NULL
rand3_Dataset3$node.label <- NULL

# write trees

write.tree(s_tree_Dataset3, file = "./RelTime/input/Dataset3/s_tree.nwk")
write.tree(inf_tree_Dataset3, file = "./RelTime/input/Dataset3/inferenceTest.treefile.nwk")
write.tree(inf_best800_Dataset3, file = "./RelTime/input/Dataset3/inference_liu_best800.nwk")
write.tree(inf_best600_Dataset3, file = "./RelTime/input/Dataset3/inference_liu_best600.nwk")
write.tree(inf_worst600_Dataset3, file = "./RelTime/input/Dataset3/inference_liu_worst600.nwk")
write.tree(rand0_Dataset3, file = "./RelTime/input/Dataset3/inference_liu_random0.nwk")
write.tree(rand1_Dataset3, file = "./RelTime/input/Dataset3/inference_liu_random1.nwk")
write.tree(rand2_Dataset3, file = "./RelTime/input/Dataset3/inference_liu_random2.nwk")
write.tree(rand3_Dataset3, file = "./RelTime/input/Dataset3/inference_liu_random3.nwk")


#Dataset 4

s_tree_Dataset4 <- root(read.tree("./AlexOriginalTrees/Dataset4/s_tree.trees"), "85", resolve.root = T)
inf_tree_Dataset4 <- root(read.tree("./AlexOriginalTrees/Dataset4/inferenceTest.treefile"), "85", resolve.root = T)
inf_best800_Dataset4 <- root(read.tree("./AlexOriginalTrees/Dataset4/inference_mcgowen_best800.txt.treefile"), "85", resolve.root = T)
inf_best600_Dataset4 <- root(read.tree("./AlexOriginalTrees/Dataset4/inference_mcgowen_best600.txt.treefile"), "85", resolve.root = T)
inf_worst600_Dataset4 <- root(read.tree("./AlexOriginalTrees/Dataset4/inference_mcgowen_worst600.txt.treefile"), "85", resolve.root = T)
rand0_Dataset4 <- root(read.tree("./AlexOriginalTrees/Dataset4/inference_mcgowen_random0.txt.treefile"), "85", resolve.root = T)
rand1_Dataset4 <- root(read.tree("./AlexOriginalTrees/Dataset4/inference_mcgowen_random1.txt.treefile"), "85", resolve.root = T)
rand2_Dataset4 <- root(read.tree("./AlexOriginalTrees/Dataset4/inference_mcgowen_random2.txt.treefile"), "85", resolve.root = T)
rand3_Dataset4 <- root(read.tree("./AlexOriginalTrees/Dataset4/inference_mcgowen_random3.txt.treefile"), "85", resolve.root = T)

#drop branch support
inf_tree_Dataset4$node.label <- NULL
inf_best800_Dataset4$node.label <- NULL
inf_best600_Dataset4$node.label <- NULL
inf_worst600_Dataset4$node.label <- NULL
rand0_Dataset4$node.label <- NULL
rand1_Dataset4$node.label <- NULL
rand2_Dataset4$node.label <- NULL
rand3_Dataset4$node.label <- NULL

# write trees
write.tree(s_tree_Dataset4, file = "./RelTime/input/Dataset4/s_tree.nwk")
write.tree(inf_tree_Dataset4, file = "./RelTime/input/Dataset4/inferenceTest.nwk")
write.tree(inf_best800_Dataset4, file = "./RelTime/input/Dataset4/inference_mcgowen_best800.nwk")
write.tree(inf_best600_Dataset4, file = "./RelTime/input/Dataset4/inference_mcgowen_best600.nwk")
write.tree(inf_worst600_Dataset4, file = "./RelTime/input/Dataset4/inference_mcgowen_worst600.nwk")
write.tree(rand0_Dataset4, file = "./RelTime/input/Dataset4/inference_mcgowen_random0.nwk")
write.tree(rand1_Dataset4, file = "./RelTime/input/Dataset4/inference_mcgowen_random1.nwk")
write.tree(rand2_Dataset4, file = "./RelTime/input/Dataset4/inference_mcgowen_random2.nwk")
write.tree(rand3_Dataset4, file = "./RelTime/input/Dataset4/inference_mcgowen_random3.nwk")


## check RelTimeTrees & make it ultrametric


## Dataset1-4


for(treefile in list.files("./RelTime/output/Dataset4/")){
  tree <- read.tree(paste("./RelTime/output/Dataset4/",treefile,sep=""))
  tree <- nnls.tree(cophenetic(tree), tree, rooted = TRUE)
  write.tree(tree, file=paste("./RelTime/output/Dataset4/",treefile,sep = ""))
}




