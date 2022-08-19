

library(ape)
library(geiger)
library(phytools)


#s_tree1 <- read.tree("./AlexOriginalTrees/Dataset1/s_tree.trees") 

inf_tree1 <- read.tree('./AlexOriginalTrees/Dataset1/inferenceTest.treefile')

inf_best800 <- read.tree('./AlexOriginalTrees/Dataset1/inference_fong_best800.txt.treefile')

inf_best600 <- read.tree('./AlexOriginalTrees/Dataset1/inference_fong_best600.txt.treefile')

inf_worst600 <- read.tree('./AlexOriginalTrees/Dataset1/inference_fong_worst600.txt.treefile')

s_treedataset2 <- read.tree("./AlexTrees/Dataset2/s_tree.trees")

i_treedataset2 <- read.tree("./AlexTrees/Dataset2/inferenceTest.treefile")


plot(inf_worst600)
plot(i_treedataset2)


inf_best600 <- root(inf_best600, "5", resolve.root = T)

plot(inf_best600)


chr_inf_inf_best600 <- chronoMPL(inf_best600)


plot(chr_inf_inf_best600)


gammaStat(chr_inf_inf_best600)


true_tree <- read.tree("./AlexOriginalTrees/Dataset1/s_tree.trees")

plot(true_tree)

plot(inf_tree1)


is.rooted(true_tree)

gammaStat(true_tree)
gammaStat(inf_tree1)
gammaStat(reltimetree)



is.ultrametric(inf_tree1)

plot(true_tree)


inf_tree1 <- root(inf_tree1, "5",resolve.root = T)


plot(inf_tree1)


#

inf_tree1 <- chronoMPL(inf_tree1)


is.ultrametric(chr_ultrametric_inf_best800)


plot(chr_ultrametric_inf_best800)


gammaStat(chr_ultrametric_inf_best800)

ltt(chr_ultrametric_inf_best800,log.lineages = T,gamma = T) + title(main = "chr_ultrametric_inftree1")



###
unrs_inftree <- read.tree("../Chapter3/AlexOriginalTrees/Dataset1/inferenceTest.treefile")

is.rooted(unrs_inftree)

plot(unrs_inftree)

unrs_inftree <- root(reltimetree,"5",resolve.root = F)

plot(unrs_inftree)

gammaStat(unrs_inftree)

plot.phylo(unrs_inftree) + title(main = "unrs_inftree")

is.ultrametric(unrs_inftree)

unrs_inftree <- chronos(unrs_inftree, lambda = 1, model = "discrete", quiet = FALSE)

plot(unrs_inftree)

gammaStat(unrs_inftree)

ltt(unrs_inftree,gamma = T,log.lineages = T) + title(main="unrs_inftree-ltt")


## Reltime

#figtree <- read.tree("./RelTime/inf_fig_alt.nwk")


#is.rooted(figtree)

#plot(figtree)


##

inf_tree1$node.label <- NULL

write.tree(inf_tree1, file="./RelTime/alexRT.nwk")



reltimetree <- read.tree("./RelTime/Newick Export.nwk")

plot(reltimetree)

is.ultrametric(reltimetree)

ltt(chr_ultrametric_inf_best800,gamma = T,log.lineages = T) + title(main="reltimetree-ltt")

ltt(true_tree,gamma = T,log.lineages = T) + title(main="reltimetree-ltt")


min(reltimetree$edge.length)


scaled <- geiger::rescale(chr_ultrametric_inf_best800, model = "depth", 1)


plot(scaled)


gammaStat(scaled)

##

unroot

