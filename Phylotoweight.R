library(phytools)
library(castor)
library(ips)
library(dplyr)
# filename <- "HW0EZAgBrbloeCRPDLSDmA_newick.txt"
#filename <-"rooted_tree.nwk"
# plot(tree,show.node=TRUE)
phylotoweight<-function(nwkfile){
  tree<-phytools::read.newick(nwkfile)
  root <- find_root(tree)
  numberofnode<-length(tree$node.label)
  numberoftip<-length(tree$tip.label)
  namelist<-c(tree$tip.label,tree$node.label)
  #tree1 <- as(tree,"phylo4")
  # .tipToRoot(tree1,300, root)
  # distances = get_pairwise_distances(tree, c(1),c(526))
  # summary(tree1)
  #
  # .tipToRoot(tree1,527, root)
  start=root
  startweight=100^30
  alllist0<-c()
  allnodelist0<-c()
  allnamelist0<-c()
  alllist0<-append(alllist0,1)
  allnodelist0<-append(allnodelist0,root)
  allnamelist0<-append(allnamelist0,namelist[root])
  for (j in 2:(numberofnode+numberoftip)){
    if (start<root){
      alllist0<-alllist0
      allnodelist0<-allnodelist0
      allnamelist0<-allnamelist0
      start=allnodelist0[j]
      startweight=alllist0[j]
    } else{
      eachdistancelist0<-c()
      nodeidlist0<-c()
      for (i in (1:length(descendants(tree, start, type = "daughter", ignore.tip = TRUE, labels = FALSE)))){
        dist<-get_pairwise_distances(tree, descendants(tree, start, type = "daughter", ignore.tip = TRUE, labels = FALSE)[i], start)
        eachdistancelist0<-append(eachdistancelist0,dist)
        nodeidlist0<-append(nodeidlist0,descendants(tree, start, type = "daughter", ignore.tip = TRUE, labels = FALSE)[i])
      }
      if (sum((eachdistancelist0))==0){
        eachdistancelist0[eachdistancelist0 == 0] <- 1
        eachdistancelist1<-startweight*(eachdistancelist0)/sum((eachdistancelist0))
        alllist0<-append(alllist0,eachdistancelist1)
        allnodelist0<-append(allnodelist0,nodeidlist0)
        allnamelist0<-append(allnamelist0,namelist[nodeidlist0])
        start=allnodelist0[j]
        startweight=alllist0[j]
      } else{
        eachdistancelist1<-startweight*(eachdistancelist0)/sum((eachdistancelist0))
        alllist0<-append(alllist0,eachdistancelist1)
        allnodelist0<-append(allnodelist0,nodeidlist0)
        allnamelist0<-append(allnamelist0,namelist[nodeidlist0])
        start=allnodelist0[j]
        startweight=alllist0[j]}
    }
  }
  phyloweight<-as.data.frame(list(allnamelist0=((allnamelist0)),alllist0=((alllist0)),allnodelist0=((allnodelist0))))
  phyloweight<-filter(phyloweight,allnodelist0 < root)
  phyloweight<-phyloweight[,1:2]
  return(phyloweight)
}
descendants

# get_pairwise_distances(tree, 39162,84691)
#
#
# tree1 <- as(tree,"phylo4")
# library(phytools)
# library(castor)
# library(ips)
# .tipToRoot(tree1, 39162, root)
# .tipToRoot
# library(treeman)
# library(tidytree)
#
# library(adephylo)
#
# library(tidytree)
#
#
#


#
# get_pairwise_distances
# descendants(tree, start, type = "daughter", ignore.tip = TRUE, labels = FALSE)
# allnodelist0
# allnamelist0
# alllist0
# numberofnode+numberoftip
#
#
# eachdistancelist0<-c()
# nodeidlist0<-c()
# for (i in (1:length(descendants(tree, start, type = "daughter", ignore.tip = TRUE, labels = FALSE)))){
#   dist<-get_pairwise_distances(tree, descendants(tree, start, type = "daughter", ignore.tip = TRUE, labels = FALSE)[i], start)
#   eachdistancelist0<-append(dist,eachdistancelist0)
#   nodeidlist0<-append(descendants(tree, start, type = "daughter", ignore.tip = TRUE, labels = FALSE)[i],nodeidlist0)
# }
# eachdistancelist1<-startweight*(eachdistancelist0)/sum((eachdistancelist0))
#
#
# root
# descendants(tree, root, type = "daughter", ignore.tip = TRUE, labels = FALSE)
# get_pairwise_distances(tree, 48118, root)
# get_pairwise_distances(tree, 48132, root)
#
# list0<-list()
# namelist0<-list()
# for (i in (1:(find_root(tree)-1))){
#   totallength<-get_pairwise_distances(tree, i, root)
#   # for (j in (1:(length(.tipToRoot(tree1,i, root))))){
#   #   f<-c(i,.tipToRoot(tree1,i, root))
#   #   tiptonode<-get_pairwise_distances(tree, f[j],f[j+1])
#   #   numberofde<-length(descendants(tree, f[j+1], type = "daughter", ignore.tip = TRUE, labels = FALSE))
#   #   totallength<-append(totallength,tiptonode/numberofde)
#   # }
#   list1<-(totallength)
#   list0<-append(list1,list0)
#   namelist0<-append(tree$tip.label[i],namelist0)
# }
# # fg<-as.data.frame(list0)
# # ghyj<-as.data.frame(namelist0)
# phyloweight<-data.frame(namelist0=t(as.data.frame(namelist0)),list0=t(as.data.frame(list0)))
#
# #
# names(ghyj)<−NULL
# list0
# tr$tip.label
#
# list0[523]
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# myrmeFile <- system.file("iqmGHzaNTnW2h6iCBGlK-w_nexus.txt", package="phylobase")
# myrme <- readNewick(file="28TQ9-Ijo5Or5DnWVoE9rg_newick.txt")
# sed 's\; \|\g' 6uccD9engKXUPq-gxGXuYQ_newick.txt
#
# plot(tree)
# x <- as_tibble(tree)
# offspring(x, 2)
#
# library(ape)
#
# library(tidytree)
# data(bird.families)
# cophenetic(bird.families)['Eurostopodidae', 'Passeridae']
#
# treeman
# library(treeman)
# calcNdBlnc(tree)
# tree$edge.list
# plot(tree,show.node=TRUE)
# root <- rootNode(tree)
# library(adephylo)
# (require(ape) & require(phylobase))
# x <- as(rtree(20),"phylo4")
# tree <- as(tree,"phylo4")
# library(phytools)
# library(castor)
#
# .tipToRoot(tree, 526, root)
# distTips(tree, tips = c(1,526), method = "patristic")
#
# distRoot(tree, tips = 1, method = "sumDD")
#
# get_pairwise_distances(tree, c(1),c(526), as_edge_counts=FALSE, check_input=TRUE)
#
# Ntips = 100
# tree1 = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree
# # pick 3 random pairs of tips or nodes
# Npairs = 3
# A = sample.int(n=(Ntips+tree$Nnode), size=Npairs, replace=FALSE)
# B = sample.int(n=(Ntips+tree$Nnode), size=Npairs, replace=FALSE)
# distances = get_pairwise_distances(tree, c(1),c(526))
# tree = root_at_node(tree, new_root_node=325, update_indices=FALSE)
#
#
# filename <- "HW0EZAgBrbloeCRPDLSDmA_newick.txt"
# tree<-phytools::read.newick(filename)
#
# find_root(tree)
# root <- find_root(tree)
# tree1 <- as(tree,"phylo4")
#
# .tipToRoot(tree1,300, root)
# distances = get_pairwise_distances(tree, c(1),c(526))
# summary(tree1)
#
# .tipToRoot(tree1,527, root)
#
# list0<-c()
# for (i in (1:(find_root(tree)-1))){
#   totallength<-c()
#   for (j in (1:(length(.tipToRoot(tree1,i, root))))){
#     f<-c(i,.tipToRoot(tree1,i, root))
#     tiptonode<-get_pairwise_distances(tree, f[j],f[j+1])
#     numberofde<-length(descendants(tree, f[j+1], type = "daughter", ignore.tip = TRUE, labels = FALSE))
#     totallength<-append(totallength,tiptonode/numberofde)
#   }
#   list1<-sum(totallength)
#   names(list1)<-i
#   list0<-append(list1,list0)
# }
#
# list0<-as.data.frame(list0)
#
# get_pairwise_distances(tree, c(1),c(526))
#
# get_pairwise_distances
#
# get_distances_between_clades_CPP
#
# library(ape)
# dnd1 <- as.dendrogram(tree)
#
# myTree <- ape::read.tree(text='((A, B), ((C, D), (E, F)));')
#
# ggplot(tree) + geom_tree() + theme_tree()
# ggtree(tree)
