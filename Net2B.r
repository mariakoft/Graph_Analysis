setwd("C:/Users/Koft/Desktop/R Networks Project") ## Bale to path tou fakelou sto opoio exeis to dataset
Sys.setlocale("LC_CTYPE","Greek")
library(igraph)
library(Matrix)
library(lattice)
library(vcd)

#########Open Network File################
listg<-read.table("CA-GrQc.txt", header=TRUE)
g0<-graph_from_data_frame(listg,dir=FALSE)
summary(g0)
##########################################################
#################1. ενικές Πληροφορίες Δικτυου############
############################################################
is.simple(g0)  #FALSE
g<-simplify(g0)   # τώρα το δίκτυο έγινε απλό και ξαναπερνάω στους προηγούμενους ελέγχους
gorder(g0)  # η τάξη του δικτύου#εναλλακτικά vcount(g)
           # είναι το πλήθος των κορυφών 5242
gsize(g0)   # ενώ το μέγεθος του δικτύου #εναλλακτικά ecount(g)
           # είναι το πλήθος των ακμών 28980
is.simple(g0)  # βλέπω αν το δίκτυο είναι απλο (δεν υπάρχουν πολλαπλές ακμές, λούπες)
              # FALSE άρα δεν είναι και θα το κάνω με την επόμενη εντολή.
  # Μάλιστα ήταν αναμενόμενο το FALSE αφού κάποιο επιστήμονες μπορεί να συνεργάστηκαν 
              # περισσότερες από μία φορές
gsize(g)     # η τάξη του δικτύου, οόπως είναι και λογικό άλλωστε παραμένει η ίδια 
             # όμως το μέγεθος του δικτύου είναι τώρα 14484
gsize(g)/gsize(g0)       # δηλαδή υπήρξαν 14484 συνεργασίες μεταξύ δύο διαφορετικών ατόμων 49.9793 % του δικτυου εφυγε
is.weighted(g)   #FALSE
is.named(g)		#TRUE
is.directed(g)  #FALSE  προφανώς όχι λόγο κατασκευής
is.connected(g) # FALSE άρα το δίκτυό μου είναι μη συνδετικό
                # άρα πολλές φορές αντί να δουλεύω σε ολόκληρο το δίκτυο θα δουλεύω σε συνιστώσες
graph.density(g)  # Η πυκνότητα του δικτύου είναι 0.001054405%, θεωρήται αραιό
		   
layg<-layout_nicely(g) 
V(g)$coord1 <-layg[ ,1]
V(g)$coord2 <-layg[ ,2]    

## PLOT the GRAPH                   
E(g)$color <-  "darkgray"      
E(g)$width <- 1.5
V(g)$size<-3

plot(g, layout=layg, 
                  edge.color="steelblue1", 
		      edge.width=1, 
			vertex.shape="circle",
		      vertex.frame.color="red1",
		      vertex.color="red1", 
		      vertex.size=1.5, 
			vertex.label=NA)  
##########################################################################
#
# 3.2 Plot the network so as to distinguish the components
#        Διαφορετικό χρώμα ανά συνιστώσα
#
##########################################################################

			##########################################################################
			#
			# 3.1  Ανάλυση Συνιστωσών και εύρεση γιγαντιαίας
			#
			##########################################################################


			clg<-clusters(g); clg                           ## components info
			clg$no                                          ## 355 components
			clg$csize                                       ## size of each component
			clg$membership                                  ## nodes in each component

			vgn<-vcount(g)                                  ## components' descriptives
			mc<-max(clg$csize); mc ;mc/vgn                  ## GC: 4158 nodes (79.32%)
			mc2<-rev(sort(clg$csize))[2]; mc2; mc2/vgn      ## 2nd order component
								## 14 nodes (0.02%)
			sort(clg$csize)
#Σχεδιασμός της γιγαντιαίας

g1<-g                                              ## drawing components
comps1<-clg$membership
colbar<-rainbow(max(comps1))
par(bg="black",new=FALSE)                           ## black backround
V(g1)$color <- colbar[comps1]
plot(g1, layout=layg, 
            edge.color="lightsteelblue1", 
		edge.width=0.1, vertex.shape="circle",
            vertex.frame.color=V(g1)$color,
            vertex.color=V(g1)$color,  vertex.size=1.5,  
		vertex.label=NA)




pdf("component_drawing.pdf")       ## save in pdf file format

g1<-g                                              
comps1<-clg$membership
colbar<-rainbow(max(comps1))
par(bg="black",new=FALSE)                           
V(g1)$color <- colbar[comps1]
plot(g1, layout=layg, 
            edge.color="lightsteelblue1", 
		edge.width=0.1, vertex.shape="circle",
            vertex.frame.color=V(g1)$color,
            vertex.color=V(g1)$color,  vertex.size=1.5,  
		vertex.label=NA)
dev.off()

   
		   
		   
###############################################
##########11. Τυχαίο Δικτυο και Σύγκριση#######
###################################################
erG <- erdos.renyi.game(5242, 28980, type= "gnm",
                        directed = FALSE, loops = FALSE)

is.simple(erG)                              ## OK !!
graph.density(erG)                          ## 0.002109683

layerG<-layout_nicely(erG)    ## get the coordinates
V(erG)$coord1 <-layerG[ ,1]
V(erG)$coord2 <-layerG[ ,2]                   ## use them for subgraphs

erG

clerG<-clusters(erG); erG                      ## components info
clerG$no                                       ## 1 συνιστώσα-άρα συνδετικό

vgnerG<-vcount(erG)                             ## components' descriptives


summary(clerG$csize)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   5242    5242    5242    5242    5242    5242  προφανως αφού είναι συνδετικο
## The real net
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    2.00    2.00   14.77    3.00 4158.00  



## Δημιουργία και των 2 γραφημάτων για σύγκριση

par(mfrow = c(1, 2))                                              
comps1<-clg$membership
colbar<-rainbow(max(comps1))
par(bg="black",new=FALSE)                           ## heterogeneity?
V(g1)$color <- colbar[comps1]
plot(g1, layout=layg, 
     edge.color="lightsteelblue1", 
     edge.width=0.1, vertex.shape="circle",
     vertex.frame.color=V(g1)$color,
     vertex.color=V(g1)$color,  vertex.size=1.5,  
     vertex.label=NA)
title(main = list("GrQc", cex = 1.5,
                  col = "white", font = 3))


g2<-erG                                             
comps2<-clerG$membership
colbar<-rainbow(max(comps2))
par(bg="black",new=FALSE)                           ## homogeneity?
V(g2)$color <- colbar[comps2]
plot(g2, layout=layerG, 
     edge.color="lightsteelblue1", 
     edge.width=0.1, vertex.shape="circle",
     vertex.frame.color=V(g2)$color,
     vertex.color=V(g2)$color,  vertex.size=1.5,  
     vertex.label=NA)
par(mfrow = c(1, 1)) 
title(main = list("Τυχαίο", cex = 1.5,
                  col = "white", font = 3))




## HSI - > inhomogeneous net with hubs   (ανομοιογένεια και ομφαλοί)

## ER  - >  homogeneous net, almost connected (ομοιογένεια σχεδόν συνδετικό)



##########################################################################
# Κατανομή βαθμών τυχαίου δικτύου
########################################################################## 

erGd <- degree(erG)
summary(erGd )      
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    9.00   11.00   11.06   13.00   24.00        
summary(degree(g))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #0.000   2.000   3.000   5.526   6.000  81.000
# standard deviation
sd(erGd ) #3.359724
sd(degree(g)) #7.918456
# coefficient of variation
cverG <-sd(erGd )/mean(erGd );cverG    ## 3.038591%
#make a data frame with degree distribution
erGdist <-table(erGd)           
vert <-as.integer(names(erGdist ))    ## names are the number of veertices
erGdistdf<- data.frame( "degrees" = vert, 
                        "number of vertices" = as.vector(erGdist ) )
erGdistdf                               ## frequency distribution

## plot the degree distribution and do a Hypothesis test

derG<-degree.distribution(erG, cumulative = FALSE)   
par(mfrow=c(1,2))
plot(0:max(degree(erG)),derG, xlab="degree", ylab="freq");
hist(erGd ,  freq=FALSE, right = FALSE, xlab="degree", ylab="freq") 

par(mfrow=c(1,1))
gfp<-goodfit(erGd ,type= "poisson",method= "ML")
summary(gfp)
#Goodness-of-fit test for poisson distribution
 #                     X^2 df  P(> X^2)
#Likelihood Ratio 17.13388 22 0.7559296

plot(gfp,main="Count data vs Poisson distribution")

#############################################################
################ 4. Ιδιότητα Μικροκόσμου ###################
##############################################################
# Η μελέτη της ιδιότητας small world βασίζεται στον υπολογισμό
# της διαμέτρου, της μέσης απόστασης αλλά και του συντελεστή σύμπλεξης,
# σε σύγκριση με τις τιμές ενός δείγματος ER τυχαίων δικτύων 
clg$membership 

gcg0<-induced.subgraph(g,V(g)[clg$membership == 1])  
## nodes in GC noted with 1                          

layg0<-cbind(V(gcg0)$coord1,
             V(gcg0)$coord2)

plot(gcg0, layout=layg0, 
     edge.color="steelblue1",              
     edge.width=1,                         
     vertex.shape="circle",
     vertex.frame.color="red1",
     vertex.color="red1", 
     vertex.size=2, 
     vertex.label=NA,
     vertex.label.dist=0.2, 
     vertex.label.font=4, 
     vertex.label.cex=0.5, 
     vertex.label.color="red")

## DIAMETER
dg<-diameter(g, unconnected = TRUE, weights = NA) ; dg   ## 17 (using the whole net)
dgcg0<-diameter(gcg0, unconnected = FALSE, weights = NA) ; dgcg0 ## 17 (using the GC)


## DIAMETER in ER
dgerG<-diameter(erG, unconnected = TRUE, weights = NA) ; dgerG   ## 6

## AVERAGE PATH LENGTH  (Average distance)
aplg0<-average.path.length(g, directed=FALSE, unconnected=TRUE)  ; aplg0  ## 6.048515 (the whole net)
aplgcg0<-average.path.length(gcg0,directed=FALSE, 
                             unconnected=FALSE) ; aplgcg0 ## 6.04938(the GC-that will be used)

## AVERAGE PATH LENGTH in ER
aplg0erG<-average.path.length(erG, directed=FALSE, 
                              unconnected=TRUE)  ; aplg0erG ##3.819717

## make the distribution of distances
splH<-path.length.hist (gcg0, directed = FALSE) 


## sum(splH$res);409*408/2 these are all distances

##      1<=DISTANCE<=17 : diameter
plot(1:17, splH$res, xlab="distance", ylab="paths",frame=TRUE)
md<-sum(splH$res %*% 1:17)/(409*408/2); md   ## compute mean distance=626.6022
## or average path length
## using theory
mdd<-median(rep(1:17, splH$res));mdd  ## compute the median=6
mdsd<-sd(rep(1:17, splH$res));mdsd    ## compute the SD=1.570
mdcv<-(mdsd/md)*100;mdcv              ## compute the coefficient of variation=0.25


# 8.3    Transitivity - Clustering Coefficient THE SMALL WORLD property
#  Η τιμή του συντελεστή σύμπλεξης είναι αρκετά πιο μικρή σε 
#  πραγματικό δίκτυο από ότι σε αντίστοιχο ER τυχαίο δίκτυο

transitivity( g, type="global", isolates= "zero")    ##  0.6298425

transitivity( g, type="average", isolates= "zero")   ## 0.5296358  (watts and strogatz)

##  SMALL WORLD
## high value compared to that of the 
## ER random network

##  ER
dgerG<-diameter(erG, unconnected = TRUE, weights = NA) ; dgerG   ## 6
##  very close  (SMALL WORLD)

transitivity( erG, type="global", isolates= "zero")    ## 0.002196077 

transitivity( erG, type="average", isolates= "zero")   ##  0.002196077 (watts and strogatz)
## expected value :  2*28980/(5242*(5242-1))  =  0.002109683

########################################################################
############# 7. Υπαρξη Μοτίβων ###################################
#####################################################################

#Αναζητώ την ύπαρξη υπογραφημάτων 3 ή 4 κορυφών, η οποία μου δίνει σημαντικές πληροφορίες για τη δομή και την εξέληξη του

graph.motifs(g,size=3) #Ερευνώ για τρίγωνα στο δίκτυο : μη συνδετικα=NA / με μια ακμη= NA / με 2 ακμες= 85087 / τρίγωνα= 48260
graph.motifs(g,size=4)  #NA     NA     NA     NA 405750     NA 553322 628366   1115  65717 329297

### USE FANNMOD   (με έλεγχο υπόθεσης για το πλήθος των Κ3 ή Κ4)

write.graph(g, "FOR_FANMOD CA.txt", format=c("edgelist"))
## then open FANNMOD and run the commands in there

################################################################################
############### 9. Μέτρα Κεντρικότητας #########################
#########################################################################
cDg<-degree(gcg0)                         ## degree centrality

cEg1<-round(evcent 
            (gcg0, scale = TRUE)$vector,6)    ## eigenvector centrality (max=1)

cEg<-round(evcent 
           (gcg0, scale = FALSE)$vector,6)   ## eigenvector centrality 
## (max=1)

cCg<-round(closeness(gcg0),6)              ## closeness centrality

cBg<-round(betweenness
           (gcg0, directed = FALSE),6)          ## betweenness centrality



##  Δημιουργία ενός data frame με τα μέτρα κεντρικότητας ανά κορυφή

nodesg<-as.numeric(V(gcg0)) 
componentg<-as.numeric(clusters(gcg0)$membership)

centralitiesg<-data.frame(x1=nodesg, x2=componentg, x3=cDg,
                          x4=cEg1,  x5=cCg, x6=cBg)
dimnames(centralitiesg)[[2]]<-c("node" , "component", "degree", 
                               "eigenvC", "closeness", "betweenness")

str(centralitiesg)                         ## variables in the data frame

## Get centrality measures for GIANT COMPONENT
gcdata.f<-subset(centralitiesg, 
                 centralitiesg$component == 1,) 


## Get centrality measures for GIANT COMPONENT
gcdata.f<-subset(centralitiesg, 
                 centralitiesg$component == 1,) 

## Degree Centrality first 3 nodes
ord.deg.centralitiesg<-centralitiesg[order(-centralitiesg$degree), ][1:3,];  
ord.deg.centralitiesg  ## 
#  		node component degree 
#21012   59         1     81        
#21281  722         1     79        
#22691   60         1     77        

## eigenvC Centrality first 3 nodes
ord.deg.centralitiesg1<-centralitiesg[order(-centralitiesg$eigenvC), ][1:3,];  
ord.deg.centralitiesg1
#		node    eigenvC 
#21012   59     1.000000 
#2741   570     0.987224 
#12365 1488     0.983995 

## closeness Centrality first 3 nodes
ord.deg.centralitiesg4<-centralitiesg[order(-centralitiesg$closeness), ][1:3,];  
ord.deg.centralitiesg4
# 		node component  closeness 
#13801 1314         1    5.9e-05    
#2654    55         1     5.7e-05 
#21012   59         1    5.7e-05    
> 

## betweenness Centrality first 3 nodes
ord.deg.centralitiesg5<-centralitiesg[order(-centralitiesg$betweenness), ][1:3,];  
ord.deg.centralitiesg5
#		node  betweenness
#13801 1314     508435.4
#9572    65     352746.5
#14599  903     349992.2

#Δημιουργία γραφημάτων   
layg0<-layg                       
gcC<-gcg0                            ## another drawing of the GC                         
par(bg="black", new=FALSE)                   ## vertex size analog to  
Meiggc<-max(centralitiesg$eigenvC)                    ## eigenvector centrality
V(gcC)$size<-(centralitiesg$eigenvC/Meiggc)*8
plot.igraph(gcC, layout=layg0, 
            #asp=1, margin=0.1,
            edge.color="lightseagreen", edge.width=0.05, 
            vertex.shape="circle",
            vertex.frame.color="lightblue",
            vertex.color="snow",  
            vertex.size=V(gcC)$size,  vertex.label=NA,xlim=c(-0.9,0.9), ylim=c(-0.5,0.8))
title(main = list("Ιδιοκεντρικότητα", cex = 1.5,
                  col = "white", font = 3))


par(bg="black",new=FALSE)                    ## vertex size analog to 
Meiggc<-max(centralitiesg$degree)                      ## degree centrality
V(gcC)$size<-(centralitiesg$degree/Meiggc)*8
plot.igraph(gcC, layout=layg0, 
            asp=1, margin=0.1,
            edge.color="lightseagreen", edge.width=0.05, 
            vertex.shape="circle",
            vertex.frame.color="lightblue",
            vertex.color="snow",  
            vertex.size=V(gcC)$size,  vertex.label=NA,xlim=c(-0.9,0.9), ylim=c(-0.5,0.8))
title(main = list("Βαθμική Κεντρικότητα", cex = 1.5,
                  col = "white", font = 3))




par(bg="black",new=FALSE)                      ## vertex size analog to     
Meiggc<-max(centralitiesg$betweenness)                   ## betweenness centrality
V(gcC)$size<-(centralitiesg$betweenness/Meiggc)*8
plot.igraph(gcC, layout=layg0, 
            asp=1, margin=0.1,
            edge.color="lightseagreen", edge.width=0.05, 
            vertex.shape="circle",
            vertex.frame.color="lightblue",
            vertex.color="snow",  
            vertex.size=V(gcC)$size,  vertex.label=NA,xlim=c(-0.9,0.9), ylim=c(-0.5,0.8))
title(main = list("Διαμεσότητα", cex = 1.5,
                  col = "white", font = 3))


par(bg="black",new=FALSE)                      ## vertex size analog to  
Meiggc<-max(centralitiesg$closeness)                     ## closeness centrality
V(gcC)$size<-(centralitiesg$closeness/Meiggc)*5
plot.igraph(gcC, layout=layg0, 
            asp=1, margin=0.1,
            edge.color="lightseagreen", edge.width=0.05, 
            vertex.shape="circle",
            vertex.frame.color="lightblue",
            vertex.color="snow",  
            vertex.size=V(gcC)$size,  vertex.label=NA,xlim=c(-0.9,0.9), ylim=c(-0.5,0.8))
title(main = list("Κεντικότητα Εγγύτητας", cex = 1.5,
                  col = "white", font = 3))


##########################################################################################
######################## 10. Κλίκες και bi-components#######################
####################################################################################

# εύρεση των κλικών του δικτύου
table(sapply(cliques(g), length))

#ο δίκτυο μου είναι πολύ μεγάλο οπότε 
#θα μάθω περισσότερα και πιο ουσιώδη από την αναζήτηση bi-conmponents
bc <- biconnected_components(g)  #  εύρεση bi-connected components 
                                 # στην igraph υπολογίζονται και οι ακμές 
bc$no	#1548
bc$components[]
# bi-connected components με τουλάχιστον 3 κορυφές
bic<-list()
k<-0
for (i in 1:bc$no) {
  if ( length(bc$components[[i]])  >= 3 ) {
    k<-k+1
    bic[k]=bc$components[i]
  }
}
bic
bil<-list()
k<-0
for (i in 1:bc$no) {
  if ( length(bc$components[[i]])  > 3 ) {
    k<-k+1
    bil[k]=bc$components[i]
  }
}
bil
#υπάρχουν πολυ bi-components τους περνάω κατευθείαν στην παρουσίαση
