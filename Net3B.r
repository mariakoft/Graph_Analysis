setwd("C:/Users/Koft/Desktop/R Networks Project") ## Bale to path tou fakelou sto opoio exeis to dataset
Sys.setlocale("LC_CTYPE","Greek")
library(igraph)
library(Matrix)
library(lattice)
library(vcd)

#########Open Network File################
g<-read_graph(file="Erdos02.net", format="pajek")
summary(g)
V(g)$name[1:10]
is.simple(g) #TRUE άρα είναι κανονικό δίκτυο χωρίς λούπες και πολλαπλές ακμές και είμαι εντάξει
gorder(g)  # η τάξη του δικτύου : 6927
           # είναι το πλήθος των κορυφών #εναλλακτικά vcount(g)
gsize(g)   # ενώ το μέγεθος του δικτύου  :11850
           # είναι το πλήθος των ακμών # εναλλακιτκά ecount(g)

layg<-layout_with_fr(g) #εξάγω τις συντεταγμένες του δικτύου για να φτιάξω γράφημα
V(g)$coord1 <-layg[ ,1]
V(g)$coord2 <-layg[ ,2]    

## PLOT the GRAPH
E(g)$color <-  "darkgray"      
E(g)$width <- 1.5
V(g)$size<-3


plot(g, layout=layg, 
          edge.color="steelblue1", 
		      edge.width=0.2, 
			vertex.shape="circle",
		      vertex.frame.color="red1",
		      vertex.color="red1", 
		      vertex.size=0.5, 
			vertex.label=NA,             
	 		vertex.label.color="grey", xlim=c(-0.5, 0.8), ylim=c(-0.4,0.6)) 
			
	##Γενικά Χαρακτηριατικά
is.directed(g)						###FALSE
is.connected(g)                     ## TRUE
vertex.connectivity(g)             ##1
edge.connectivity(g) 				##1
graph.density(g)                     ##0.0004939929 << 0.25


###############################################
##########11. Τυχαίο Δικτυο και Σύγκριση#######
###################################################
erG <- erdos.renyi.game(6927, 11850, type= "gnm",
                        directed = FALSE, loops = FALSE)

is.simple(erG)                              ## OK !!
graph.density(erG)                          ## 0.0004939929  οσο και το αρχικό μας δίκτυο

layerG<-layout_nicely(erG)    ## εξαγωγή συντεταγμένων
V(erG)$coord1 <-layerG[ ,1]
V(erG)$coord2 <-layerG[ ,2]                   

erG
is.connected(erG)			#FALSE άρα ψάχνω για συνιστώσες

clerG<-clusters(erG); erG                      
clerG$no                                       ## 258 συνιστώσες
clerG$csize                                    ## 6650 κόμβοι στην 1η (γιγαντιαια) συνιστωσα
clerG$membership                               ## κομβοι ανα συνιστωσα

vgnerG<-vcount(erG)                             ## components' descriptives

##  the answers depent on each random graph
mcerG<-max(clerG$csize); mcerG ;mcerG/vgnerG    ## GC: 6650 nodes (96%) η γιγαντιαία συνιστώσα
mc2erG<-rev(sort(clerG$csize))[2]; mc2erG;mc2erG/vgnerG 	  ## 3 κόμβοι για την 2η σε σειρά συνιστώσα	 ## (0.0004330879%)  η πυκνοτητα της 2ης συνιστώσας

summary(clerG$csize)
##  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  1.00    1.00    1.00   26.85    1.00 6650.00 

## Για το αρχικό μου δίκτυο δεν εχω που ειναι συνδετικο, προφανως:
##Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 6927    6927    6927    6927    6927    6927 



##Φτιάχνω τα δύο γραφήματα των δικτύων δίπλα δίπλα για σύγκριση
par(mfrow = c(1, 2))                                              
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
			vertex.label=NA,xlim=c(-0.8,0.8), ylim=c(0.1,0.6))  
title(main = list("Erdos02", cex = 1.5,
                  col = "white", font = 3))

g2<-erG                 
E(g2)$color <-  "darkgray"      
E(g2)$width <- 1.5
V(g2)$size<-3

plot(g2, layout=layerG, 
                  edge.color="steelblue1", 
		      edge.width=1, 
			vertex.shape="circle",
		      vertex.frame.color="red1",
		      vertex.color="red1", 
		      vertex.size=1.5, 
			vertex.label=NA,xlim=c(-0.8,0.8), ylim=c(0.1,0.6))  
par(mfrow = c(1, 1)) 
title(main = list("        Τυχαίο", cex = 1.5,
                  col = "white", font = 3))


## HSI - > inhomogeneous net with hubs   (ανομοιογένεια και ομφαλοί)

## ER  - >  homogeneous net, almost connected (ομοιογένεια σχεδόν συνδετικό)



##########################################################################
# Κατανομή βαθμών Τυχαίου Δικτύου και Σύγκριση#######
########################################################################## 
erGd <- degree(erG)
summary(erGd )  #τυχαίο δίκτυο                
##Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##0.000   2.000   3.000   3.421   5.000  12.000 
summary(degree(g)) #το δικτυο μου
##Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##1.000   1.000   1.000   3.421   2.000 507.000 
  
# standard deviation
sd(erGd )		#1.840168 τυχαιο δικτυο
sd(degree(g))	#11.68919	το δικτυο μου
# coefficient of variation
cverG <-sd(erGd )/mean(erGd );cverG    ## 0.5378415%  
#make a data frame with degree distribution
erGdist <-table(erGd)           
vert <-as.integer(names(erGdist ))    ## names are the number of veertices
erGdistdf<- data.frame( "degrees" = vert, 
                        "number of vertices" = as.vector(erGdist ) )
erGdistdf                               ## frequency distribution

## Γραφικη Αναπαράσταση και Έλεγχος Υποθεσης για την κατανομη των βαθμων του τυχαίου

derG<-degree.distribution(erG, cumulative = FALSE)   
par(mfrow=c(1,2))
plot(0:max(degree(erG)),derG, xlab="degree", ylab="freq");
hist(erGd ,  freq=FALSE, right = FALSE, xlab="degree", ylab="freq") 

par(mfrow=c(1,1))
gfp<-goodfit(erGd ,type= "poisson",method= "ML")
summary(gfp)
##Goodness-of-fit test for poisson distribution

             ##         X^2 df  P(> X^2)
##Likelihood Ratio 9.450215 11 0.5804094 
# Δεν απορρίπτεται η υπόθεση ό,τι η κατανομή των βαθμών στο τυχαίο δίκτυο είναι Poissom
plot(gfp,main="Count data vs Poisson distribution")


#############################################################
################ 4. Ιδιότητα Μικροκόσμου ###################
##############################################################
# Η μελέτη της ιδιότητας small world βασίζεται στον υπολογισμό
#της διαμέτρου, της μέσης απόστασης αλλά και του συντελεστή σύμπλεξης,
# σε σύγκριση με τις τιμές ενός  ER 
## DIAMETER
dg<-diameter(g, unconnected = TRUE, weights = NA) ; dg   ## 4 (using the whole net)
## DIAMETER in ER
dgerG<-diameter(erG, unconnected = TRUE, weights = NA) ; dgerG   ## 15

## AVERAGE PATH LENGTH  (Average distance)
aplg0<-average.path.length(g, directed=FALSE, unconnected=TRUE)  ; aplg0 ## 3.776441 (the whole net)
## AVERAGE PATH LENGTH in ER
aplg0erG<-average.path.length(erG, directed=FALSE, 
                              unconnected=TRUE)  ; aplg0erG ## 7.250466
## δεδομένου ότι το τυχαίο δίκτυο έχει την ιδιότητα μικροκόσμου και οι αποστάσεις στο κανονικό δίκτυο είναι ακόμη πιο μικρές,
## φαίνεται πολύ έντονα η ιδιότητα μικροκόσμου στο δίκτυο μας

## Κατανομή των αποστάσεων
splH<-path.length.hist (g, directed = FALSE) 

## sum(splH$res);409*408/2 these are all distances

##      1<=DISTANCE<=4 : diameter
plot(1:4, splH$res, xlab="distance", ylab="paths",frame=TRUE)
md<-sum(splH$res %*% 1:4)/splH$res; md    
## or average path length
## using theory
mdd<-median(rep(1:4, splH$res));mdd  ## compute the median=4
mdsd<-sd(rep(1:4, splH$res));mdsd    ## compute the SD=0.4616965
mdcv<-(mdsd/md)*100;mdcv              ## compute the coefficient of variation=0.04252357


# Ελέγχω Συντελεστές Σύμπλεξης και Μεταβατικότητα 
transitivity( g, type="global", isolates= "zero")    ## 0.03570461

transitivity( g, type="average", isolates= "zero")   ## 0.1239001 (watts and strogatz)

##  ER(6927, 11850)
dgerG<-diameter(erG, unconnected = TRUE, weights = NA) ; dgerG   ## 15

transitivity( erG, type="global", isolates= "zero")    ##  0.0004453241 
transitivity( erG, type="average", isolates= "zero")   
##  0.0002977766  (watts and strogatz)
## expected value :  2*11850/(6927*(6927-1))  =  0.0004939929

##Επαληθεύεται η ιδιότητα του μιρκοκόσμου

########################################################################
############# 7. Υπαρξη Μοτίβων ###################################
#####################################################################

#Αναζητώ την ύπαρξη υπογραφημάτων 3 ή 4 κορυφών, η οποία μου δίνει σημαντικές πληροφορίες για τη δομή και την εξέληξη του

graph.motifs(g,size=3) #Ερευνώ για τρίγωνα στο δίκτυο : μη συνδετικα=NA / με μια ακμη= NA / με 2 ακμες= 483949 / τρίγωνα= 5973
graph.motifs(g,size=4) #NA       NA       NA       NA  	33000237       NA  8445001  1420445    13999  46989     3323
### USE FANNMOD   (με έλεγχο υπόθεσης για το πλήθος των Κ3 ή Κ4)
write.graph(g, "FOR_FANMOD ERDOS.txt", format=c("edgelist")) #εξαγωγή αρχείου ακμών για να το εισάγω στο fanmod

################################################################################
############### 9. Μέτρα Κεντρικότητας #########################
#########################################################################

cDg<-degree(g)                         ## Βαθμική Κεντρικότητα

cEg1<-round(evcent 
            (g, scale = TRUE)$vector,6)    ## Ιδιοκεντρικότητα (max=1)

cEg<-round(evcent 
           (g, scale = FALSE)$vector,6)   ## Ιδιοκεντρικοτητα μη κανονικοποιημένη 

cCg<-round(closeness(g),6)              ## Εγγύτητα

cBg<-round(betweenness
           (g, directed = FALSE),6)          ## Διαμεσότητα



##  Δημιουργία ενός data frame με τα μέτρα κεντρικότητας ανά κορυφή

nodesg<-as.numeric(V(g)) 
componentg<-as.numeric(clusters(g)$membership)

centralitiesg<-data.frame(x1=nodesg, x2=componentg, x3=cDg,
                          x4=cEg1, x5=cCg, x6=cBg)
dimnames(centralitiesg)[[2]]<-c("node" , "component", "degree",  
                                "eigenvC",  
                                "closeness", "betweenness")

str(centralitiesg)                         ## variables in the data frame


## Εξαγωγή των 3 κορυφών με τη μεγαλύτερη βαθμική κεντρικότητα
ord.deg.centralitiesg<-centralitiesg[order(-centralitiesg$degree), ][1:3,];  
ord.deg.centralitiesg  ## 
# 		 node component degree  eigenvC closeness betweenness
#6927 6927         1    507 1.000000   7.5e-05  19352996.8
#186   186         1    297 0.300845   5.3e-05   1389653.4
#10     10         1    215 0.367622   5.4e-05    821278.5   

## Εξαγωγή των 3 κορυφών με τη μεγαλύτερη ιδιοκεντρικότητα

ord.deg.centralitiesg1<-centralitiesg[order(-centralitiesg$eigenvC), ][1:3,];  
ord.deg.centralitiesg1
#node component degree  eigenvC closeness betweenness
#6927 6927         1    507 1.000000   7.5e-05  19352996.8
#10     10         1    215 0.367622   5.4e-05    821278.5
#186   186         1    297 0.300845   5.3e-05   1389653.4

## Εξαγωγή των 3 κορυφών με τη μεγαλύτερη εγγύτητα
ord.deg.centralitiesg4<-centralitiesg[order(-centralitiesg$closeness), ][1:3,];  
ord.deg.centralitiesg4
# node component degree  eigenvC closeness betweenness
#6927 6927         1    507 1.000000   7.5e-05  19352996.8
#10     10         1    215 0.367622   5.4e-05    821278.5
#164   164         1    132 0.286753   5.4e-05    503964.7

## Εξαγωγή των 3 κορυφών με τη μεγαλύτερη διαμεσότητα
ord.deg.centralitiesg5<-centralitiesg[order(-centralitiesg$betweenness), ][1:3,];  
ord.deg.centralitiesg5
#	node component degree  eigenvC closeness betweenness
#6927 6927         1    507 1.000000   7.5e-05  19352996.8
#186   186         1    297 0.300845   5.3e-05   1389653.4
#10     10         1    215 0.367622   5.4e-05    821278.5

#Δημιουργία γραφημάτων  
layg0<-layg                      
gcC<-g                            
par(bg="black", new=FALSE)                   ## vertex size analog to  
Meiggc<-max(centralitiesg$eigenvC)                    ## eigenvector centrality
V(gcC)$size<-(centralitiesg$eigenvC/Meiggc)*8
plot.igraph(gcC, layout=layg0, 
            #asp=1, margin=0.1,
            edge.color="lightseagreen", edge.width=0.05, 
            vertex.shape="circle",
            vertex.frame.color="lightblue",
            vertex.color="snow",  
            vertex.size=V(gcC)$size,  vertex.label=NA,xlim=c(-0.9,0.9), ylim=c(-0.3,0.8))
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
            vertex.size=V(gcC)$size,  vertex.label=NA,xlim=c(-0.9,0.9), ylim=c(-0.3,0.8))
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
            vertex.size=V(gcC)$size,  vertex.label=NA,xlim=c(-0.9,0.9), ylim=c(-0.3,0.8))
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
            vertex.size=V(gcC)$size,  vertex.label=NA,xlim=c(-0.9,0.9), ylim=c(-0.3,0.8))
title(main = list("Κεντικότητα Εγγύτητας", cex = 1.5,
                  col = "white", font = 3))

############## Μελέτη της γειτονιάς του Erdos μιας και ειναι ο κεντρικότερος κομβος
gcC<-g
iniv<-subset(centralitiesg, 
			centralitiesg$eigenvC == 
				max(centralitiesg$eigenvC))$node
iniv<-V(gcC)[iniv]

neieigenvC<-V(gcC)[ nei( iniv)]
neieigenvC
eisub<-induced.subgraph(gcC,c(iniv,neieigenvC))
eisub



par(mfrow = c(1, 1))
laysub<-cbind(V(eisub)$coord1,
		 V(eisub)$coord2)

eisub
par(bg="black",new=FALSE)                          
eisubplot<-plot.igraph(eisub,   
		layout=laysub, 
		vertex.color="snow" ,                                 
			vertex.size=0.3, 
			edge.color =  "lightseagreen", 
			edge.width=0.2, 
			vertex.label=V(eisub)$names, 
			vertex.label.color="white", 
			vertex.label.dist=0.5, vertex.label.font=1, 
			vertex.label.cex=0.4, ylim=c(-0.4,0.5)) 

				  
##########################################################################################
######################## 10. Κλίκες και bi-components#######################
####################################################################################

# εύρεση των κλικών του δικτύου
table(sapply(cliques(g), length))
 #1     2     3     4     5     6     7     8 
 #6927 11850  5973  3323  1402   391    64     5 
#ο δίκτυο μου είναι πολύ μεγάλο οπότε 
#θα μάθω περισσότερα και πιο ουσιώδη από την αναζήτηση bi-conmponents
bc <- biconnected_components(g)  #  εύρεση bi-connected components 
                                 # στην igraph υπολογίζονται και οι ακμές 
bc$no	#4782
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
#3/6927 vertices:
#  507  156 6927
