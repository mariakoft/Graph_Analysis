setwd("~/R/win-library/3.4") ## βάλε  το path του φακέλου που έχεις το dataset 
library(igraph)
library(Matrix)
library(lattice)
library(vcd)
#########Open Network File################
name1<-read.table("country.txt", header = TRUE)
name1[1:3,]
> mat1<-read.table("ws0 - Food and live animals.txt", header = FALSE)
> mat2<-as.matrix(mat1)
isSymmetric(mat2) # mode = "directed"
[1] FALSE
colnames(mat2)<-name1$name
> net1<- graph_from_adjacency_matrix(
      mat2, mode = "directed",
     weighted = NULL,
     diag = FALSE,
     add.colnames = NULL,
     add.rownames = NA)

##########################################################
#
#	1. Γενικές Πληροφορίες Δικτυου
#
############################################################
is.simple(net1) #TRUE άρα το δίκτυο είναι απλό, δηλαδλή δεν έχουμε λούπες
is.connected(net1) #TRUE άρα το δίκτυο είναι συνδεδεμένο
is.directed(net1) #TRUE άρα το δίκτυο είναι κατευθυνόμενο
is.weighted(net1) #FALSE άρα το δίκτυο είναι μη σταθμισμένο
is.named(net1) #TRUE
vcount(net1)#24 υπολογίζουμε τον αριθμό των κορυφών
is_bipartite(net1)  # αν υπάρχει attribute: "type" για τις κορυφές
girth(as.undirected(net1, edge.attr.comb = "ignore")) #μήκος μικρότερου κύκλου
ecount(net1) #307 υπολογίζουμε τον αριμό των ακμών
vertex.connectivity(net1)#1 υπολογίζουμε τον συνδετικό αριθμό των κορυφών, που σημαίνει ότι 
							ο ελάχιστος αριθμός κόμβων που αν διαγραφεί κάνει το δίκτυο μη-συνδετικό είναι ένας κόμβος
edge.connectivity(net1)#1 υπολογίζουμε τον συνδετικό αριθμό των ακμών, που σημαίνει ότι αν αφαιρεθεί 1  ακμή, 
						  το δίκτυο γίνεται μη-συνδετικό
graph.density(net1)#0.5561594 υπολογίζουμε την πυκνότητακνότητα του δικτύου, άρα το δίκτυο΄ είναι πολύ πυκνο
	plot(net1)
V(net1)$color<-"red1"
E(net1)$color <- "steelblue1"
E(net1)$width <- 1
V(net1)$size<-5
E(net1)$arrow.size=2  #αναπαράσταση αρχικού δικτύου

##########################################################################
#
#   2. Ακολουθίες βαθμών έσω-έξω
#
##########################################################################
outD<-degree(net1, mode="out")
outD #ακολουθία των έξω βαθμών του κάθε κόμβου
inD<-degree(net1, mode="in")
> inD #ακολουθία των έσω βαθμών του κάθε κόμβου
sum(outD)==sum(inD)
##########################################################################
#
#   3. Εύρεση αμιβαίων σχέσεων
#
##########################################################################
mutnet1<-which_mutual(net1)
mutnet1 #εύρεση αμιβαίων σχέσεων
net1wm<-net1-edges(E(net1)[mutnet1])
net1wm #υπογράφημα με αφαίρεση ακμών που δεν είναι αμοιβαίες σχέσεις
net1M<-simplify(as.undirected(net1,mode="mutual"),remove.multiple = TRUE)
summary(net1M)# το δίκτυο των αμιβαίων σχέσεων
summary(net1) # το αρχικό δίκτυο
plot(net1M)
V(net1M)$color<-"red1"
E(net1M)$color <- "steelblue1"
E(net1M)$width <- 1
V(net1M)$size<-1.5 #αναπαράσταση δικτύου των αμιβαίων σχέσεων
cn<-clique_num(net1M)
cn #έρευνα για μέγιστη κλίκα, υπολογίζουμε τον αριμό τοων μόμβων τις κλίκας
Vmc<-cliques(net1M, min=cn, max=cn)
Vmc# λίστα με τα ονόματα των κόμβων που ανήκουν στην κλίκας
#Φτιάχνουμε υπογράφημα που δεν περιέχει τις κλίκες στο δίκτυο 
cq<-as_ids(Vmc[[1]])
cq#παίρνουμε τις κορυφές τις κλίκας
neicq<-as_ids(neighbors(net1, cq,mode="all"))#γείτονες τις κλίκας
sec<-unique(neicq)#πέρνουμε κάθε άτομο μία φορά
net1T<-net1-vertices(cq) #αφαιρούμε τις κλίκες από το αρχικό δίκτυο
net1T<-net1M-vertices(cq) #αφαιρούμε τις κλίκες από το  δίκτυο αμοιβα σχέσεων
plot(net1T) #αναπαράσταση δικτύου χωρίς την παρουσία τις κλίκας
##########################################################################
#
#   4. Αποστάσεις στο δικτυο
#
##########################################################################
> gnet1<-net1
> gnet1lay <- layout_with_fr(gnet1)  # συντεταγμένες κορυφών
> V(gnet1)$x <- gnet1lay [ ,1]    # μόνιμο χαρακτηριστικό των κορυφών (x)
> V(gnet1)$y <- gnet1lay [ ,2]    # μόνιμο χαρακτηριστικό των κορυφών (y)
> V(gnet1)$color<- rainbow(7,alpha=.9)[5]
> E(gnet1)$arrow.size=0.2
> E(gnet1)$color="darkgray"
> E(gnet1)$distance<-1/(0.5+E(gnet1)$weight) #απόσταση αντόστροφη του βάρους
> summary(gnet1) # ποια νέα attributes προστέθηκαν
 plot(gnet1 , #layout=gFlay,
      vertex.shape="circle",
      vertex.frame.color= V(gnet1)$color,
      vertex.size=10,
      vertex.label=V(gnet1)$name,
      vertex.label.font=2,
      vertex.label.cex=1,
      vertex.label.dist=2,
      vertex.label.degree=-pi/2) #αναπαράσταση δικτύου με νέα χαρακτηριστικά
Dtgnet11<-distances(gnet1, v = V(gnet1), mode = "out", weights =
                      E(gnet1)$distance,
                  algorithm = "automatic")	  #πίνακας αποστάσεων out με βάρη
Dtgnet12<-distances(gnet1, v = V(gnet1), mode = "in", weights =
                      E(gnet1)$distance,
                  algorithm = "automatic") #πίνκας αποστάσεων in με βάρη
Dtgnet11==t(Dtgnet12) # TRUE
dmgnet1<-diameter(gnet1, directed = TRUE,
                unconnected = F,
                weights = E(gnet1)$distance)
 dmgnet1 #υπολογισμός σταθμισμένης διαμέτρου
  mdnet1<-mean_distance(net1, directed = T, unconnected = F) # μέση απόσταση
 mdnet1
rgnet11<-radius(gnet1, mode="out") # υπολογισμός ακτίνας κατευθυνόμενου out
rgnet11
rgnet12<-radius(gnet1, mode="in") # υπολογισμός ακτίνας κατευθυνόμενου in
rgnet12  
ecgnet11<-eccentricity(net1, vids=V(net1), mode="out")
 ecgnet11 #εκκεντρότητα κορυφών out
ecgnet12<-eccentricity(net1, vids=V(net1), mode="in")
 ecgnet12 #εκκεντρότητα κορυφών in
Cgnet11 <- as_ids(V(net1)[ecgnet11 == min(ecgnet11)])#το κέντρο του δικτύου out
Cgnet12 <- as_ids(V(net1)[ecgnet12 == min(ecgnet12)])#το κέντρο του δικτύου in
V(net1)[Cgnet11]$color<-rainbow(7,alpha=1)[2]
 plot(net1) #αναπαρασταση δικτυου με τονισμενο το κέντρο out
##########################################################################
#
#   4. Μέτρα Κεντρικότητας
#
##########################################################################	
				# ΒΑΘΜΙΚΗ ΚΕΝΤΡΙΚΟΤΗΤΑ
dcnet1<-degree(net1, mode = "all", loops = F,
              normalized = F)			  #ο βαθμός κάθε κορυφής
dcnet1<-degree(net1, mode = "all", loops = F,
              normalized = TRUE)			#κανονικοποιημένος
dcnet1NS<-sort(round(dcnet1N,3))
dcnet1NS # στρογγυλοποίηση σε 3 δεκαδικά και διάταξη σε αύξουσα σειρά
dcnet1OUT<-degree(net1, mode = "out", loops = F,
                 normalized = F)
dcnet1NOUT<-degree(net1, mode = "out", loops = F,
dcnet1NOUTS<-sort(round(dcnet1NOUT,3))
                normalized = TRUE)
dcnet1NOUTS  #κανονικοποιημένος βαθμός κεντρικότητας out, στρογγυλοποιημένος σε 3 δεκαδικα ψηφία 
			  #διάταξη σε αύξουσα σειρα
dcnet1NOUT<-degree(net1, mode = "out", loops = F,
                    normalized = TRUE)
dcnet1IN<-degree(net1, mode = "in", loops = F,
                  normalized = F)
dcnet1NIN<-degree(net1, mode = "in", loops = F,
                   normalized = TRUE)
 dcnet1NINS<-sort(round(dcnet1NIN,3))
 dcnet1NINS		#κανονικοποιημένος βαθμός κεντρικότητας out, στρογγυλοποιημένος σε 3 δεκαδικα ψηφία 
				#διάταξη σε αύξουσα σειρα

min_d <- min(nchar) # χρησιμοποιηθεί για τον καθορισμό
max_d <- max(nchar) # του μεγέθους των κορυφών
rscl <- ((high-low)*(nchar-min_d))/(max_d-min_d)+low
rscl
rescale <- function(nchar,low,high) {# η συνάρτηση αυτή θα
min_d <- min(nchar) # χρησιμοποιηθεί για τον καθορισμό
max_d <- max(nchar) # του μεγέθους των κορυφών
rscl <- ((high-low)*(nchar-min_d))/(max_d-min_d)+low
rscl
sizdc<-rescale(dcnet1N,1,15)#προσαρμογή του μεγέθους των κορυφών στη βαθμική κεντρικότητα
plot(net1 , 
edge.width=1, 
vertex.shape="circle", 
vertex.frame.color=V(net1)$color, 
vertex.size=sizdc, 
vertex.label=V(net1)$name, 
vertex.label.font=2, 
vertex.label.cex=1, 
vertex.label.dist=2, 
vertex.label.degree=-pi/2) #αναπαρασταση δικτύου με προσαρμογή της κεντρικότας στο μέγεθος των κορυφών

			# ΑΛΓΟΡΙΘΜΟΣ Hubs και Authority
hbnet1<-hub_score(net1, scale = TRUE, weights = NULL)
sort(round(hbnet1$vector,3)) #αλγοριθμος Hubs ,  σε αύξουσα σειρά, 3 δεκαδικά
asnet1<-authority_score(net1, scale = TRUE, weights = NULL)
sort(round(asnet1$vector,3)) #αλγοριθμος Authority, σε αύξουσα σειρά, 3 δεκαδικά
sizhb<-rescale(hbnet1$vector,1,15) 
sizhb
sizas<-rescale(asnet1$vector,1,15)
sizas #προσαρμογή μεγέθους κορυφών
plot(net1 ,
vertex.shape="circle", 
vertex.frame.color=V(net1)$color, 
vertex.size=sizhb, 
vertex.label=V(net1)$name, 
vertex.label.font=2, 
vertex.label.cex=1, 
vertex.label.dist=2, 
vertex.label.degree=-pi/2,
main="Hubs") #plot Hubs
plot(net1 , 
vertex.shape="circle", 
vertex.frame.color=V(net1)$color, 
vertex.size=sizas, 
vertex.label=V(net1)$name, 
vertex.label.font=2, 
vertex.label.cex=1,
vertex.label.dist=2, 
vertex.label.degree=-pi/2,
main="Authorities") #plot Authority
			##PageRank
prnet1<-page_rank(net1, algo = "arpack",
directed = TRUE, damping = 0.85,
personalized = NULL, weights = NULL)
prnet1
prnet1N<-prnet1$vector 
prnet1N
sort(round(prnet1N,3))
sizpr<-rescale(prnet1N,1,15) 
sizpr #προσαρμογή μεγέθους κορυφών σε κεντρικοτητα
plot(net1 , 
vertex.shape="circle", 
vertex.frame.color=V(net1)$color, 
vertex.size=sizpr, 
vertex.label=V(net1)$name, 
vertex.label.font=2, 
vertex.label.cex=1, 
vertex.label.dist=2, 
vertex.label.degree=-pi/2,
main="PageRank")
			#ΚΕΝΤΡΙΚΟΤΗΤΑ ΕΓΓΥΤΗΤΑΣ - ΚΑΤΕΥΘΥΝΟΜΕΝΟ
 ccnet11<-round(closeness(net1, mode = "in", weights = E(net1)$distance, normalized = T), 3)
 ccnet11
 sccnet11<-sort(ccnet11)
  sccnet11 
 sizccD<-rescale(ccnet11,1,15)
 sizccD #προσαρμογή μεθέθους κορυφών σε κεντρικοτητα
 plot(net1 , 
      vertex.frame.color=V(net1)$color, 
      vertex.label=V(net1)$name, 
      vertex.label.font=2, 
      vertex.label.cex=1, 
      vertex.label.dist=2, 
      vertex.label.degree=-pi/2,
      main="In-Closeness Centr")
	  
			#ΔΙΑΜΕΣΟΤΗΤΑ#
bcnet11<-round(betweenness(net1,
directed = TRUE,
weights = NULL,
nobigint = TRUE,
normalized = TRUE), 4)
sbcnet11<-sort(bcnet11)
bcnet11
sbcnet11
siznet11D<-rescale(bcnet11,1,15)
siznet11D #προσαρμογή μεγέθους κορυφών σε κεντρικοτητα
plot(net1 , 
vertex.shape="circle", 
vertex.frame.color=V(net1)$color, 
vertex.size=siznet11D, 
vertex.label=V(net1)$name, 
vertex.label.font=2, 
vertex.label.cex=1, 
vertex.label.dist=2, 
vertex.label.degree=-pi/2,
main="In-Betweenness Centr")


