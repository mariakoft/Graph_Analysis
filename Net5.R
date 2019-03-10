library(igraph)
library(Matrix)
library(lattice)
library(vcd)
setwd("C:/Users/User/Dropbox/Networks 2017-18/data-assignment-instructions 2017-18/data-relations(11)/knowledge nets(3)")
# βάλε  το path του φακέλου που έχεις το dataset 
ts <-read.table("Seeks meta-knowledge from.txt", check.names = TRUE)
class(ts)

ms<-as.matrix(ts) #μετατροπή σε πίνακα
class(ms)
row.names(ms)
colnames(ms)
isSymmetric(ms) #FALSE άρα κατευθυνόμενο (σε συμφωνία και με το σχεδιασμό της 
                # της έρευνας)
gS<-graph_from_adjacency_matrix(ms, 
                                mode = "directed",
                                weighted = TRUE,
                                diag = TRUE,
                                add.colnames = NULL,
                                add.rownames = NA)

##### Γενικές Πληροφορίες Δικτυου ######

print_all(gS) # αναλυτική παράθεση στοιχείων
summary(gS) # η σύντομη περιγραφή
gS
is_simple(gS) # TRUE 
gorder(gS) # πλήθος κορυφών (τάξη δικτύου)
gsize(gS) # πλήθος ακμών (μέγεθος δικτύου)


is.directed(gS)   #TRUE άρα το δίκτυο είναι κατευθυνόμενο

is.connected(gS, mode = "weak")  # εξετάζουμε τη συνδετικότητα 
                                 # χωρίς να λάβουμε υπόψη την κατεύθυνση
                                 # TRUE

is.connected(gS, mode = "strong")  # το ίδιο λαμβάνοντας υπόψη την κατεύθυνση
                                   # FALSE

is.weighted(gS) #έχει οριστεί από εμάς attribute: "weight"

is.named(gS) # έχει ονόματα στις κορυφές

is_bipartite(gS)  # αν υπάρχει attribute: "type" για τις κορυφές
girth(as.undirected(gS, edge.attr.comb = "ignore")) #μήκος μικρότερου κύκλου

vcount(gS)                   #υπολογίζουμε τον αριθμό των κορυφών
ecount(gS)                  #υπολογίζουμε τον αριμό των ακμών
ecount(gS)==vcount(gS)-1  # είναι δένδρο? (θεωρία : πλήθος ακμών = πλήθος κορυφών -1)

vertex.connectivity(gS)# υπολογίζουμε τον συνδετικό αριθμό των κορυφών, που σημαίνει ότι 
                       #ο ελάχιστος αριθμός κόμβων που αν διαγραφεί κάνει το δίκτυο μη-συνδετικό είναι ένας κόμβος
edge.connectivity(gS)  # υπολογίζουμε τον συνδετικό αριθμό των ακμών, που σημαίνει ότι αν αφαιρεθεί 1  ακμή,
                       #το δίκτυο γίνεται μη-συνδετικό
graph.density(gS)      #υπολογίζουμε την πυκνότητα του δικτύου: 0.3875, άρα το δίκτυο΄ είναι πυκνο
# Πληροφορίες για τις ακμές και τις κορυφές #
V(gS)
vdf<-as_data_frame(gS, what = "vertices")  #Δημιουργεί data frame με πληροφορίες για τις κορυφές
vdf
E(gS) #Για τις ακμές (με τη μορφή διανύσματος χαρακτήρων)


outD<-degree(gS, mode = "out")  # ακολουθίες των έξω βαθμών
outD
inD<-degree(gS, mode = "in")  # ακολουθίες των έσω βαθμών
inD
sum(outD)==sum(inD)
plot(gS,edge.arrow.size=0.2)


mutgS<-which_mutual(gS)  #εύρεση αμοιβαίων σχέσεων
gSwm<-gSm-edges(E(gS)[mutgS]) #υπογράφημα με αφαίρεση ακμών που δεν ειναι αμοιβαίες σχέσεις

gSm<-simplify(as.undirected(gS, mode = "mutual"), remove.multiple = TRUE)
summary(gSm)
plot(gSm) #αναπαρασταση του δικτυου των αμοιβαιων σχέσεων

cn<-clique_num(gSm) 
cn

vmc<-largest.cliques(gSm)
vmc
class(vmc[[1]])
cq1<-as_ids(vmc[[1]])
cq2<-as_ids(vmc[[2]])
cq3<-as_ids(vmc[[3]])
cq4<-as_ids(vmc[[4]])
cq5<-as_ids(vmc[[5]])
cq6<-as_ids(vmc[[6]])
gSt<-gS-vertices(c(cq1,cq2,cq3,cq4,cq5,cq6))
plot.igraph(gSt, edge.arrow.size=.3)
title("Δίκτυο εκτός των Μέγιστων Κλικών")

gSt2<-gSm-vertices(c(cq1,cq2,cq3,cq4,cq5,cq6))
plot.igraph(gSt2, edge.arrow.size=.3)


par(mfrow=c(1,3))
plot.igraph(gS, edge.arrow.size=.1)
title("Αρχικό Δίκτυο")
plot.igraph(gSm, layout=layout_with_fr, weights="weight")
title("Δίκτυο Αμοιβαίων Σχέσεων")
plot.igraph(gSwm, edge.arrow.size=0.2)
title("Δίκτυο εκτός Αμοιβαίων Σχέσεων")
par(mfrow=c(1,1))






####################################### 4o ##########################################


gSlay <- layout_with_fr(gS)  # συντεταγμένες κορυφών
V(gS)$color<- rainbow(7,alpha=.9)[5]
E(gS)$arrow.size=0.2
E(gS)$color="darkgray"
E(gS)$distance<-1/(0.5+E(gS)$weight)
summary(gS)      # ποια νέα attributes προστέθηκαν


plot(gS , #layout=gSlay,             
     vertex.shape="circle",            
     vertex.frame.color= V(gS)$color,     
     vertex.size=10,                  
     vertex.label=V(gS)$name,         
     vertex.label.font=2,		         
     vertex.label.cex=1,          
     vertex.label.dist=2,       
     vertex.label.degree=-pi/2)		  


# Η απόσταση που χρησιμοποιούμε είναι η αντίστροφη του βάρους
DtgS1<-distances(gS, v = V(gS), mode = "out", weights = E(gS)$distance, 
                 algorithm = "automatic")  # πίνακας αποστάσεων out
DtgS1
DtgS2<-distances(gS, v = V(gS), mode = "in", weights = E(gS)$distance, 
                 algorithm = "automatic")  # πίνακας αποστάσεων in
DtgS2
isSymmetric(DtgS1)
isSymmetric(DtgS2)  # η σχέση μεταξύ των δυο πινάκων (ισχυρά συνδετικό)
DtgS1==t(DtgS2)  #TRUE  
sum(DtgS1==t(DtgS2))==vcount(gS)^2 # Πλήθος(TRUE)=n^2



####### τα σημαντικότερα μέτρα από τη θεωρία γραφημάτων #######

rgS1<-radius(gS, mode="out") # υπολογισμός ακτίνας κατευθυνόμενου (out)
rgS1
rgS2<-radius(gS, mode="in") # υπολογισμός ακτίνας κατευθυνόμενου (in)
rgS2                        # δεν δέχεται βάρη !!!!

#υπολογισμός διαμέτρου και μέσης απόστασης χωρίς βάρη
dmgS2<-diameter(gS, directed = TRUE, 
                unconnected = F, 
                weights = NA)       # υπολογισμός  διαμέτρου
dmgS2
mdgS<-mean_distance(gS, directed = T, unconnected = F) # μέση απόσταση
mdgS


ecgS1<-eccentricity(gS, vids=V(gS), mode="out") # εκκεντρότητα κορυφών (out)
ecgS1                                           # χωρίς στάθμιση
CgS1 <- as_ids(V(gS)[ecgS1 == min(ecgS1)]) # το κέντρο του δικτύου (out)
CgS1 

ecgS2<-eccentricity(gS, vids=V(gS), mode="in") # εκκεντρότητα κορυφών (in)
ecgS2                                           # χωρίς στάθμιση
CgS2 <- as_ids(V(gS)[ecgS2 == min(ecgS2)]) # το κέντρο του δικτύου (in)
CgS2 

V(gS)[CgS1]$color<-rainbow(7,alpha=1)[2] # διαφορετικό χρώμα στο κέντρο (out)
plot(gS , #layout=gSlay,             
     vertex.shape="circle",            
     vertex.frame.color= V(gS)$color,     
     vertex.size=10,                  
     vertex.label=V(gS)$name,         
     vertex.label.font=2,		         
     vertex.label.cex=1,          
     vertex.label.dist=2,       
     vertex.label.degree=-pi/2)	
V(gS)[CgS1]$color<-rainbow(7,alpha=.9)[5]

V(gS)[CgS2]$color<-rainbow(7,alpha=1)[2] # διαφορετικό χρώμα στο κέντρο (in)
plot(gS , #layout=gSlay,             
     vertex.shape="circle",            
     vertex.frame.color= V(gS)$color,     
     vertex.size=10,                  
     vertex.label=V(gS)$name,         
     vertex.label.font=2,		         
     vertex.label.cex=1,          
     vertex.label.dist=2,       
     vertex.label.degree=-pi/2)	
V(gS)[CgS2]$color<-rainbow(7,alpha=.9)[5]



###### ελάχιστο δεντρο ζεύξης ######

gsMST<-mst(gS, weights = E(gS)$distance, algorithm="prim")
gsMST

par(mfrow=c(1,2))
plot(gS , #layout=gFlay,             
     vertex.shape="circle",            
     vertex.frame.color= V(gS)$color,     
     vertex.size=10,                  
     vertex.label=V(gS)$name,         
     vertex.label.font=2,		         
     vertex.label.cex=1,          
     vertex.label.dist=2,       
     vertex.label.degree=-pi/2)	
title("Το αρχικό δίκτυο")

plot(gsMST,  #layout=gFlay,             
     vertex.shape="circle",            
     vertex.frame.color= V(gS)$color,     
     vertex.size=10,                  
     vertex.label=V(gS)$name,         
     vertex.label.font=2,		         
     vertex.label.cex=1,          
     vertex.label.dist=2,       
     vertex.label.degree=-pi/2)
title("Ελάχιστο δένδρο ζεύξης")

par(mfrow=c(1,1))




###### κεντρικότητες ######


###### βαθμική in, out και συνολικά ####
dcgS<-degree(gS, mode = "all", loops = F,
              normalized = F)
dcgS

dcgSN<-degree(gS, mode = "all", loops = F,
              normalized = TRUE)
dcgSN

dcgSout<-degree(gS, mode = "out", loops = F,
             normalized = F)
dcgSout
dcgSNout<-degree(gS, mode = "out", loops = F,
              normalized = TRUE)
dcgSNout
dcgSin<-degree(gS, mode = "in", loops = F,
             normalized = F)
dcgSin
dcgSNin<-degree(gS, mode = "in", loops = F,
              normalized = TRUE)
dcgSNin

dcgSNS<-sort(round(dcgSN,3)) 
dcgSNS
dcgSNSout<-sort(round(dcgSNout,3)) 
dcgSNSout
dcgSNSin<-sort(round(dcgSNin,3)) 
dcgSNSin

par(mfrow=c(1,1))
rescale <- function(nchar,low,high) {  # προσαρμογή του μεγέθους των 
   min_d <- min(nchar)                  # κορυφών στη βαθμική κεντρικότητα
   max_d <- max(nchar)
   rscl <- ((high-low)*(nchar-min_d))/(max_d-min_d)+low  
   rscl
}
sizdc<-rescale(dcgSN,1,15)
sizdc


V(gS)[sizdc]$color<-rainbow(7,alpha=.9)[5]

sizout<-rescale(outD,1,15)
sizin<-rescale(inD,1,15)

####γραφημα για τον βαθμό out###
plot(gS , #layout=gSlay,             # συντεταγμένες
      edge.width=0.3,                 # πάχος ακμών
      edge.arrow.size=0.2,
      vertex.shape="circle",           # σχήμα κορυφών 
      vertex.frame.color=V(gS)$color,  # χρώμα περιφέρειας κορυφών
      vertex.size=sizout,                  # μέγεθος κορυφών
      vertex.label=V(gS)$name,         # ετικέτες κορυφών
      vertex.label.font=2,		         # γραμματοσειρά ετικετών
      vertex.label.cex=1,          # μέγεθος γραμματοσειράς
      vertex.label.dist=2,       # απόσταση ετικέτας από την κορυφή
      vertex.label.degree=-pi/2)

####γραφημα για τον βαθμό in###
plot(gS , #layout=gSlay,             # συντεταγμένες
      edge.width=0.3,                 # πάχος ακμών
      edge.arrow.size=0.2,
      vertex.shape="circle",           # σχήμα κορυφών 
      vertex.frame.color=V(gS)$color,  # χρώμα περιφέρειας κορυφών
      vertex.size=sizin,                  # μέγεθος κορυφών
      vertex.label=V(gS)$name,         # ετικέτες κορυφών
      vertex.label.font=2,		         # γραμματοσειρά ετικετών
      vertex.label.cex=1,          # μέγεθος γραμματοσειράς
      vertex.label.dist=2,       # απόσταση ετικέτας από την κορυφή
      vertex.label.degree=-pi/2)


####γραφημα για τον βαθμό in και out μαζι###
plot(gS , #layout=gSlay,             # συντεταγμένες
                  edge.width=0.3,                 # πάχος ακμών
                 edge.arrow.size=0.2,
                  vertex.shape="circle",           # σχήμα κορυφών 
                  vertex.frame.color=V(gS)$color,  # χρώμα περιφέρειας κορυφών
                  vertex.size=sizdc,                  # μέγεθος κορυφών
                  vertex.label=V(gS)$name,         # ετικέτες κορυφών
                  vertex.label.font=2,		         # γραμματοσειρά ετικετών
                  vertex.label.cex=1,          # μέγεθος γραμματοσειράς
                  vertex.label.dist=2,       # απόσταση ετικέτας από την κορυφή
                  vertex.label.degree=-pi/2)


######## hubs & authorities ########
hbgS<-hub_score(gS, scale = TRUE, weights = NULL)
sort(round(hbgS$vector,3))

asgS<-authority_score(gS, scale = TRUE, weights = NULL)
sort(round(asgS$vector,3))
sizhb<-rescale(hbgS$vector,1,15) # προσαρμογή μεγέθους κορυφών
sizhb
sizas<-rescale(asgS$vector,1,15)
sizas

par(mfrow=c(1,2))
V(gS)[CgS1]$color<-rainbow(7,alpha=.9)[5]
plot(gS , #layout=gSlay,             # συντεταγμένες
     edge.arrow.size=0.02,
      vertex.shape="circle",           # σχήμα κορυφών 
      vertex.frame.color=V(gS)$color,  # χρώμα περιφέρειας κορυφών
     vertex.size=sizhb,                  # μέγεθος κορυφών
      vertex.label=V(gS)$name,         # ετικέτες κορυφών
      vertex.label.font=2,		         # γραμματοσειρά ετικετών
     vertex.label.cex=1,          # μέγεθος γραμματοσειράς
      vertex.label.dist=2,       # απόσταση ετικέτας από την κορυφή
      vertex.label.degree=-pi/2,
      main="Hubs")
 
plot(gS , #layout=gSlay,             # συντεταγμένες
     edge.arrow.size=0.02,
      vertex.shape="circle",           # σχήμα κορυφών 
      vertex.frame.color=V(gS)$color,  # χρώμα περιφέρειας κορυφών
      vertex.size=sizas,                  # μέγεθος κορυφών
      vertex.label=V(gS)$name,         # ετικέτες κορυφών
      vertex.label.font=2,		         # γραμματοσειρά ετικετών
      vertex.label.cex=1,          # μέγεθος γραμματοσειράς
      vertex.label.dist=2,       # απόσταση ετικέτας από την κορυφή
      vertex.label.degree=-pi/2,
      main="Authorities")
par(mfrow=c(1,1))



########## page rank #########
prgS<-page_rank(gS, algo = "arpack",
                 directed = TRUE, damping = 0.85, 
                 personalized = NULL, weights = NULL)
prgS
prggN<-prgS$vector
prggN
sort(round(prggN,3))
sizpr<-rescale(prggN,1,15)
sizpr
par(mfrow=c(1,1))
plot(gS , #layout=gSlay,             # συντεταγμένες
     edge.arrow.size=0.2,
     vertex.shape="circle",           # σχήμα κορυφών 
     vertex.frame.color=V(gS)$color,  # χρώμα περιφέρειας κορυφών
     vertex.size=sizpr,                  # μέγεθος κορυφών
     vertex.label=V(gS)$name,         # ετικέτες κορυφών
     vertex.label.font=2,		         # γραμματοσειρά ετικετών
     vertex.label.cex=1,          # μέγεθος γραμματοσειράς
     vertex.label.dist=2,       # απόσταση ετικέτας από την κορυφή
     vertex.label.degree=-pi/2,
     main="PageRank")



######### closeness ########
ccgSin<-round(closeness(gS, mode = "in",
                      weights = E(gS)$distance, normalized = F), 3)
ccgSin
################################################################# 
sccgSin<-sort(ccgSin)
sccgSin
sizccDin<-rescale(ccgSin,1,15)
sizccDin
par(mfrow=c(1,1))
plot(gS , #layout=gSlay,             # συντεταγμένες
     edge.arrow.size=0.2,
     vertex.shape="circle",           # σχήμα κορυφών 
     vertex.frame.color=V(gS)$color,  # χρώμα περιφέρειας κορυφών
     vertex.size=sizccDin,                  # μέγεθος κορυφών
     vertex.label=V(gS)$name,         # ετικέτες κορυφών
     vertex.label.font=2,		         # γραμματοσειρά ετικετών
     vertex.label.cex=1,          # μέγεθος γραμματοσειράς
     vertex.label.dist=2,       # απόσταση ετικέτας από την κορυφή
     vertex.label.degree=-pi/2,
     main="In-Closeness Centr")
par(mfrow=c(1,1))


ccgSout<-round(closeness(gS, mode = "out",
                      weights = E(gS)$distance, normalized = F), 3)
ccgSout
################################################################# 
sccgSout<-sort(ccgSout)
sccgSout
sizccDout<-rescale(ccgSout,1,15)
sizccDout
par(mfrow=c(1,1))
plot(gS , #layout=gSlay,             # συντεταγμένες
     edge.arrow.size=0.2,
     vertex.shape="circle",           # σχήμα κορυφών 
     vertex.frame.color=V(gS)$color,  # χρώμα περιφέρειας κορυφών
     vertex.size=sizccDout,                  # μέγεθος κορυφών
     vertex.label=V(gS)$name,         # ετικέτες κορυφών
     vertex.label.font=2,		         # γραμματοσειρά ετικετών
     vertex.label.cex=1,          # μέγεθος γραμματοσειράς
     vertex.label.dist=2,       # απόσταση ετικέτας από την κορυφή
     vertex.label.degree=-pi/2,
     main="Out-Closeness Centr")
par(mfrow=c(1,1))


######### betweenness ########

bcgSin<-round(betweenness(gS, 
                        directed = TRUE, 
                        weights = NULL,
                        nobigint = TRUE, 
                        normalized = TRUE), 4)
sbcgSin<-sort(bcgSin)
bcgSin
sbcgSin

sizbcDin<-rescale(bcgSin,1,15)
sizbcDin
par(mfrow=c(1,1))
plot(gS , #layout=gSlay,             # συντεταγμένες
     vertex.shape="circle",           # σχήμα κορυφών 
     vertex.frame.color=V(gS)$color,  # χρώμα περιφέρειας κορυφών
     vertex.size=sizbcDin,                  # μέγεθος κορυφών
     vertex.label=V(gS)$name,         # ετικέτες κορυφών
     vertex.label.font=2,		         # γραμματοσειρά ετικετών
     vertex.label.cex=1,          # μέγεθος γραμματοσειράς
     vertex.label.dist=2,       # απόσταση ετικέτας από την κορυφή
     vertex.label.degree=-pi/2,
     main="Betweenness Centr")
par(mfrow=c(1,1))


