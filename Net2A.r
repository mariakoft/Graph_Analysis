# gia ellhnika
# Sys.setlocale("LC_CTYPE","Greek")
# Δίκτυο 2, Φουδούλη Αναστασία
setwd("C:/Users/xx/Desktop/δικτυο2.1/diktyo2")

#194,

# CA-GrQc, δίκτυο συνεργασίας
# επιστήμενες στον τομέα γενικής σχετικότητας και κοσμολογίας

# βιβλιοθήκες που θα χρησιμοποιηθούν
library(igraph)
library(Matrix)
library(lattice)
library(vcd)

#############################################################################
####  1. Τάξη, μέγεθος, συνδετικότητα, συνιστώσες, βαθμός και εκτιμήσεις ####
####     κατανομών      #####################################################
#############################################################################

#κατασκευή δικτύου απο αρχείο .txt 
listg <-read.table("CA-GrQc.txt", header = TRUE)
g0<-graph_from_data_frame(listg, directed = FALSE)
summary(g0)   #ένα μη κατευθυνόμενο γράφημα 5242 κορυφών και 28980 ακμών

is.connected(g0)   #όχι, δεν είναι συνδετικό
is.named(g0)       #ναι, οι κορυφές έχουν ονόματα
is.directed(g0)    #όχι, άλλωστε ζήτησα και directed=FALSE
is.simple(g0)      #όχι,άρα έχει πολλάπλές ακμές ή,και λούπες
is.weighted(g0)    #όχι, κάθε ακμή δηλώνει μία συνεργασία

#Κάνω αρχικά το δίκτυο απλό 
g<-simplify(g0)
is.simple(g)    #TRUE
vcount(g)                         #5242 κορυφές, είναι η τάξη του δικτύου, εναλλακτικά gorder(g)
ecount(g); ecount(g)/ecount(g0)   #14484 παρέμειναν μετά την εντολή simplify, το 49,9793% του αρχικού δικτύου
                                  # είναι το μέγεθος του δικτύου, ενναλακτικά gsize(g)
is.connected(g)   # FALSE
graph.density(g)  # αραιό 0,001054
is.weighted(g)    # FALSE

#για να προχωρήσω σε μία γραφική αναπαράσταση του δικτύου πρώτα βρίσκω συντεταγμένες 
# layg<-layout_with_graphopt(g)
# save(layg, file="CoordsforGrQc")
load("CoordsforGrQc")
V(g)$coord1 <-layg[ ,1]
V(g)$coord2 <-layg[ ,2]
par(bg="black", new=FALSE)
plot(g, layout=layg, edge.color="gray50", edge.width=0.01, vertex.shape="circle",
     vertex.frame.color="firebrick4",vertex.color="firebrick4", vertex.size=0.6, 
     vertex.label=NA)

#zoom in στο plot
plot(g, layout=layg, edge.color="gray50", edge.width=0.01, vertex.shape="circle",
     vertex.frame.color="firebrick4",vertex.color="firebrick4", vertex.size=0.6, 
     vertex.label=NA, xlim=c(-0.4, 0.4), ylim=c(0.3,0.9))

#Το δίκτυό μου επιπλέον δεν έιναι συνδετικό, άρα θα βρω τις συνιστώσες του 
c<-clusters(g)
c$no                              # 355 συνιστώσες
c$csize                           # το μέγεθος κάθε μίας από τις 355 συνιστώσες
tail(sort(c$csize))               # οι 6 μεγαλύτερες συνιστώσες (σε αύξουσα σειρά)
c$membership                      # και σε ποιά συνιστώσα ανήκει η κάθε κορυφή του δικτύου

#και θέλω αρχικά να τις δω γραφικά
g1<-g   
cmemb<-c$membership
colorbar<-rainbow(c$no)
par(bg="black",new=FALSE)                         
V(g1)$color <- colorbar[cmemb]
plot(g1, layout=layg, edge.color="lavender", edge.width=0.01, vertex.shape="circle",
     vertex.frame.color=V(g1)$color, vertex.color=V(g1)$color,  vertex.size=1, vertex.label=NA)
	 
	 
#με ενδιαφέρουν οι μεγαλύτερες σε μέγεθος συνιστώσες                               
max(c$csize); max(c$csize)/vcount(g)          #η γιγάντια συνιστώσα περιλαμβάνει 4158 κορυφές του δικτύου
#που αποτελόυν το 79,93% των σνολικών κορυφών του δικτύου
rev(sort(c$csize))[2]; rev(sort(c$csize))[2]/vcount(g)    #η δεύτερη σε μέγεθος συνιστώσα περιπλαμβάνει 14
#κορυφές, μόλις το 0,26% του δικτύου
sort(c$csize)  #βλέπω το μέγεθος κάθε συνιστώσας σε αύξουσα σειρά
summary(c$csize)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    #πολύ λοξή κατανομή
#1.00    2.00    2.00   14.77    3.00 4158.00 

#η κατανομή μοιάζει να είναι Power law, θα το δω γραφικά και θα το ελέγξω
cd<-cluster.distribution(g, cumulative = FALSE)   #frequency distribution of components
ccd<-cluster.distribution(g, cumulative = TRUE)   #cumulative frequency distribution of components

#γραφικά λοιπόν, σε διάγραμματα κατανομής                   
par(mfrow = c(1, 3))
par(bg="white", new=FALSE)
plot(cd, xlab="size of component", ylab="freq")
plot(cd, xlab="Log size of component", ylab="Log freq", log="xy")
plot(ccd, xlab="Log size of component", ylab="Log CCD", log="xy")
par(mfrow=c(1,1))

#μοιάζει να είναι δυναμοκατανομή , και γι αυτό προχωράω στον έλεγχο (cx^(-a))
fitpl1<-power.law.fit(c$csize, xmin=NULL, start=2, force.continuous=FALSE,
                      implementation= "plfit" )   
fitpl1$alpha   # 3.262853
fitpl1$xmin    # 3
fitpl1$KS.p    # 0.9947234>0.1 άρα γίνεται δεκτός ο έλεγχος ότι το μέγεθος των συνιστωσών
# ακολουθεί δυναμοκατανομή
# χρησιμοποιώ λοιπόν ότι βρήκα στη fitpl2                                          
fitpl2 <-power.law.fit(c$csize, xmin=3, start=2, force.continuous=FALSE,
                       implementation= "R.mle")  
summary(fitpl2)       #alpha=3.262864,  std.error=0.1750832

pl1<-function(x) x^(-3.262864)       # component distribution
pl2<-function(x) x^(-3.262864+1)     # ccd (comulative component distribution)


#και θα το δω και γραφικά προσθέτοντας τις συναρτήσεις εκτίμησης
par(mfrow = c(1, 3))
D<-max(c$csize);D  #4158
plot(cd, xlab="size of component", ylab="freq");
plot(pl1, 3, D,n=2000,add=TRUE,type="l")
plot(cd, xlab="Log size of component", ylab="Log freq", log="xy");
plot(pl1, 1.2, D-1000,n=2000,add=TRUE,type="l")
plot(ccd, xlab="Log size of component", ylab="Log CCD", log="xy");
plot(pl2, 3, D,n=2000,add=TRUE,type="l")
par(mfrow = c(1, 1))
#φαίνεται λοιπόν και γραφικά ότι έχει γίνει μία καλή εκτίμηση


######################################################################################
# θα κοιτάξω στη συνέχεια την κατανομή των βαθμών του δικτύου g
summary(degree(g))      #ξανά λοξή κατανομή
mean(degree(g))         #5.526135
sd(degree(g))           #7.918456
sd(degree(g))/mean(degree(g))  #143.291% συντελεστής μεταβλητότητας

#η κατανομή που ακολουθείται
de <- as.numeric(degree(g))
D<-max(de)     #81, ο μεγαλύτερος βαθμός
dd<-degree.distribution(g, cumulative = FALSE)      #degree distribution
cdd<-degree.distribution(g, cumulative = TRUE)      #comulative degree distribution

#και γραφικά
par(bg="white",new=FALSE) 
par(mfrow = c(1, 3))
plot(dd, xlab="degree", ylab="freq")
plot(dd, xlab="Log degree", ylab="Log freq", log="xy")
plot(cdd, xlab="Log degree", ylab="Log CCD", log="xy")
par(mfrow=c(1,1))

#θυμίζει ελαφρώς power law αλλά μάλλον δεν είναι
#παρόλα αυτά προχωράω σε έλεγχο
fitpl3<-power.law.fit(de1, xmin=NULL, start=2, force.continuous=FALSE,
                      implementation= "plfit" )    
fitpl3$alpha     #2.112951
fitpl3$xmin      #3
fitpl3$KS.p      # 2.456303e-05
# το p που προκείπτει έιναι πολύ μοκρο <0,05, άρα η υπόθεση ότι οι κατανομή των βαθμών είναι power law
# απορίπτεται (σύμφωνα με το help)

#δεύτερη δοκιμή για κατανομή poisson
gf<-goodfit(de, type="poisson",method="ML") 
summary(gf)  #απορίπτεται επίσης

#τρίτη δοκιμή για lognormal
ks<-ks.test(log(de), "pnorm", mean(log(de)), sd(log(de)))
summary(ks)
# ties should not be present for the Kolmogorov-Smirnov test
# το τεστ δεν εφαρμόζεται γιατί έχω εμφάνιση των ίδιων τιμών

par(mfrow=c(1,2))
plot(dd, xlab="degree", ylab="freq", main="Degree Distribution")
hist(de, main="Histogram of degrees",xlab="Degree")
par(mfrow=c(1,1))
# με τους παραπάνω ελέγχους δεν είναι δυνατόν να γίνει μία καλή εκτίμηση της κατανομής 
# που ακολουθούν οι βαθμοί

########################################################################
# δεν περνάω σε ελέγχους με βάρη αφού
is.weighted(g)    #FALSE

#############################################################
###  2.Αποστάσεις (διάμετρος, μέση απόσταση)       ##########
#############################################################

c$membership 
#όλες οι κορυφές που φαίνονται ανήκουν στη πρώτη συνιστώσα, άρα είναι η γιγαντια συνιστώσα
gc<-induced.subgraph(g,V(g)[c$membership == 1])  
gc

#Θέλω να δω μία πρώτη μορφή της
# laygc<-layout.graphopt(gc)
# save(laygc, file="CoordsforGCGrQc")

load("CoordsforGCGrQc")
V(gc)$coord1 <-laygc[ ,1]
V(gc)$coord2 <-laygc[ ,2]
par(bg="black",new=FALSE)
plot(gc, layout=laygc, ede.color="grey67", edge.width=0.02, vertex.shape="circle",
     vertex.frame.color="firebrick4", vertex.color="firebrick4", vertex.size=1, 
     vertex.label=NA)

#ελέγχω συνδετικό και γραμμοσυνδετικό αριθμό
vertex.connectivity(gc)  #1
edge.connectivity(gc)    #1

#διάμετρος δικτύου και gc
#Διάμετρο ενός δικτύου ανομάζουμαι την μεγαλύετη από τις αποστάσεις του 
#Πρακτικά πρόκειται για την κοντινότερη απόσταση μεταξύ των δύο πιο απομακρισμένων κορυφών του δικτύου
dg<-diameter(g, unconnected = TRUE, weights = NA) ; dg     #17, για το δίκτυο
dgc<-diameter(gc, unconnected = FALSE, weights = NA) ; dgc    #17 και για τη gc

#μέση απόσταση
aplg<-average.path.length(g, directed=FALSE, unconnected=TRUE); aplg   #6.048515 για το δίκτυο
aplgc<-average.path.length(gc,directed = FALSE, unconnected = FALSE); aplgc   #6.04938 για τη gc

#για να δουμε και γραφικά την κατανομή των αποστάσεων
#distance_table calculates a histogram, by calculating the shortest path length between each pair of vertices.
dt<-distance_table (gc, directed = FALSE) 

## sum(dt$res);vcount(gc)*(vcount(gc)-1)/2 these are all distances

##      1<=DISTANCE<=17 : διάμετρος
par(bg="white")
plot(1:17, dt$res, xlab="distance", ylab="paths",frame=TRUE)

#η μέση απόσταση ακολουθεί κανονικη κατανομη ~N(mean(plh$res),sd(plh$res))
ks.test(dt$res,"pnorm",mean(dt$res),sd(dt$res))
#	 One-sample Kolmogorov-Smirnov test
# data:  plh$res
# D = 0.29159, p-value = 0.08949

####################################################################
### 3. Συντελεστές σύμπλεξης, Μεταβατικότητα, Ιεραρχική δομή    ####
####################################################################

# Συντελεστής σύμπλεξης (clustering coefficient) μια κορυφής i είναι ένα μέτρο που υπολογίζει κατά πόσον
# οι γείτονες αυτής της κορυφής απέχουν από το να αποτελούν "κλίκα".
# Ορίζεται βάσει της σχέσης:
# c(i)=τ(i)/ρ(i)
# όπου τ(i)= το πλήθος των τριγώνων που μία κορυφή τους είναι το i
#  και ρ(i)= το πλήθος των τριάδων με μία κορυφή το i
# πεδίο τιμών του c(i) είναι το [0,1]
# 0 είναι "αστέρι" ==> τότε οι γείτονες της κορυφής i ασηματίζουν ένα πλήρως ασυνδετικό δίκτυο
# 1 όταν το δίκτυο είναι κλίκα ==> οι γείτονες του i σχηματίζουν ένα πλήρες γράφημα

# Συντελεστής σύμπλεξης του γραφήματος είναι ο μέσος όρος των συντελεστών σύμπλεξης όλων των κορυφών
ltr<-transitivity(g,type="local",isolates="zero") #ο συντελεστής σύμπλεξης κάθε κορυφής 
ccg<-sum(ltr)/vcount(g)
ccg          # 0.5296358 ο συτελεστής σύμπλεξης

# ο συντελεστής μεταβατικότητας θα μας δείξει σε ποιό ποσοστό αν ο ν ερευνητής έχει δημοσιεύσει με
# τους ερευνητές κ και μ, οι δύο τελευταίοι έχουν συνεργαστεί και μεταξύ τους
# (αν έχουν συνεργαστεί τότε αποτελούν οι τρεις μαζί μια μεταβατική τριάδα)
# Ο συντελεστής μεταβατικότητας για ένα δίκτυο υπολογίζεται από τον τύπο
# Τ(G)=πλήθος μεταβατικών τριάδων/πλήθος δυνατών συνδετικών τριάδων
#     =3xπλήθος τριγώνων του G/πλήθος δυνατών συνδετικών τριάδων του G
transitivity(g,type="global", isolates="zero")   #0.6288945 η μεταβατικότητα των κορυφών 

# θέλω να ελέγξω αν υπάρχει μία "ιεραρχική δομή" στο δίκτυο της γ.σ. 
# Σε κάποιες περιπτώσεις έχει παρατηρήθεί ότι σε λογαριθμικούς άξονες
# όταν αυξάνεται ο βαθμός μειώνεται ο συντελεστής σύμπλεξης, κάτι 
# που κάποιοι ερευνητές ερμηνεύουν ως ύπαρξη ιεραρχικής σχέσης 
# μεταξύ των κορυφών του δικτύου (ειδικότερα όταν η κλίση είναι σχεδόν -1).

# το περνάω σαν πίνακα
dltr<-matrix(c(degree(g), ltr), nrow=vcount(g), ncol=2, byrow=FALSE)
min(dltr)  #  είναι 0 και δενν το θέλω, διώχνω μηδενικές τιμές
dltr<-dltr[dltr[,2]!= 0, ] ;  min(dltr) #0.008333333

# κατασκευάζω διάγραμμα διασποράς με λογαριθμικούς άξονες
plot(dltr[,1], dltr[,2], xlab="Log degree", ylab="Log Clustering Coefficient",
       log="xy", frame=TRUE) 

# και θα κάνω γραμμική παλινδρόμηση
linltr<-lm(log(dltr[,2]) ~ log(dltr[,1]) +1) 
summary(linltr)

# Residuals:
#  Min      1Q  Median      3Q     Max 
#-3.7687 -0.3818  0.1530  0.3218  1.4444 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     0.13564    0.02162   6.274   3.9e-10 ***
# log(dltr[, 1]) -0.41639   0.01234  -33.739  < 2e-16 ***
#  ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.6183 on 3853 degrees of freedom
#Multiple R-squared:  0.2281,	Adjusted R-squared:  0.2279 
#F-statistic:  1138 on 1 and 3853 DF,  p-value: < 2.2e-16

#το p-value είναι πάρα πολύ μικρό και δεν θα έπρεπε να δεχτώ προφανώς μία γραμμική σχέση
#όμως θα την τυπώσω πάνω στο διάγραμμα για να φανεί 
#η κλίση είναι μεγαλύτερη του -1, οπότε δεν μπορώ να υποθέσω οτι υπάρχει ιεραρχική δομή

#γραφικά λοιπόν				
plot(dltr[,1], dltr[,2], xlab="Log degree", ylab="Log Clustering Coefficient",
      log="xy", frame=TRUE) ;
lines(dltr[,1], exp( 0.13564 -0.41639*log(dltr[,1])))  


#####################################################
####  4. Ιδιότητα μικροκόσμου       #################
####  Κοφτερού Μαρία, Φουδούλη Αναστασία   ##########
#####################################################

# 10ς τρόπος
# κατασκευάζω αρχικά τυχαίο δίκτυο με ίδιο πλήθος ακμων και κορυφών
# Το τυχαίο δίκτυο θεωρήται πάντα ότι έχει την ιδίοτητα του μικροκόσμου
set.seed(1)  #για να έχω το ίδιο δίκτυο
erG<-erdos.renyi.game(vcount(g), ecount(g), type= "gnm",
                      directed = FALSE, loops = FALSE)

is.simple(erG) 
is.connected(erG)    # FALSE
graph.density(erG)   #0.001054405
graph.density(erG)==graph.density(g)   #TRUE

# μία ένδειξη για την ιδιότητα μικροκόσμου είναι η απόσταση στο δίκτυό μου
# καθώς και η μέση απόσταση να είναι μικρότερες από τις αντίστοιχες για το τυχαίο δίκτυο
derG<-diameter(erG, unconnected=TRUE, weights=NA); derG  #11
derG>dg   #FALSE


aplerG<-average.path.length(erG, unconnected=TRUE, directed=FALSE); aplerG #5.187879
aplerG>aplg   #FALSE
# άρα με βάση τα παραπάνω , το αρχικό δίκτυο δεν αποτελεί μικρόκοσμο

# ένα ακόμα μέτρο είναι ο συντελεστής watts & strogatz
wsg<-transitivity(g, type="average", isolates= "zero");wsg    #0.5296358
wserG<-transitivity(erG, type="average", isolates= "zero");wserG #0.001080615
# και πάλι η τιμή του πραγματικού δικτύου είναι μεγαλύτερη άρα ακόμη μία ένδειξη
# ότι το δίκτυο δεν είναι μικρόκοσμος

#2ος τρόπος
# βρίσκοντας τις κανονικοποιημένες τιμές 
### αργαει αρκετα, ειδικά για τη διαμετρο, μπορει ταχύτερα να γίνει σε 100 δίκτυα
# μπορει για ταχύτητα να γίνει σε 100
erL<-list()        #φτιάχνω κενή λίστα
for (i in 1:1000) {
  erL[[i]]<-diameter(erdos.renyi.game(5242,14496,type="gnm",directed=FALSE,loop=FALSE),  #βρίσκω διάμετρο σε 1000
                     directed=FALSE,unconnected=TRUE,weights=NA)                         #τυχαία δίκτυα
}
erV<-unlist(erL,recursive=TRUE)    #κάνω τη λιστα διάνυσμα
summary(erV)
mean(erV)    #10.19

erL2<-list()
for (i in 1:1000){
  erL2[[i]]<-transitivity(erdos.renyi.game(5242,14496,type="gnm",directed=FALSE,loop=FALSE),
                          type="average", isolates="zero")
}
erV2<-unlist(erL2,recursive=TRUE)
summary(erV2)
mean(erV2)  #0.001042283

l1<-dg/mean(erV); l1       #1.668302
c1<-wsg/mean(erV2); c1     #490.1245
#Και οι δύο τιμές θα έπρενα να είναι ~1 για να μπορώ να υποθέση ιδιότητα μικροκόσμου. 
#Αφού λοιπόν δεν είναι το δίκτυο μου δεν έχει αυτή την ιδιότητα
#Η κανονικοποίηση γίνεται με διαίρεση της διαμέτρου με τη μέση τιμή
#της διαμέτρου πολλών τυχαίων δυκτίων που δημιουργούσα
#ομοίως και για τον συντελεστή


########################################################################
### 5. Συντελεστές ομοιότητας μεταξύ των κορυφών (ως προς το βαθμό)  ###
########################################################################

# Μέσω του συντελεστή ομοιότητας (ως προς το βαθμό) υπολογίζουμε το ποσοστό των κορυφών μεγάλου
# βαθμού που συνδέονται με κορυφές όμοιου βαθμού
# θα δούμε αν το δίκτυο είναι assortative ή disassortative γραφικά
assortativity.degree(g)   #  0.6593246, assortative
# άρα πράγματι κορυφές μεγάλου βαθμού συνδέονται με όμοιες κορυφές

# υπολογισμός του ANND(average nearest neighbor degree)
# ο αλγοριθμος υπολογίζει το μέσο βαθμό των γειτόνων κάθε κορυφής και το μέσο βαθμό των
# γειτόνων κορυφών ίσου βαθμού
knng<-graph.knn(g,weights = NA)
knng$knnk       #A numeric vector, its length is the maximum (total) vertex degree in the graph. 
                #The first element is the average nearest neighbor degree of vertices with degree one, etc.
                # πάρα πολλές τιμές NaN, αφού κορυφλες βαθμού που δεν εμφανίζεται στο δίκτυο παράγουν
                # τιμή NaN.
knnkg<-na.omit(knng$knnk)
degreesg0<-as.numeric(names(table(degree(g))))
degreesg<-degreesg0[degreesg0!=0]

#γραμμική παλινδρόμησηση
fitknn<-lm(knnkg~degreesg+1)
summary(fitknn)
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-18.4528  -4.7844  -0.5214   5.7146  19.5956 

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  6.62890    1.97162   3.362  0.00132 ** 
#degreesg     0.57877    0.04838  11.963  < 2e-16 ***
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 8.302 on 63 degrees of freedom
#Multiple R-squared:  0.6944,	Adjusted R-squared:  0.6895 
#F-statistic: 143.1 on 1 and 63 DF,  p-value: < 2.2e-16

#γενικά εξαιτίας της μικρή τιμής p (λόγω της μεγάλης διασποράς των τιμών) 
#δεν μπορω να υποθέσω γραμμική συσχέτηση
#παρόλα αυτά, θα το τυπώσω γραφικά για να φανεί οτι το δίκτυο είναι assortative
par(bg='white')
plot(degreesg,knnkg, xlab="degree", ylab="knnkg",frame=TRUE); 
lines(degreesg, 6.6289 +0.57877*degreesg )  

    

#######################################################################################
### 6. Το φαινόμενο Rich Club για το 5% και το 2,5% των κορυφών μεγαλύτερου βαθμού  ####
########################################################################################

# ουσιαστικά ψάχνω να βρω αν υπάρχουν ακμές στο δίκτυο μου μεταξύ κορυφών μεγάλου βαθμού
# και έτσι δημιουργήται η εικόνα "Rich Club"
# ελέγχω στη γιγάντια σινιστώσα


# Πρέπει αρχικά να βρώ κορυφές ποιου βαθμού και πάνω αποτελούν το 5% των πιο υψηλοβαθμων κορυφών
deg<-degree(g)
head(sort(deg)); tail(sort(deg))

vcount(g)*0.05    # 262
# άρα θελώ να πάρω κορυφές βαθμού μεγαλύτερου ή ίσου του βαθμού της 262ης κορυφής
# αν τοποθετήσω τους βαθμούς των κορυφών σε φθίνουσα σειρά
rev(sort(deg))[262]   # 20

#και για να φανεί γραφικά
par(bg="black")
par(mfrow=c(1,2))
rcgc20<-induced.subgraph(gc,V(gc)[degree(gc)>=20])
rcgc20   #283 κορυφές    3456 ακμές
lay20<-layout_with_graphopt(rcgc20)
plot(rcgc20, layout=lay20, edge.color="lightcyan4", edge.width=0.05, 
     vertex.shape="circle", vertex.frame.color="brown4", vertex.color="brown4", 
     vertex.size=(degree(gc)/max(degree(gc)))*4+0.8, vertex.label=NA)
title(main=list("Degree>=20", col="white"))

rcgc20<-induced.subgraph(gc,V(gc)[degree(gc)>20])
rcgc20   #255 κορυφές   3216 ακμές
lay20<-layout_with_graphopt(rcgc20)
plot(rcgc20, layout=lay20, edge.color="lightcyan4", edge.width=0.05, 
     vertex.shape="circle", vertex.frame.color="brown4", vertex.color="brown4", 
     vertex.size=(degree(gc)/max(degree(gc)))*4+0.8, vertex.label=NA )
title(main=list("Degree>20", col="white"))
par(mfrow=c(1,1))

#rich club community
com<-fastgreedy.community (rcgc20, merges=TRUE, modularity=TRUE, weights=NULL)
com         #groups: 13, mod: 0.73                   
modularity(com) #0.7273917                    
membership(com)

par(bg='black')
gcom<-rainbow(max(com$membership))
plot(rcgc20, layout=lay20,asp=1, margin=0.2, edge.color="grey45", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gcom[as.numeric(membership(com))],
     vertex.color=gcom[as.numeric(membership(com))], vertex.size=(degree(gc)/max(degree(gc)))*4+0.6, 
     vertex.label=NA)
title(list("rich club 20 communities", col="white"))

# για το 2,5%
vcount(g)*0.025  # 131
rev(sort(deg))[131]  # 32

#και γραφικά
par(bg="black")
par(mfrow=c(1,2))
rcgc32<-induced.subgraph(gc,V(gc)[degree(gc)>=32])
rcgc32    #132 2280
lay32<-layout_with_graphopt(rcgc32)
plot(rcgc32, layout=lay32, edge.color="lightcyan4", edge.width=0.05, 
     vertex.shape="circle", vertex.frame.color="brown4", vertex.color="brown4", 
     vertex.size=(degree(gc)/max(degree(gc)))*4+0.8, vertex.label=NA)
title(main=list("Degree>=32", col="white"))

rcgc32<-induced.subgraph(gc,V(gc)[degree(gc)>32])
rcgc32     #129 2278
lay32<-layout_with_graphopt(rcgc32)
plot(rcgc32, layout=lay32, edge.color="lightcyan4", edge.width=0.05, 
     vertex.shape="circle", vertex.frame.color="brown4", vertex.color="brown4", 
     vertex.size=(degree(gc)/max(degree(gc)))*4+0.8, vertex.label=NA)
title(main=list("Degree>32", col="white"))
par(mfrow=c(1,1))

#rich club communities again
com<-fastgreedy.community (rcgc32,merges=TRUE, modularity=TRUE, weights=NULL)
com         #groups: 8, mod: 0.64                   
modularity(com) #0.6429877                    
membership(com)

par(bg='black')
gcom<-rainbow(max(com$membership))
plot(rcgc32, layout=lay32,asp=1, margin=0.2, edge.color="grey45", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gcom[as.numeric(membership(com))],
     vertex.color=gcom[as.numeric(membership(com))], vertex.size=(degree(gc)/max(degree(gc)))*4+0.6, 
     vertex.label=NA)
title(list("rich club 32 communities", col="white"))




# σε κάθε περίπτωση καθως το πληθος των κορυφων και κυριως των ακμων είναι πολύ μεγάλο
# δεν μπορεί γραφικά να γίνει καλή εκτίμηση, αλλά μοιάζει να ισχύει το φαινόμενο rich club.

# ενναλακτικά μπορούμε να δούμε τον συντελεστή rich club του δικτύου
# μέσω της εντολής rich_club_coef(g,k=1, weighted=NULL)
# όπου k ο βαθμός που θέλω να είναι μεγαλύτερος
# της βιβλιοθήκης brainGraph
# https://en.wikipedia.org/wiki/Rich-club_coefficient
# φ(k)=2Εk/(Nk(Nk-1))
# όπου φ ο συντελεστής, Εk ο αριθμός των ακμών που συνδέουν κορυφές βαθμού μεγαλύτερου 
# του k, και Νkο αριθμός των κορυφών βαθμού μεγαλύτερου του k
# εντολή επιστρέφει τον συντελεστή και το γράφημα με τις αντίστοιχες κορυφές και ακμές
library(brainGraph)
rc20<-rich_club_coeff(g, k=20, weighted = NULL)
rc20$phi    # 0.09930523
rc32<-rich_club_coeff(g,32,weighted =NULL )
rc32$phi    # 0.2759205


#rich_club_coeff(g,k=40,weighted = NULL)
#$phi
#[1] 0.5019026
#rich_club_coeff(g,k=50,weighted=NULL)
#$phi
#[1] 0.9090909
# όσο λοιπόν αυξάνεται ο βαθμός τελικά, τόσο αυξάνεται και ο συντελεστής
# μπορούμε λοιπόν πλέον με μεγαλύτερη σιγουριά να υποθέσουμε ότι το φαινόμενο 
# rich club ισχύει.

# ιδανικά θα υπολογίζαμε το normalized rich club coefficient
# Ομως τα τυχαια δίκτυα που δημιουργούνται ακολουθούν κατανομή Poisson στους βαθμούς
# και δεν εμφανίζουν τόσο μεγάλους βαθμούς όπως κ=20,32 κλπ
# έτσι ο συντελεστής δίνει $phi =ΝaN



##########################################################################
###  8. Εύρεση κοινοτήτων μέσα στο δίκτυο (ή στη γιγάντια συνιστώσα)   ###
##########################################################################

#δεν θα χρησιμοποιηθούν ολοι οι αλγόριθμοι εύρεσης κοινοτήτων της igraph
#καθώς ορισμένοι είτε χρειάζονταν πολύ χρόνο, ή πολύ μνήμη, είτε αναφέρονται σε 
#σταθμισμένα δίκτυα.
#από τους υπόλοιπους επιλέχτηκαν 4 καθώς το δίκτυο είναι μεγάλο και τα αποτελέσματα 
#γραφικά δεν φαίνονταν καλά με κανέναν.
#Επιπλέον αν και οι κορυφές του γραφήματος είναι ονοματισμένες, αυτό έχει γίνει με αριθμούς
#που ο καθένας αντιστοιχεί σε κάποιον ερευνητή αλλά δεν έχω πρόσβαση σε αυτήν τη γνώση
#Ούτε έχω κάποια άλλη επιπλέον πληροφορία για το δίκτυο ώστε να μπορώ να δικαιολογήσω κάποια κοινότητα

# οι δύο πρωτοι αλγόριθμοι για τη γιγάντια συνιστώσα
#1.Community structure via greedy optimization of modularity
# This function tries to find dense subgraph, also called communities in graphs
# via directly optimizing a modularity score.

com<-fastgreedy.community (gc,merges=TRUE, modularity=TRUE, weights=NULL)
com         #groups: 74, mod: 0.8                   
modularity(com) #0.7997737                    
membership(com)

# default plot για communities
par(bg="white")
plot(com, gc, vertex.label=NA)
# δεν δίνει καλά αποτελέσματα και δεν θα χρησιμοποιήται

# Plot αλλιώς
par(bg='black')
gcom<-rainbow(max(com$membership))
plot(gc, layout=laygc,asp=1, margin=0.2, edge.color="grey45", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gcom[as.numeric(membership(com))],
     vertex.color=gcom[as.numeric(membership(com))], vertex.size=(degree(gc)/max(degree(gc)))*4+0.6, 
     vertex.label=NA, xlim=c(-0.2, 0.9), ylim=c(-0.4, 0.4))

	 
#2. Finding communities in graphs based on statistical meachanic
#This function tries to find communities in graphs
#via a spin-glass model and simulated annealing.

cs<-cluster_spinglass(gc, weights=NA) 
cs       #groups: 25, mod: 0.84
modularity(cs)				#0.839957

gcs<-rainbow(max(cs$membership))
plot(gc, layout=laygc, asp=1, margin=0.2, edge.color="grey45", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gcs[as.numeric(membership(cs))],
     vertex.color=gcs[as.numeric(membership(cs))], vertex.size=(degree(gc)/max(degree(gc)))*4+0.6, 
     vertex.label=NA, xlim=c(-0.2, 0.9), ylim=c(-0.4, 0.4))
	 
#γραφικά μαζί
par(mfrow=c(1,2))
gcom<-rainbow(max(com$membership))
plot(gc, layout=laygc,asp=1, margin=0.2, edge.color="grey45", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gcom[as.numeric(membership(com))],
     vertex.color=gcom[as.numeric(membership(com))], vertex.size=(degree(gc)/max(degree(gc)))*4+0.6, 
     vertex.label=NA, xlim=c(-0.2, 0.9), ylim=c(-0.4, 0.4))
title(main=list("FastGreedy for gc", col="white", cex=1))

plot(gc, layout=laygc, asp=1, margin=0.2, edge.color="grey45", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gcs[as.numeric(membership(cs))],
     vertex.color=gcs[as.numeric(membership(cs))], vertex.size=(degree(gc)/max(degree(gc)))*4+0.6, 
     vertex.label=NA, xlim=c(-0.2, 0.9), ylim=c(-0.4, 0.4))
title(main=list("ClusterSpinglass for gc", col="white", cex=1))
par(mfrow=c(1,1))
	 

#Οι δύο επόμενοι αλγόριθμοι για ολόκληρο το δίκτυο
#3.Community strucure via short random walks
# This function tries to find densely connected subgraphs,
# also called communities in a graph via random walks. 
# The idea is that short random walks tend to stay in the same community.

wc <- cluster_walktrap(g,weights = NULL) 
wc <- cluster_walktrap(g,weights = NULL) 
wc		#groups: 815, mod: 0.78
modularity(wc)		#0.7823644
membership(wc)

gwc<-rainbow(max(wc$membership))
plot(g, layout=layg,asp=1, margin=0.2, edge.color="grey45", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gwc[as.numeric(membership(wc))],
     vertex.color=gwc[as.numeric(membership(wc))], vertex.size=(degree(g)/max(degree(g)))*4+0.6, 
     vertex.label=NA, xlim=c(-0.2, 0.9), ylim=c(-0.4, 0.4))



#4. Finding community structure by multi-level optimization of modularity
# This function implements the multi-level modularity optimization algorithm for finding
# community structure. It is based on the modularity measure and a hierarchial approach
cl<-cluster_louvain(g)
cl    #groups: 395, mod: 0.86
modularity(cl)    # 0.8609414

gcl<-rainbow(max(cl$membership))
plot(g, layout=layg, asp=1, margin=0.2, edge.color="lightsteelblue4", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gcl[as.numeric(membership(cl))],
     vertex.color=gcl[as.numeric(membership(cl))], vertex.size=(degree(g)/max(degree(g)))*4+0.6, 
     vertex.label=NA, xlim=c(-0.2, 0.9), ylim=c(-0.4, 0.4))

#γραφικά μαζί
par(mfrow=c(1,2))
plot(g, layout=layg,asp=1, margin=0.2, edge.color="grey45", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gwc[as.numeric(membership(wc))],
     vertex.color=gwc[as.numeric(membership(wc))], vertex.size=(degree(g)/max(degree(g)))*4+0.6, 
     vertex.label=NA, xlim=c(-0.2, 0.9), ylim=c(-0.4, 0.4))
title(main=list("ClusterWalktrap for g", col="white", cex=1))

plot(g, layout=layg, asp=1, margin=0.2, edge.color="lightsteelblue4", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gcl[as.numeric(membership(cl))],
     vertex.color=gcl[as.numeric(membership(cl))], vertex.size=(degree(g)/max(degree(g)))*4+0.6, 
     vertex.label=NA, xlim=c(-0.2, 0.9), ylim=c(-0.4, 0.4))
title(main=list("ClusterLouvain for g", col="white", cex=1))
par(mfrow=c(1,1))


# και περνάω σε σύγκριση μεταξύ των αποτελεσμάτων κάθε ζευγαριού αλγορίθμων
# έυρεσης κοινοτήτων με το δείκτη NMI (Normalized Mutual Information)

compare(com,cs,method = "nmi") # 0.5509663
compare(cl,wc,method="nmi")    # 0.7904059
# και τα δύο νούμερα που προκείπτουν ερμηνεύονται ως ομοιότητα στις κοινότητες
# στη δεύτερη περίπτωση που αφορόυσε ολόκληρο το δίκτυο η ομοιότητα είναι ακόμα μεγαλύτερη
# παρόλου που η διαφορα σε πλήθος κοινοτήτων ήταν πάρα πολύ μεγάλη ~500
dev.off()

###############################################################################################
###############################################################################################



