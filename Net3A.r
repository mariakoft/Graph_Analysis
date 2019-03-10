## Erdos net- Pajek
# Βακρατσά Θεοδώρα

setwd("C:/Users/xx/Desktop/diktyo3")## Bale to path tou fakelou sto opoio exeis to dataset
Sys.setlocale("LC_CTYPE","Greek")
library(igraph)
library(Matrix)
library(lattice)
library(vcd)

#########Open Network File################3
g<-read_graph(file="Erdos02.net", format="pajek")
summary(g)
V(g)$name[1:10]
#εισαγωγικά
is.named(g)     #true
is.directed(g)  #false
is.connected(g) #true ,είναι συνδετικό δίκτυο
is.simple(g)    #true, δεν έχει λούπες

vcount(g)    #6927 κορυφές
ecount(g)    #11850 ακμές

graph.density(g)    #0.0004939...πολύ αραιό

#και θέλω να το δω αρχικά γραφικά
layg<-layout_with_fr(g)
plot(g, layout=layg, edge.size=0.5, edge.color="grey67", vertex.shape="circle", vertex.size=1,
     vertex.frame.color="firebrick4",vertex.color="firebrick4", vertex.label=NA)

#zoom in  sto plot
plot(g,layout=layg, edge.size=0.5, edge.color="grey67", vertex.shape="circle", vertex.size=1,
     vertex.frame.color="firebrick4",vertex.color="firebrick4", vertex.label=NA,
     xlim=c(-0.2, 0.2), ylim=c(0.3,0.6))

#δεν χρειάζεται να ψάξω για συνιστώσες επειδή το δίκτυό μου είναι συνδετικό
#αρα αποτελείται από μία συνιστώσα

#########################################################################################
#περνάμε οπότε στον έλεγχο των βαθμών του δικτύου
summary(degree(g))    #πολύ λοξή κατανομή
mean(degree(g))   #3.421395
sd(degree(g))     #11.68919
sd(degree(g))/mean(degree(g))  #συντελεστής μεταβλητότητας 341,6499%

#φαίνεται ότι μοιάζει με δυναμοκατανομή, θα τη δω γραφικά και θα ελέγξω
de<-degree(g)
D<-max(de); D  #507

#και βρίσκω τις κατανομές
dd<-degree.distribution(g,cumulative =FALSE )   #συνάρτηαση πυκτότητας
cdd<-degree.distribution(g,cumulative = TRUE)   #cumulative degree distribution

#και θα τις δω γραφικά
par(mfrow=c(1,3))
plot(dd,xlab = "degree", ylab="Frequency")
plot(dd,xlab = "Log degree", ylab = "Log Freq", log="xy")
plot(cdd, xlab = "Log degree", ylab="Log CDD", log="xy")
par(mfrow=c(1,1))

#μοιάζει να είναι power law κατανομή
plf1<-power.law.fit(de,xmin=NULL,start=2, implementation = "plfit",force.continuous = FALSE)
plf1$alpha
plf1$xmin
plf1$KS.p   #1.356228e-05<0.05 άρα δεν ακολουθει power law


#############################################################################
#επόμενος έλεγχος θα ήταν τα βάρη των ακμών στο δίκτυο E(g2)$weight
#αλλά δεν έχω δίκτυο με βάρη
is.weighted(g2)    #FALSE


#διάμετρος του δικτύου
dg<-diameter(g, unconnected = FALSE, weights = NA) ; dg   #4 
                                                          #λογικό, εξαιτίας της κατασκευής του δικτύου
#μέσο μήκος μονοπατιού
aplg<-average.path.length(g, unconnected = FALSE, directed = FALSE); aplg   #3.776441

#θα δω και την κατανομή των αποστάσεων
plhg<-path.length.hist (g, directed = FALSE) 
#sum(plhg2$res)=vcount(g)*(vcount(g)-1)/2 =23988201

#γραφικά
plot(1:dg, plhg$res, xlab="distance", ylab="paths",frame=TRUE)
#kapoio sxolio?


#################################################################
# μεταβατικότητα, συντελεστής σύμπλεξης
#παρόλο που συχνά συγχέονται αποτελούν δύο διαφορετικά μέτρα
#η πρώτη, σε ένα δύκτυο συνεργασιών, όπως αυτό εδώ μετράει κατά πόσο αν οι B και Γ έχουν συνεργαστεί
#με τον Α συνεργάζονται και μεταξύ τους, δηλαδη υπαρχει η ακμή {Β,Γ}
#Για ολόκληρο το γράφημα G ορίζεται ως:
#T(G)=πλήθος μεταβατικών τριάδων/πλήθος δυνατών συνδετικών τριάδων
#    =3x πλήθος τριγώνων του G/πλήθος δυνατών συνδετικών τριάδων
#    =3τ(G)/ρ(G)
transitivity(g, type="global", isolates= "zero")     #0.03570461


#ο συνετελεστής σύμπλεξης για μία κορυφή ν του δικτύου μετρά κατά πόσον οι γείτονές της
#απέχουν από το να αποτελούν κλίκα
#ορίζεται ως c(ν)=πλήθος τριγώνων με κορυφή το v/πλήθος τριάδων με κορυφή το ν
#και για το δίκτυο C(ν)=Σc(ν)/Πλήθος κορυφών δικτύου

ltg<-transitivity(g, type="local", isolates= "zero")  #o σ.σ. για κάθε κορυφή
ltg    #δεν έχει NaN τιμές
ccg<-sum(ltg)/vcount(g); ccg    # 0.1239001

#για να βρω ιεραρχική δομή, θα το δω αρχικά γραφικά
#τοποθετώ τις παραπάνω τιμές σε πίνακα
ltm<-matrix(c(degree(g), ltg), nrow=vcount(g), ncol=2, byrow=FALSE)
#"διώχνω" μηδενικές τιμές
ltm0<-ltm[ltm[,2]!= 0, ] ;   min(ltm0)   #0.001403509


# scatterplot of degree vs CC
plot(ltm0[,1], ltm0[,2], log='xy', xlab="Log degree", ylab="Log CC",frame=TRUE) 

# do the linear regression
lr<-lm(log(ltm0[,2]) ~ log(ltm0[,1]) +1) 

summary(lr)
#constant (=intercept)=0.64347 και 
#slope (=log(ltm0[,1]))=-1.10142

#και γραφικά το scatterplot μαζί με τη γραμμική παλινδρόμηση
plot(ltm0[,1], ltm0[,2], log='xy', xlab="degree", ylab="CC",frame=TRUE);
lines(ltm0[,1], exp(0.6434722  -1.1014233*log(ltm0[,1])))  
#παρατηρώ οτι κορυφές μεγάλου βαθμού έχουν μικρό συντελεστή σύμπλεξης
#άρα δεν αποτελούν μέλη κλικας


########################################################################
### 5. Συντελεστής ομοιότητας μεταξύ των κορυφών (ως προς το βαθμό)  ###
########################################################################

# Μέσω του συντελεστή ομοιότητας (ως προς το βαθμό) υπολογίζουμε το ποσοστό των κορυφών μεγάλου
# βαθμού που συνδέονται με κορυφές όμοιου βαθμού
# θα δούμε αν το δίκτυο είναι assortative ή disassortative γραφικά
assortativity.degree(g)   #  -0.1155774 disassortative
# άρα τελικά οι κορυφές μεγάλου βαθμού τείνουν να μην συνδέονται με όμοιές τους κορυφές

# υπολογισμός του ANND(average nearest neighbor degree)
# ο αλγοριθμος υπολογίζει το μέσο βαθμό των γειτόνων κάθε κορυφής και το μέσο βαθμό των
# γειτόνων κορυφών ίσου βαθμού
knng<-graph.knn(g,weights = NA)
knng$knnk   # πάρα πολλές τιμές NaN
knnkg<-na.omit(knng$knnk)
degreesg0<-as.numeric(names(table(degree(g))))
degreesg<-degreesg0[degreesg0!=0]

#γραμμική παλινδρόμησηση
fitknn<-lm(knnkg~degreesg+1)
summary(fitknn)
# constant=38.66231, slope=-0.11101 

#γραφικά
plot(degreesg,knnkg, xlab="degree", ylab="knnkg",frame=TRUE); 
lines(degreesg, 38.66231 -0.11101 *degreesg )      
# δηλαδή κορυφές μεγάλου βαθμού τείνουν να συνδέονται περισσότερο με κορυφές 
# με μικρότερο μέσο βαθμό


#######################################################################################
### 6. Το φαινόμενο Rich Club για το 5% και το 2,5% των κορυφών μεγαλύτερου βαθμού  ####
########################################################################################

# ουσιαστικά ψάχνω να βρω αν υπάρχουν ακμές στο δίκτυο μου μεταξύ κορυφών μεγάλου βαθμού
# και έτσι δημιουργήται η εικόνα "Rich Club"
# ελέγχω στη γιγάντια σινιστώσα


# _Πρέπει αρχικά να βρώ κορυφές ποιου βαθμού και πάνω αποτελούν το 5% των πιο υψηλοβαθμων κορυφών
deg<-degree(g)
head(sort(deg)); tail(sort(deg))

vcount(g)*0.05    # 346
# άρα θελώ να πάρω κορυφές βαθμού μεγαλύτερου ή ίσου του βαθμού της 208ης κορυφής
# αν τοποθετήσω τους βαθμούς των κορυφών σε φθίνουσα σειρά
rev(sort(deg))[346]   # 12
#δηλαδή το 5% των κορυφών μεγαλύτερου βαθμού, έχουν βαθμό >= 12

#και για να φανεί γραφικά
rcg12<-induced.subgraph(g,V(g)[degree(g)>=12])
rcg12
lay12<-layout_with_fr(rcg12)
plot(rcg12, layout=lay12, edge.color="grey67", edge.width=0.25, 
     vertex.shape="circle", vertex.frame.color="brown4", vertex.color="brown4", 
     vertex.size=(degree(g)/max(degree(g)))*10+2, vertex.label=NA, main="βαθμός(g)>=12")


# για το 2,5%
vcount(g)*0.025  # 173
rev(sort(deg))[173]  # 27

#και γραφικά
rcg27<-induced.subgraph(g,V(g)[degree(g)>=27])
rcg27
lay27<-layout_with_fr(rcg27)
plot(rcg27, layout=lay27, edge.color="grey67", edge.width=0.25, 
     vertex.shape="circle", vertex.frame.color="brown4", vertex.color="brown4", 
     vertex.size=(degree(g)/max(degree(g)))*10+2, vertex.label=NA, main="βαθμος>=27")

# υπάρχουν πολλές συνδέσει και μπορώ να υποθέσω οτι το φαινόμενο rich club επαληθεύεται
# όμως μπορεί να φαίνονται εξαιτίας του Paul Erdos στο δίκτυο
# θα δούμε τα γραφήματα αφού τον αφαιρέσουμε
gwithout<-delete.vertices(g, 6927)
par(mfrow=c(1,2))
rcg12<-induced.subgraph(gwithout,V(gwithout)[degree(gwithout)>=12])
lay12<-layout_with_fr(rcg12)
plot(rcg12, layout=lay12, edge.color="grey67", edge.width=0.25, 
     vertex.shape="circle", vertex.frame.color="brown4", vertex.color="brown4", 
     vertex.size=(degree(gwithout)/max(degree(gwithout)))*10+2, vertex.label=NA)
title(list("Rich club 12 w/out Erdos", col="brown4"))


rcg27<-induced.subgraph(gwithout,V(gwithout)[degree(gwithout)>=27])
lay27<-layout_with_fr(rcg27)
plot(rcg27, layout=lay27, edge.color="grey67", edge.width=0.25, 
     vertex.shape="circle", vertex.frame.color="brown4", vertex.color="brown4", 
     vertex.size=(degree(gwithout)/max(degree(gwithout)))*10+2, vertex.label=NA)
title(list("Rich club 27 w/out Erdos", col="brown4"))
par(mfrow=c(1,1))


#σε κάθε περίπτωση φαίνονται πολλές κορυφές υψηλού βαθμού να έχουν μεταξύ τους ακμές
#άρα το φαινόμενο rich club υπάρχει


##########################################################################
######### 8. Ύπαρξη κοινοτήτων μέσα στο δίκτυο (ή στη γιγάντια συνιστώσα)
##########################################################################


com<-fastgreedy.community (g,merges=TRUE, modularity=TRUE, weights=NULL)
com                                          ## groups: 58, mod: 0.67
modularity(com)                              ## Q=0.6735519

## PLOT the communities using igraph function
plot(com, g, vertex.label=NA)							


## plot διαφορετικά
gcom<-rainbow(58)
plot(g, layout=layg,asp=1, margin=0.2, edge.color="grey35", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gcom[as.numeric(membership(com))],
     vertex.color=gcom[as.numeric(membership(com))], vertex.size=(degree(g)/max(degree(g)))*5, 
     vertex.label=NA, xlim=c(0.0, 0.9), ylim=c(-0.3, 0.3))


#This function tries to find densely connected subgraphs,
#also called communities in a graph via random walks. 
#The idea is that short random walks tend to stay in the same community.
# Larger edge weights correspond to stronger connections
wc <- cluster_walktrap(g) # ignore direction, get edge weights
wc		#groups: 244, mod: 0.63
modularity(wc)		#0.6271699
plot(wc, g, vertex.label=NA)
members4<-as.numeric(membership(wc))
V(g)$wc<- members4
assortativity_nominal(g, types=V(g)$wc)  #0.6857343

gwc<-rainbow(244)
par(bg="black", new=FALSE)
plot(g, layout=layg,asp=1, margin=0.2, edge.color="lightsteelblue4", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gwc[as.numeric(membership(wc))],
     vertex.color=gwc[as.numeric(membership(wc))], vertex.size=(degree(g)/max(degree(g)))*5+0.6, 
     vertex.label=NA, xlim=c(0.0, 0.9), ylim=c(-0.3, 0.3))



#This function tries to find communities in graphs
#via a spin-glass model and simulated annealing.
# Larger edge weights correspond to stronger connections
cs<-cluster_spinglass(g) # ignore directionm, get edge weights oups
cs       #groups: 25, mod: 0.71
modularity(cs)				#0.7102484
par(bg="white",new=FALSE)
plot(cs, g, vertex.label=NA)
members5<-as.numeric(membership(cs))
V(g)$cs<- members5
assortativity_nominal(g, types=V(g)$cs)		#0.7850835
###2nd plot for cs###
gcs<-rainbow(25)
par(bg="black", new=FALSE)
plot(g, layout=layg,asp=1, margin=0.2, edge.color="lightsteelblue4", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gcs[as.numeric(membership(cs))],
     vertex.color=gcs[as.numeric(membership(cs))], vertex.size=(degree(g)/max(degree(g)))*8+0.6, 
     vertex.label=NA, xlim=c(0.0, 0.9), ylim=c(-0.3, 0.3))



ci<-cluster_infomap(g) # weights on edges and vertices
ci     # groups: 314, mod: 0.62
modularity(ci)   #0.6203635		 
par(bg="white", new=FALSE)
plot(ci, g, vertex.label=NA)
members7<-as.numeric(membership(ci))
V(g)$ci<- members7
assortativity_nominal(g, types=V(g)$ci)   #0.6252543

gci<-rainbow(314)
par(bg="black", new=FALSE)
plot(g, layout=layg,asp=1, margin=0.2, edge.color="lightsteelblue4", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gci[as.numeric(membership(ci))],
     vertex.color=gci[as.numeric(membership(ci))], vertex.size=(degree(g)/max(degree(g)))*8+1, 
     vertex.label=NA, xlim=c(0.0, 0.9), ylim=c(-0.3, 0.3))



# This is a fast, nearly linear time algorithm for 
#detecting community structure in networks. In works 
#by labeling the vertices with unique labels and then 
#updating the labels by majority voting in the neighborhood 
# of the vertex.
# Larger edge weights correspond to stronger connections.
clp<-cluster_label_prop(g) # weighted but undirected
clp    # groups: 725, mod: 0.79
names(clp)
clp$algorithm
modularity(clp)    #0.7949811
membership(clp)
sizes(clp)
par(bg="white", new=FALSE)
plot(clp, g) 
members8<-as.numeric(membership(clp))
V(g)$clp<- members8
assortativity_nominal(g, types=V(g)$clp)

gclp<-rainbow(725)
par(bg="black", new=FALSE)
plot(g, layout=layg,asp=1, margin=0.2, edge.color="lightsteelblue4", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gclp[as.numeric(membership(clp))],
     vertex.color=gclp[as.numeric(membership(clp))], vertex.size=(degree(g)/max(degree(g)))*8+1, 
     vertex.label=NA, xlim=c(0.0, 0.9), ylim=c(-0.3, 0.3))


#This function implements the multi-level modularity 
#optimization algorithm for finding community structure, 
#see references below. It is based on the modularity measure 
#and a hierarchial approach
cll<-cluster_louvain(g)
cll    #groups: 395, mod: 0.86
names(cll)
cll$algorithm     
modularity(cll)    #0.8609414
membership(cll)
sizes(cll)
par(bg="white", new=FALSE)
plot(cll, g, vertex.label=NA)
members9<-as.numeric(membership(cll))
V(g)$cll<- members9
assortativity_nominal(g, types=V(g)$cll)   #0.8977026

gcll<-rainbow(395)
par(bg="black", new=FALSE)
plot(g, layout=layg,asp=1, margin=0.2, edge.color="lightsteelblue4", edge.width=0.01, 
     vertex.shape="circle", vertex.frame.color=gcll[as.numeric(membership(cll))],
     vertex.color=gcll[as.numeric(membership(cll))], vertex.size=(degree(g)/max(degree(g)))*8+1, 
     vertex.label=NA, xlim=c(0.0, 0.9), ylim=c(-0.3, 0.3))



