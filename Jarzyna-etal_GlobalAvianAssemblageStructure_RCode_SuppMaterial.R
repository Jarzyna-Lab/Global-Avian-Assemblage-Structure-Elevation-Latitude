##### "Global functional and phylogenetic structure of avian assemblages across elevation and latitude"
##### Supplementary Material (R code)
##### Direct all correspondence to Marta A. Jarzyna at jarzyna.1@osu.edu

require(data.table)
require(FD)
require(ape)
require(phytools)
require(picante)

###################################################
###################################################STEP 1: Functional & phylogenetic diversity calculation
###################################################
#### Read in data on species distributions
adt <- load(file="data/points_all_species.rda") #dowload from https://mol.org/downloads/
adt.mo <- subset(adt, adt$moId != "NA")


#### Dendrogram-based functional diversity and distinctness
## Create master functional tree (load trait data for all birds: EltonTraits 1.0, 
## available at https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/13-1917.1)
sp.traits <- read.csv("data/Master_Funcdiv.csv", header=TRUE)
rownames(sp.traits) <- sp.traits[,2]
sp <- as.character(sp.traits[,2])
sp.ids <- sp.traits[,1:3]

## Select traits of interest
funcdat <- subset(sp.traits, select = -c(WJSpecID,Scientific,English,Diet_Breadth_Levin, Diet_breadthLevinsStanderdized,log_dietBreathLevin))
## Assign weights to traits, the entire dietary (and foraging) axis is to be weighted equally to body mass, pelagic, or nocturnality 
Weights=c(rep(0.2/7,7),rep(0.2/7,7),0.2,0.2,0.2) 
 
## Create functional dissimilarity matrix using gower distance and weights
funcdist <- gowdis(funcdat, ord = "podani", w = Weights)
 
## Create dengrogram based on UPGMA clustering method
SpecUPGMA <- hclust(funcdist,  method = "average")
SpecUPGMA$labels <- sp
functree <- as.phylo(SpecUPGMA)

## Quantify functional divesity (FD) and mean local FD distinctness (FDS) for each point location
poIDs <- unique(adt.mo$poId)
FD <- matrix(NA,length(poIDs),1)
FDS <- list()
FDSm <- matrix(NA,length(poIDs),5)

for (k in 1:length(poIDs)) {
  comm <- subset(adt.mo, poId == poIDs[k])
  nspec <- nrow(comm)
  comm <- setDF(comm)
  
  sp.id <- SpecUPGMA$labels
  rmTip <- as.character(sp.id[! sp.id %in% comm[,1]]) #Remove species tips for species that are not present
  subTree <- drop.tip(functree, tip=rmTip)  # Remove those branches
  FD[k,] <- sum(subTree$edge.length) # FD = Summing up branch lengths
  
  fddist <- evol.distinct(subTree, type = "fair.proportion",scale = FALSE, use.branch.lengths = TRUE) # functional distinctness
  
  OrigOrder.df <- data.frame(Species=comm[,1])
  merged <- merge(OrigOrder.df, fddist, sort=FALSE, all.x=TRUE)
  merged$Species <- factor(merged$Species, level=OrigOrder.df$Species)	
  merged <- with(merged, merged[order(Species),])
  dist.mat <- merged[,2]
  
  FDS[[k]] <- dist.mat
  FDSm[k,1] <- mean(dist.mat) # mean FDS for the assemblage
  FDSm[k,2] <- psych::geometric.mean(dist.mat)
  FDSm[k,3] <- quantile(dist.mat, probs = 0.05)
  FDSm[k,4] <- quantile(dist.mat, probs = 0.5)
  FDSm [k,5] <- quantile(dist.mat, probs = 0.95)
}

FDid <- as.data.frame(cbind(FD,poIDs))
FDid$longitude <- 0
FDid$latitude <- 0
FDid$elevation <- 0
FDid$moIDs <- 0

FDSid <- as.data.frame(cbind(FDSm,poIDs))
FDSid$longitude <- 0
FDSid$latitude <- 0
FDSid$elevation <- 0
FDSid$moIDs <- 0

for (k in 1:length(poIDs)){
  comm <- subset(adt.mo, poId == poIDs[k])
  FDid[k,3:6] <- comm[1,3:6]
  FDSid[k,6:9] <- comm[1,3:6]
}

colnames(FDid) <- c("FD","poIDs","longitude","latitude","elevation","moIDs")
saveRDS(FDid, file="FD.rds")
colnames(FDdistid) <- c("mean","gmean","median","lci","uci","poIDs","longitude","latitude","elevation","moIDs")
saveRDS(FDdistid, file="FDSm.rds")


#### Hypervolume-based functional diversity
sp.traits <- read.csv("data/Master_Funcdiv.csv", header=TRUE)
rownames(sp.traits) <- sp.traits[,1]
sp <- as.character(sp.traits[,1])
sp.ids <- sp.traits[,1:3]
funcdat <- subset(sp.traits, select = -c(WJSpecID,Scientific,English,Diet_Breadth_Levin, Diet_breadthLevinsStanderdized,log_dietBreathLevin))

## Principal Coordinate Analysis
pcoa.calc <- pcoa(funcdist, correction="lingoes")
nbdim <- ncol(pcoa.calc$vectors)
eigen <- pcoa.calc$values
vect <- pcoa.calc$vectors[,1:nbdim]
rownames(vect) <- rownames(funcdat)
colnames(vect) <- paste("PC",1:nbdim,sep="")

FDH <- matrix(NA,length(poIDs),3)

for (i in 1:length(poIDs)){
  comm <- adt.mo[adt.mo$poId ==  poIDs[i],]
  sp.i <- unique(comm$spId)
  vect.i <- vect[rownames(vect) %in% sp.i,1:2]
  fdh.i <- hypervolume::hypervolume(vect.i, method = "gaussian")
  saveRDS(fdh.i, file=paste0("fdh_2D_",poIDs[i],".rds"))
  
  centr.i <- hypervolume::get_centroid(fdh.i) # Centroid of the hypervolume
  vol.i <- hypervolume::get_volume(fdh.i) # Volume of the hypervolume
  FDH[i,1:2] <- centr.i
  FDH[i,3] <- vol.i
}

FDHid <- as.data.frame(cbind(FDH,poIDs))
FDHid$longitude <- 0
FDHid$latitude <- 0
FDHid$elevation <- 0
FDHid$moIDs <- 0

for (k in 1:length(poIDs)){
  comm <- subset(adt.mo, poId == poIDs[k])
  FDHid[k,5:8] <- comm[1,3:6]
}

colnames(FDHid) <- c("FDH.centr1","FDH.centr2","FDH.vol","poIDs","longitude","latitude","elevation","moIDs")
saveRDS(FDHid, file="FDHm.rds")


#### Phylogenetic divesity (PD)
phylotree <- read.nexus("output.nex") #downloaded from birdtree.org
x1 <- gsub('([[:punct:]])|\\s+','_',adt.mo$spp)
adt.mo <- cbind(adt.mo, x1)

poIDs <- unique(adt.mo$poId)
PD <- matrix(NA,length(poIDs),100) # Phylogenetic diversity
PDSm <- matrix(NA,length(poIDs),100) # Mean phylogenetic distinctness
PDS <- list() # To hold species-level phylogenetic distinctness
PDSi <- list()

for (i in 1:100){	
  pdi <- functree[[i]]
  
  for (k in 1:length(poIDs)) {
    comm <- subset(adt.mo, poId == poIDs[k])
    nspec <- nrow(comm)
    comm <- setDF(comm)
    sp.id <- pdi$tip.label
    
    rmTip <- as.character(sp.id[! sp.id %in% comm$x1]) # Remove missing tips 
    subTree <- drop.tip(pdi, tip=rmTip)  # remove those branches
    PD[k,i] <- sum(subTree$edge.length) # Sum the branch length
    di <- evol.distinct(subTree, type = "fair.proportion",scale = FALSE, use.branch.lengths = TRUE)
    PDSm[k,i] <- psych::geometric.mean(di[,2])
    PDSi[[k]] <- di[,2]
  }
  
  PDS[[i]] <- PDSi
}


PDid <- as.data.frame(cbind(PD,poIDs))
PDid$longitude <- 0
PDid$latitude <- 0
PDid$elevation <- 0
PDid$moIDs <- 0
PDSid <- as.data.frame(cbind(PDSm,poIDs))
PDSid$longitude <- 0
PDSid$latitude <- 0
PDSid$elevation <- 0
PDSid$moIDs <- 0

for (k in 1:length(poIDs)){
  comm <- subset(adt.mo, poId == poIDs[k])
  PDid[k,102:105] <- comm[1,3:6]
  PDSid[k,102:105] <- comm[1,3:6]
}

saveRDS(PDid, file="PD.rds")
saveRDS(PDSid, file="PDSm.rds")
saveRDS(PDS, file="PDS.rds")



###################################################
###################################################STEP 2: Null models (Example calculation for FD)
###################################################
#### Randomize species presences
sp.all <- load(file="data/points_all_species.rda")
adt.mo <- subset(adt, adt$moId != "NA")
moIDs <- unique(adt.mo$moId)
poIDs <- unique(adt.mo$poId)
require(dplyr)

for (i in 1:100){
  comm.null.all <- list()
  for (k in 1:length(moIDs)) {
    comm.null.mo <- list()
    comm.mo <- subset(adt.mo, moId == moIDs[k])
    comm.mo.sub <- comm.mo %>%
      distinct(spp, .keep_all= TRUE)
    poIDs.mo <- unique(comm.mo$poId)
    
    for (j in 1:length(poIDs.mo)) {
      comm.po <- subset(comm.mo, poId == poIDs.mo[j])
      nspec <- nrow(comm.po)
      comm.null <- sample_n(comm.mo.sub, nspec)
      comm.null[2:7] <- comm.po[,2:7]
      comm.null.mo <- rbind(comm.null.mo,comm.null)
    }
    comm.null.all <- rbind(comm.null.all,comm.null.mo)
  }
  saveRDS(comm.null.all,file=paste0("nullout/null_mountsp_",i,"_test.rds"))
}


#### Calculate null expectation of functional diversity
poIDs <- unique(adt.mo$poId)
FDnull <- matrix(NA,length(poIDs),100)

for (i in 1:100){
  n1 <- readRDS(file=paste0("nullout/null_mountsp_",i,".rds"))
  
  for (k in 1:length(poIDs)) {
    comm <- subset(adt.mo, poId == poIDs[k])
    nspec <- nrow(comm)
    comm <- setDF(comm)
    
    sp.id <- SpecUPGMA$labels
    rmTip <- as.character(sp.id[! sp.id %in% comm[,1]]) #Remove species tips for species that are not present
    subTree <- drop.tip(functree, tip=rmTip)  # Remove those branches
    FDnull[k,i] <- sum(subTree$edge.length) # Sum up the branch length
  }}

FDnullid <- as.data.frame(cbind(FD,poIDs))
FDnullid$longitude <- 0
FDnullid$latitude <- 0
FDnullid$elevation <- 0
FDnullid$moIDs <- 0

for (k in 1:length(poIDs)){
  comm <- subset(adt.mo, poId == poIDs[k])
  FDnullid[k,102:105] <- comm[1,3:6]
}

FDnullid <- cbind(FDnullid[101:105],FDnullid[1:100])
saveRDS(FDnullid, file="FDnull.rds")

#### Calculate null expectation of functional duversity
FD <- readRDS(file="FD.rds")
FDnull <- readRDS(file="FDnull.rds")
FDnull.n <- cbind(FDnull[,1:5],FD[,1],FDnull[,6:ncol(FDnull)])

## Rank observed FD
ranks.n <- matrix(NA,length(poIDs),1)
FDnull.n$rank <- NA

for (i in 1:length(poIDs)){
  sel <- as.matrix(FDnull.n[,6:(ncol(FDnull.n)-1)])
  rank.sim <- apply(sel,1,rank)[,1]
  FDnull.n[i,ncol(FDnull.n)] <- rank.sim[1]
}

FDnull.n$pval <- FDnull.n$rank/101

## Quantify SES
FDnull.n$ses <- NA
for (i in 1:nrow(FDnull.n)){
  FDnull.n[i,ncol(FDnull.n)] <- ((as.numeric(FDnull.n[i,6]) - mean(as.numeric(FDnull.n[i,7:(ncol(FDnull.n)-3)])))/sd(as.numeric(FDnull.n[i,7:(ncol(FDnull.n)-3)])))
}


###################################################
###################################################STEP 3: Calculate contribution of individual traits
###################################################
#### Get species traits
## available at https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/13-1917.1)
sp.traits <- read.csv("data/Master_Funcdiv.csv", header=TRUE)
comm.traits <- matrix(NA, length(poIDs), 16)

for (k in 1:length(poIDs)) {
  comm <- subset(adt.mo, poId == poIDs[k])
  nspec <- nrow(comm)
  comm <- setDF(comm)
  
  sp.tr <- sp.traits[sp.traits[,2] %in% comm[,8],]
  comm <- comm[comm[,8] %in% sp.tr[,2],]
  sp.tr <- sp.tr[order(sp.tr[,2]),] 
  comm <- comm[order(comm[,8]),] 
  comm.tr <- cbind(comm, sp.tr[,4:20])
  
  inv <- mean(comm.tr$Inv)
  vert <- mean(comm.tr$Vert)
  scav <- mean(comm.tr$Scav)
  fruit <- mean(comm.tr$Fruits)
  nect <- mean(comm.tr$Nectar)
  seed <- mean(comm.tr$Seeds)
  plant <- mean(comm.tr$PlantOther)
  
  wbs <- mean(comm.tr$Height_watbelowsurface)
  was <- mean(comm.tr$Height_wataroundsurface)
  ground <- mean(comm.tr$Height_ground)
  under <- mean(comm.tr$Height_understory)
  midh <- mean(comm.tr$Height_midhigh)
  canop <- mean(comm.tr$Height_canopy)
  aerial <- mean(comm.tr$Height_aerial)
  
  noct <- mean(comm.tr$Nocturnal_cat)
  bmass <- mean(comm.tr$log10SpecGenFinalMass)
  
  comm.traits[k,] <- c(inv,vert,scav,fruit,nect,seed,plant,wbs,was,ground,under,midh,canop,aerial,noct,bmass)
}

saveRDS(comm.traits, file="Traits.rds")


###################################################
###################################################STEP 4: Calculate trait space overlap along elevational gradient
###################################################
FD <- readRDS("FD.rds")
ids <- FD[,c(2,5,6)]
ids$moIDs <- as.numeric(ids$moIDs)
ids$elevation <- as.numeric(ids$elevation)
ids$poIDs <- as.numeric(ids$poIDs)
ids <- ids[order(ids$moIDs, ids$elevation),]
moIDs <- as.numeric(unique(ids$moIDs))

ids <- ids %>%
  dplyr::mutate(elevcat = case_when(
    elevation < 500 ~ 1,
    elevation >= 500 & elevation < 1000 ~ 2,
    elevation >=1000 & elevation <1500 ~ 3,
    elevation >=1500 & elevation <2000 ~ 4,					
    elevation >=2000 & elevation <2500 ~ 5,
    elevation >=2500 & elevation <3000 ~ 6,
    elevation >=3000 & elevation <3500 ~ 7,
    elevation >=3500 & elevation <4000 ~ 8,
    elevation >=4000 & elevation <4500 ~ 9,
    elevation >=4500 & elevation <5000 ~ 10,
    elevation >=5000 & elevation <5500 ~ 11,
    elevation >=5500 & elevation <6000 ~ 12,
    elevation >=6000 & elevation <6500 ~ 13,					
    TRUE ~ 14)
  ) 

overlaps <- list() # List to hold overlap in trait space for each pair of points
for (i in 1:length(moIDs)){
  idsi <- ids[ids$moIDs == moIDs[i],] # Choose a mountain region
  elevs <- unique(idsi$elevcat)
  
  for (k in 2:length(elevs)){
    ids1 <- idsi[idsi$elevcat == elevs[k-1],]
    ids2 <- idsi[idsi$elevcat == elevs[k],]
    poids1 <- ids1[,1]
    poids2 <- ids2[,1]
    
    over.m <- matrix(NA, max(length(poids1),length(poids2)),11)
    
    if (length(poids1)<=1) {
      p1 <- rep(poids1,max(length(poids1),length(poids2)))
    } else {
      p1 <- sample(poids1, max(length(poids1),length(poids2)), replace = TRUE)
    }
    
    if (length(poids2)<=1) {
      p2 <- rep(poids2,max(length(poids1),length(poids2)))
    } else {
      p2 <- sample(poids2, max(length(poids1),length(poids2)), replace = TRUE)
    }
    
    for (j in 1:max(length(poids1),length(poids2))){
      vh1 <- readRDS(file=paste0("fdh_2D_",p1[j],".rds"))
      vh2 <- readRDS(file=paste0("fdh_2D_",p2[j],".rds"))
      hv_set <- hypervolume::hypervolume_set(vh1, vh2, check.memory=FALSE)
      stat <- hypervolume::hypervolume_overlap_statistics(hv_set)
      
      over.m[j,1] <- p1[j]
      over.m[j,2] <- p2[j]
      over.m[j,3] <- ids1[ids1[,1]==p1[j],2]
      over.m[j,4] <- ids2[ids2[,1]==p2[j],2]
      over.m[j,5] <- ids1[ids1[,1]==p1[j],3]
      over.m[j,6] <- ids1[ids1[,1]==p1[j],4]
      over.m[j,7] <- ids2[ids2[,1]==p2[j],4]
      over.m[j,8:11] <- stat
    }
    over.elev <- rbind(over.elev,over.m)
  }
  overlaps <- rbind(overlaps,over.elev)
}

saveRDS(overlaps, file="overlap.rds")

###################################################
###################################################STEP 5: INLA models (example INLA model for FD)
###################################################

#### INLA global models
moIDs <- unique(FD$moIDs)
n.loc <- length(unique(FD$moIDs))
FD$moID1 <- as.numeric(FD$moIDs)
FD$moID2 <- FD$moID1
FD$moID3 <- FD$moID1
FD$FDlog <- log(FD$FD) # Log-transform functional diversity

# Linear relationship
form <- FDlog ~ elevation + f(moID1, model = "iid") + f(moID2, elevation, copy="moID1")
result <- inla(form, data=FD, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE))
saveRDS(result, file="FD_elevation_global_1.rds")

# Mid-elevation peak
form <- logfd ~ elev + I(elev^2) + f(moID1, model = "iid") + f(moID2, elev, copy="moID1") + f(moID3, I(elev^2), copy="moID1")
result <- inla(form, data=fd, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE))
saveRDS(result, file="FD_elevation_global_2.rds")


### INLA models with latitude as an interaction term
# Linear relationship
form <- FDlog ~ elevation + latitude + elevation:latitude + f(moID1, model = "iid") + f(moID2, elevation, copy="moID1")
result <- inla(form, data=FD, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE))
saveRDS(result, file="FD_elevation-lat_global_1.rds")

# Mid-elevation peak
form <- FDlog ~ elevation + I(elevation^2) + latitude + elevation:latitude + I(elevation^2):latitude + f(moID1, model = "iid") + f(moID2, elevation, copy="moID1") + f(moID3, I(elevation^2), copy="moID1")
result <- inla(form, data=FD, family="Gaussian", control.compute = list(dic = TRUE, waic = TRUE))
saveRDS(result, file="FD_elevation-lat_global_2.rds")


