###########################################################################
# Program: exhaled.r                                                      #
# Ripped from: nothing, tabula rasa                                       #
# Purpose: identify sig. diffs in 'waves' btwn healthy, CF and asthma     #
#          patients                                                       #
# Prerequisite(s):                                                        #
# /home/karl/Dropbox/Cardioformatics/CTMM management/projects/            #
# Author: Karl Brand                                                      #
# Last update: 2013-05-17                                                 #
# Last edit: ## pick up here ##                                           #
###########################################################################

## preliminaries (set paths, load data, workspaces and required packages) 
root <- file.path("H:/MEP/Measurements/Final Data TUDelft")

##
library(limma)

## load data
use <- "raw"
## use <- "smoothed"
if (use == "raw") {
    dat.dir <- file.path(root, "raw")
    file.nams <- dir(dat.dir)
    ## i <- 1 ## debugging
    files <- list()
    for (i in 1:length(file.nams)) {
        files[[i]] <- read.table(file = file.path(dat.dir,
                                     file.nams[i]),
                                 sep = "", skip = 1, strip = TRUE,
                                 colClasses = "numeric")
        names(files)[i] <- substr(file.nams[i], 1, 20)
        ## make new column names
        col.nams <- matrix(
            read.table(file = file.path(dat.dir,
                           file.nams[i]),
                       sep = "\t", header = FALSE, nrows = 1, strip = TRUE)
        )
        col.nams. <- sub("/", ".", x = col.nams)
        grp       <- substr(file.nams[i], 1, 7)
        grp       <- sub("Group_0", "healthy", x = grp)
        grp       <- sub("Group_1", "asthma", x = grp)
        grp       <- sub("Group_2", "CF", x = grp)
        names(files[[i]]) <- paste(grp, col.nams., sep = ".")
    }
    ## inspect:
    ## head(files[[1]])
    ## dim(files[[1]])
} else if(use == "smoothed") {
    dat.dir <- file.path(root, "smoothed")
    file.nams <- dir(dat.dir)
    file.nams
    ## i <- 1 ## debugging
    files <- list()
    for (i in 1:length(file.nams)) {
        files[[i]] <- read.table(file = file.path(dat.dir,
                                     file.nams[i]),
                                 sep = "", skip = 1, strip = TRUE)
        names(files)[i] <- substr(file.nams[i], 1, 25)
        ## make new column names
        col.nams <- matrix(
            read.table(file = file.path(dat.dir,
                           file.nams[i]),
                       sep = "\t", header = FALSE, nrows = 1, strip = TRUE)
        )
        col.nams. <- sub("/", ".", x = col.nams)
        names(files[[i]]) <- paste("healthy", col.nams., sep = ".")
    }
    ## inspect:
    ## names(files)
    ## head(files[[1]])
    ## dim(files[[1]])
}
str(files)

## plot density
## lapply(files, function(x) {
##     for (i in 1:length(x)) {
##         plot(density(x[, i]))
##     }
## }
##        )

head(files[[1]])
sum(unlist(lapply(files, ncol)))
names(files)
unlist(lapply(files, names))
str(files)


## bind together as one 'expression set'
ese <- do.call(cbind, files)
## make nicer column names
names(ese) <- unlist(lapply(files, names))
## convert to matrix
eset <- as.matrix(ese)
class(eset)
## inspect data minima and maxima
min(eset)
max(eset)

## boxplot each sample,  raw data
x11(width = 11.69, height = 8.27)
## pdf(file.path(root, "results", "plots", "PDFs", paste(use, "boxplot.alldata.pdf")),
##     width = 11.69, height = 8.27)
## png(file.path(root, "results", "plots", "PNGs",  paste(use, "boxplot.alldata.png")),
##     width=11.69, height=8.27, units="in", pointsize=12,
##     bg="white",  res=600)
boxplot(ese) ## to do: make sample names vertical
dev.off()

## density plot each sample, raw data
x11(xpos= 2000)
## pdf(file.path(root, "results", "plots", "PDFs", paste(use, "desnity.alldata.pdf")),
##     width = 8.27, height = 8.27)
plot(x = NULL, type = "n",
     ## xlim = c(min(eset), max(eset)), ylim = c(0, 40),
     xlim = c(-0.1, 0.15), ylim = c(0, 40),
     xlab = NA, ylab = "Density")
for (i in 1:ncol(eset)) {
    ## i <- 1
    lines(density(eset[, i]), col = "grey")
}
## overlay density of all data
lines(density(eset), col = "black")
dev.off()

## row.names(ese.l) <- grp0[, 1]
head(ese)
dim(ese.l)

## ## based on boxplot inspection:
## ## Will need to exclude "asthma.005.02" "asthma.019.02" "asthma.021.02"
## ## since they're medians are much higher than the others
## x11(width = 11.27, height = 8.69)
## ## pdf(file.path(root, "results", "plots", "PDFs", paste(use, "boxplot.rem.pdf")),
## ##     width = 11.69, height = 8.27)
## ## png(file.path(root, "results", "plots", "PNGs", paste(use, "boxplot.rem.png")),
## ##     width=11.69, height=8.27, units="in", pointsize=12,
## ##     bg="white",  res=600)
## boxplot(ese[-c(38, 52, 54)])
## dev.off()
## names(ese[c(38, 52, 54)])

## ## remove them
## eset <- ese #[-c(38, 52, 54)]


## ## inspect data structure
## x11(xpos = 1000)
## plot(density(
##     log2(as.matrix(grp0[, -1]))
##     ))
## dev.off()
## x11(xpos = 1000)
## plot(density(
##     log2(as.matrix(grp1[, -1]))
##     ))
## dev.off()
## x11(xpos = 1000)
## plot(density(
##     log2(as.matrix(grp2[, -1]))
##     ))
## dev.off()

## generally ok - roughly gaussian

## load libraries
library("limma")

## set up factors 
File <- dimnames(eset)[[2]]
File.sp <- strsplit(File, split="[.]")
Treatment <- unlist(lapply(File.sp, "[", 1))
Patient   <- unlist(lapply(File.sp, "[", 2))
Time <- unlist(lapply(File.sp, "[", 3))
Target <- paste(Treatment, Time, sep=".")
## make targets 'file'
targets <- data.frame("FileName" = File,    "Treatment" = Treatment,
                      "Patient"  = Patient, "Time" = Time,
                      "Target"   = Target)
targets

## specify levels for factors
Treatment <- factor(Treatment, levels=c("healthy", "asthma", "CF"))
Time <- factor(Time)

## inspect correlations within patients
## corfit0 <- duplicateCorrelation(eset[, 1:36], ndups = 4,
##                                 block = rep(1:4, times = 9))
## corfit1 <- duplicateCorrelation(eset[, 37:54], ndups = 2,
##                                 block = rep(1:2, times = 9))
## corfit2 <- duplicateCorrelation(eset[, 55:64], ndups = 2,
##                                 block = rep(1:2, times = 5))
## corfit0$consensus ## [1] 0.04446755
## corfit1$consensus ## [1] -0.01083684
## corfit2$consensus ## [1] -0.07068767
## all very small so (i assume) can be ignored (?)

## Analysis 1 - all data at times 1 and 2
## first remove all the time points 3 and 4
rem <- which(targets$Time == "03" | targets$Time == "04")
## now exclude "asthma.005.02" "asthma.019.02" "asthma.021.02"
## exc <- which(targets$FileName == "asthma.005.02" |
##              targets$FileName == "asthma.019.02" |
##              targets$FileName == "asthma.021.02")

eset.1 <- eset[, -rem]
dim(eset.1)
targs.1 <- targets[-rem, ]
dim(targs.1)
## inspect to check removal worked
ncol(eset.1) == nrow(targs.1)
cbind(colnames(eset.1),  targs.1$Target)

# # # # # # #  set up design using Buffy's parameterization # # # # # # # #
## design.1 <- model.matrix(~ Treatment + Antagomir)
Target.fac.1 <- factor(targs.1$Target,
                       levels = c("healthy.01", "healthy.02",
                           "asthma.01", "asthma.02",
                           "CF.01", "CF.02"))
design.1 <- model.matrix(~0 + Target.fac.1)#, data=targets)
as.data.frame( design.1)

## give nicer column names
colnames(design.1) <- c("healthy.01", "healthy.02",
                           "asthma.01", "asthma.02",
                           "CF.01", "CF.02")
design.1

## fit the model and make contrasts
fit.1 <- lmFit(eset.1, design.1)
cont.mat.1 <- makeContrasts(
    "healthy.vs.asthma"        = (healthy.01 + healthy.02) -
                                 (asthma.01 + asthma.02),
    "healthy.01.vs.asthma.01"  = healthy.01 - asthma.01,     ## keep
    "healthy.02.vs.asthma.02"  = healthy.02 - asthma.02,
    "healthy.vs.CF"            = (healthy.01 + healthy.02) -
                                 (CF.01 + CF.02),
    "healthy.01.vs.CF.01"      = healthy.01 - CF.01,         ## keep
    "healthy.02.vs.CF.02"      = healthy.02 - CF.02,
    "asthma.vs.CF"             = (asthma.01 + asthma.02) -
                                 (CF.01 + CF.02),
    "healthy.01.vs.healthy.02" = healthy.01 - healthy.02,    ## keep
    "asthma.01.vs.asthma.02"   = asthma.01 - asthma.02,      ## keep
    "CF.01.vs.CF.02"           = CF.01 - CF.02,              ## keep
    "t01.vs.t02"               = (healthy.01 + asthma.01 + CF.01) -
                                 (healthy.02 + asthma.02 + CF.02),
    levels = design.1
    )
cont.mat.1

fit.1.2 <- contrasts.fit(fit.1, cont.mat.1)
fit.1.3 <- eBayes(fit.1.2)

## ## healthy.vs.asthma ## ##
tt.healthy.vs.asthma  <- topTable(fit.1.3, coef="healthy.vs.asthma",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.healthy.vs.asthma)
tt.healthy.vs.asthma[order(tt.healthy.vs.asthma$P.Value), ][1:20, ]
sum(tt.healthy.vs.asthma$adj.P.Val < 0.05) # 1255 # NA # 546
hist(tt.healthy.vs.asthma$adj.P.Val)
head(tt.healthy.vs.asthma)
## write.csv(tt.healthy.vs.asthma,
##           file = file.path(root, "results", "tables",
##                "healthy.vs.asthma.alldata.csv"), row.names = FALSE)

## "healthy.1.vs.asthma.1" = healthy.1 - asthma.1,
tt.healthy.1.vs.asthma.1  <- topTable(fit.1.3, coef="healthy.01.vs.asthma.01",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.healthy.1.vs.asthma.1)
tt.healthy.1.vs.asthma.1[order(tt.healthy.1.vs.asthma.1$P.Value), ][1:20, ]
sum(tt.healthy.1.vs.asthma.1$adj.P.Val < 0.05) # 68 # NA # 56
hist(tt.healthy.1.vs.asthma.1$adj.P.Val)
## write.csv(tt.healthy.1.vs.asthma.1,
##           file = file.path(root, "results", "tables",
##                "healthy.1.vs.asthma.1.alldata.csv"), row.names = FALSE)

## "healthy.2.vs.asthma.2" = healthy.2 - asthma.2,
tt.healthy.2.vs.asthma.2  <- topTable(fit.1.3, coef="healthy.02.vs.asthma.02",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.healthy.2.vs.asthma.2)
tt.healthy.2.vs.asthma.2[order(tt.healthy.2.vs.asthma.2$P.Value), ][1:20, ]
sum(tt.healthy.2.vs.asthma.2$adj.P.Val < 0.05) # 319 # NA # 14
hist(tt.healthy.2.vs.asthma.2$adj.P.Val)
## write.csv(tt.healthy.2.vs.asthma.2,
##           file = file.path(root, "results", "tables",
##                "healthy.2.vs.asthma.2.alldata.csv"), row.names = FALSE)

## ## healthy.vs.CF ## ##
tt.healthy.vs.CF  <- topTable(fit.1.3, coef="healthy.vs.CF",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.healthy.vs.CF)
tt.healthy.vs.CF[order(tt.healthy.vs.CF$P.Value), ][1:20, ]
sum(tt.healthy.vs.CF$adj.P.Val < 0.05) # 36 # NA # 353
hist(tt.healthy.vs.CF$adj.P.Val)
## write.csv(tt.healthy.vs.CF,
##           file = file.path(root, "results", "tables",
##                "healthy.vs.CF.alldata.csv"), row.names = FALSE)

## ## healthy.1.vs.CF.1 ## ##
tt.healthy.1.vs.CF.1  <- topTable(fit.1.3, coef="healthy.01.vs.CF.01",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.healthy.1.vs.CF.1)
tt.healthy.1.vs.CF.1[order(tt.healthy.1.vs.CF.1$P.Value), ][1:20, ]
sum(tt.healthy.1.vs.CF.1$adj.P.Val < 0.05) # [1] 0 #    # 59
sum(tt.healthy.1.vs.CF.1$adj.P.Val < 0.10) # [1] 0
hist(tt.healthy.1.vs.CF.1$adj.P.Val)
## write.csv(tt.healthy.vs.CF,
##           file = file.path(root, "results", "tables",
##                "healthy.1.vs.CF.1.alldata.csv"), row.names = FALSE)

## ## healthy.2.vs.CF.2 ## ##
tt.healthy.2.vs.CF.2  <- topTable(fit.1.3, coef="healthy.02.vs.CF.02",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.healthy.2.vs.CF.2)
tt.healthy.2.vs.CF.2[order(tt.healthy.2.vs.CF.2$P.Value), ][1:20, ]
sum(tt.healthy.2.vs.CF.2$adj.P.Val < 0.05) # [1] 0  #    # 5
sum(tt.healthy.2.vs.CF.2$adj.P.Val < 0.10) # [1] 18 #    # 18
hist(tt.healthy.2.vs.CF.2$adj.P.Val)
## write.csv(tt.healthy.2.vs.CF.2,
##           file = file.path(root, "results", "tables",
##                "healthy.2.vs.CF.2.alldata.csv"), row.names = FALSE)

## ## "asthma.vs.CF"      = (asthma.1 + asthma.2) - (CF.1 + CF.2),
tt.asthma.vs.CF  <- topTable(fit.1.3, coef="asthma.vs.CF",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.asthma.vs.CF)
tt.asthma.vs.CF[order(tt.asthma.vs.CF$P.Value), ][1:20, ]
sum(tt.asthma.vs.CF$adj.P.Val < 0.05) # [1] 0
sum(tt.asthma.vs.CF$adj.P.Val < 0.10) # [1] 0
hist(tt.asthma.vs.CF$adj.P.Val)
hist(tt.asthma.vs.CF$P.Value)

## ## "healthy.1.vs.healthy.2" = healthy.1 - healthy.2,
tt.healthy.1.vs.healthy.2  <- topTable(fit.1.3, coef="healthy.01.vs.healthy.02",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.healthy.1.vs.healthy.2)
tt.healthy.1.vs.healthy.2[order(tt.healthy.1.vs.healthy.2$P.Value), ][1:20, ]
sum(tt.healthy.1.vs.healthy.2$adj.P.Val < 0.05) # [1] 0 #     # 2
sum(tt.healthy.1.vs.healthy.2$adj.P.Val < 0.10) # [1] 0 #     # 3
hist(tt.healthy.1.vs.healthy.2$adj.P.Val)

## ## "asthma.1.vs.asthma.2"   = asthma.1 - asthma.2,
tt.asthma.1.vs.asthma.2  <- topTable(fit.1.3, coef="asthma.01.vs.asthma.02",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.asthma.1.vs.asthma.2)
tt.asthma.1.vs.asthma.2[order(tt.asthma.1.vs.asthma.2$P.Value), ][1:20, ]
sum(tt.asthma.1.vs.asthma.2$adj.P.Val < 0.05) # [1] 0
sum(tt.asthma.1.vs.asthma.2$adj.P.Val < 0.10) # [1] 0
hist(tt.asthma.1.vs.asthma.2$adj.P.Val)

## ## "CF.1.vs.CF.2"      = CF.1 - CF.2,
tt.CF.1.vs.CF.2  <- topTable(fit.1.3, coef="CF.01.vs.CF.02",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.CF.1.vs.CF.2)
tt.CF.1.vs.CF.2[order(tt.CF.1.vs.CF.2$P.Value), ][1:20, ]
sum(tt.CF.1.vs.CF.2$adj.P.Val < 0.05) # [1] 0 #     # 1
sum(tt.CF.1.vs.CF.2$adj.P.Val < 0.10) # [1] 0 #     # 1
hist(tt.CF.1.vs.CF.2$adj.P.Val)

## ## "t1.vs.t2"   = (healthy.1 + asthma.1 + CF.1) - (healthy.2 + asthma.2 + CF.2),
tt.t1.vs.t2  <- topTable(fit.1.3, coef="t01.vs.t02",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.t1.vs.t2)
tt.t1.vs.t2[order(tt.t1.vs.t2$P.Value), ][1:20, ]
sum(tt.t1.vs.t2$adj.P.Val < 0.05) # [1] 0
sum(tt.t1.vs.t2$adj.P.Val < 0.10) # [1] 0
hist(tt.t1.vs.t2$adj.P.Val)

## PCA plot ripped from "pca_h.clutersing_QC.r"
## can only get subsequent NA replacement working on data frame
rep.nas <- function(x, i) {
  grp <- x[, i]
  grp.df <- data.frame(x[, i])
## row means
  rmens <- rowMeans(grp.df, na.rm=TRUE)
  grp.df[is.na(grp)] <-
    rmens[c(unlist(sapply(grp.df, function(y) which(is.na(y)) )))]
  return(grp.df)
}

head(eset.1)

## select groups of interest
eset.1.samp <- dimnames(eset.1)[[2]]
eset.1.samp
healthy.T1 <- which(targs.1$Target == "healthy.01")
healthy.T2 <- which(targs.1$Target == "healthy.02")
health.samp <- c(healthy.T1, healthy.T2)

asthma.T1 <- which(targs.1$Target == "asthma.01")
asthma.T2 <- which(targs.1$Target == "asthma.02")
asthma.samp <- c(asthma.T1, asthma.T2)

CF.T1 <- which(targs.1$Target == "CF.01")
CF.T2 <- which(targs.1$Target == "CF.02")
CF.samp <- c(CF.T1, CF.T2)


## head(eset.1)
## plot(density(na.omit(eset.1[, 13])))

## PCA 3 main groups ######################################################
grp1 <- rep.nas(eset.1, i = health.samp)
## check it worked
tail(grp1[, 1:4])
tail(eset.1[, 1:4])                       # works
grp2 <- rep.nas(eset.1, i = asthma.samp)
grp3 <- rep.nas(eset.1, i = CF.samp)
## bind all gropus back together
eset.1.nas <-na.omit(cbind(grp1, grp2, grp3))
nrow(eset.1.nas)
pca <- prcomp(x=t(eset.1.nas))

## pdf(file=file.path(root, "results", "plots", "PDFs",
##       paste("van_mastrigt.PCA.plot", use, ".pdf")))
x11(xpos=2000)
plot(pca$x[,1:2], type="n")
points(pca$x[,1:2], bg=c(rep("red", length(health.samp)),
                         rep("blue", length(asthma.samp)),
                         rep("green", length(CF.samp))),
                        col=c(rep("red", length(health.samp)),
                              rep("blue", length(asthma.samp)),
                              rep("green", length(CF.samp))),
                            pch=21)
## add a legend
legend("topleft",
       legend=c("healthy", "asthma", "CF"),
       fill=c("red", "blue", "green"))
## add a title
title(paste("PCA for", use, "data", sep=" "))
#dev.off()

## PCA 6 groups ######################################################
grp1 <- rep.nas(eset.1, i = healthy.T1)
## check it worked
tail(grp1[, 1:4])
tail(eset.1[, 1:4])                       # works
grp2 <- rep.nas(eset.1, i = healthy.T1)
grp3 <- rep.nas(eset.1, i = asthma.T1)
grp4 <- rep.nas(eset.1, i = asthma.T2)
grp5 <- rep.nas(eset.1, i = CF.T1)
grp6 <- rep.nas(eset.1, i = CF.T2)
## bind all gropus back together
eset.1.nas <-na.omit(cbind(grp1, grp2, grp3, gr4, grp5, grp6,))
nrow(eset.1.nas)
pca <- prcomp(x=t(eset.1.nas))

## pdf(file=file.path(root, "results", "plots", "PDFs",
##       paste("van_mastrigt.PCA.plot", use, ".pdf")))
x11(xpos=2000)
plot(pca$x[,1:2], type="n")
points(pca$x[,1:2], bg=c(rep("red", length(healthy.T1)),
                         rep("red4", length(healthy.T2)),
                         rep("blue", length(asthma.T1)),
                         rep("lightblue", length(asthma.T2)),
                         rep("green", length(CF.T1)),
                         rep("lightgreen", length(CF.T1))),
                        
                        col=c(rep("red", length(healthy.T1)),
                            rep("red4", length(healthy.T2)),
                            rep("blue", length(asthma.T1)),
                            rep("lightblue", length(asthma.T2)),
                            rep("green", length(CF.T1)),
                            rep("lightgreen", length(CF.T1))),
                            pch=21)
## add a legend
legend("topleft",
       legend=c("healthy.T1", "healthy.T2",
                "asthma.T1", "asthma.T2",
                "CF.T1", "CF.T2"),
       fill=c("red", "red4",
              "blue", "lightblue",
              "green", "lightgreen"))
## add a title
title(paste("PCA for", use, "data", sep=" "))
dev.off()



gc()
