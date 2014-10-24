###########################################################################
# Program: exhaled.r                                                      #
# Ripped from: nothing, tabula rasa                                       #
# Purpose: identify sig. diffs in 'waves' btwn healthy, CF and asthma     #
#          patients                                                       #
# Prerequisite(s):                                                        #
# /home/karl/Dropbox/Cardioformatics/CTMM management/projects/            #
# Author: Karl Brand                                                      #
# Last update: 2014-05-05                                                 #
# Last edit: ## pick up here ##                                           #
###########################################################################

## Brand Still to do:
## -reproducibility btwn time points on the updated dataset (done)
## -correlations btwn time points (done)
## -label the PCA with smaple names (pending)
## -forward normalisation (quantile) description (pending)

## Changes by Olav 
## - col.nams: removed matrix() and added ColClasses="character" argument to read.table. This fixes the names of the samples.
## - Changed paths, and set use="raw"


## starting again?
rm(list = ls())

## preliminaries (set paths, load data, workspaces and required® packages) 
root <- file.path("L:","IST","OP","scratch","Olav Grouwstra","Measurements")
datadate <- "26-8-2014"
dataof <- paste("data of",datadate,sep=" ")
zeros <- "no zeros" #Empty if zeros are at end of the vectors
nozeros <- paste("data",zeros,sep=" ")
###########################################################################
#                                 load data                               #
###########################################################################

## choose to use raw or smoothed data
use <- "raw"
## use <- "smoothed"

## choose a scaling method
## scal.meth <- "MAD"
scal.meth <- "quantile"
## scal.meth <- "none"


## import
if (use == "raw") {
    dat.dir <- file.path(root, "CO2H2O calibrated data",dataof,nozeros)
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
        col.nams <- read.table(file = file.path(dat.dir,
                           file.nams[i]),
                       sep = "\t", header = FALSE, nrows = 1, strip = TRUE,
                       colClasses="character")
      
        col.nams. <- gsub("/", ".", x = col.nams)
        grp       <- substr(file.nams[i], 1, 2)
        grp       <- sub("G0", "healthy", x = grp)
        grp       <- sub("G1", "asthma", x = grp)
        grp       <- sub("G2", "CF", x = grp)
        names(files[[i]]) <- paste(grp, col.nams., sep = ".")
    }
    ## inspect:
    ## head(files[[1]])
    ## dim(files[[1]])
} else if(use == "smoothed") {
    dat.dir <- file.path(root, "CO2H2O calibrated data",dataof,nozeros)
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
        grp       <- sub("G0", "healthy", x = grp)
        grp       <- sub("G1", "asthma", x = grp)
        grp       <- sub("G2", "CF", x = grp)
        names(files[[i]]) <- paste(grp, col.nams., sep = ".")
    }
    ## inspect:
    ## names(files)
    ## head(files[[1]])
    ## dim(files[[1]])
}

## inspect imported data
str(files)
head(files[[1]])
sum(unlist(lapply(files, ncol)))
names(files)
unlist(lapply(files, names))

###########################################################################
#      bind together as one 'expression set' and visualise raw data       #
###########################################################################
ese <- do.call(cbind, files)
head(ese)
## make nicer column names
names(ese) <- unlist(lapply(files, names))
## convert to matrix
eset <- as.matrix(ese)
str(eset)
class(eset)
## inspect data minima and maxima
min(eset)
max(eset)

# ###########################################################################
# #               density and boxplots of unscaled/normed data              #
# ###########################################################################
# 
# boxplots <- function(ese, eset, pdf = FALSE, png = FALSE) {
#     ## boxplot each sample, raw data # # # # # # # # # # # # # # # # # # #  
#     if (pdf == TRUE) {
#         pdf(file.path(root, "results", "plots", "PDFs",
#                       paste(use, "boxplot.alldata.pdf")),
#             width = 11.69, height = 8.27)
#         boxplot(ese,
#                 main = paste("Untransformed, unscaled", use, "data"))
#         dev.off()
#     } else if (png == TRUE) {
#         png(file.path(root, "results", "plots", "PNGs",
#                       paste(use, "boxplot.alldata.png")),
#             width=11.69, height=8.27, units="in", pointsize=12,
#             bg="white",  res=600)
#         boxplot(ese,
#                 main = paste("Untransformed, unscaled", use, "data"))
#         dev.off()
#     } else {
#         x11(width = 11.69, height = 8.27)
#         boxplot(ese,
#                 main = paste("Untransformed, unscaled", use, "data"))
#     }
#     
#     if (use == "raw") {
#         ## density plot each sample, raw data # # # # # # # # # # # # # # #
#         if (pdf == TRUE) {
#             pdf(file.path(root, "results", "plots", "PDFs",
#                           paste(use, "density.alldata.pdf")),
#                 width = 8.27, height = 8.27)
#             plot(x = NULL, type = "n",
#                  ## xlim = c(min(eset), max(eset)), ylim = c(0, 40),
#                  xlim = c(-0.1, 0.15), ylim = c(0, 40),
#                  xlab = NA, ylab = "Density")
#             for (i in 1:ncol(eset)) {
#                 ## i <- 1
#                 lines(density(eset[, i]), col = "grey")
#             }
#             ## overlay density of all data
#             lines(density(eset), col = "black")
#             dev.off()
#             ## generally ok - gaussian, if a little right skewed
#         } else {
#             x11(xpos = 2000)
#             plot(x = NULL, type = "n",
#                  ## xlim = c(min(eset), max(eset)), ylim = c(0, 40),
#                  xlim = c(-0.1, 0.15), ylim = c(0, 40),
#                  xlab = NA, ylab = "Density")
#             for (i in 1:ncol(eset)) {
#                 ## i <- 1
#                 lines(density(eset[, i]), col = "grey")
#             }
#             ## overlay density of all data
#             lines(density(eset), col = "black")
#         }
#         
#     } else if (use == "smoothed") {
#         if (pdf == TRUE) {
#             ## density plot each sample, smoothed data  # # # # # # # # # #
#             pdf(file.path(root, "results", "plots", "PDFs",
#                           paste(use, "density.alldata.pdf")),
#                 width = 8.27, height = 8.27)
#             plot(x = NULL, type = "n",
#                  ## xlim = c(min(eset), max(eset)), ylim = c(0, 40),
#                  ## xlim = c(-15, 0), ylim = c(0, 0.50),
#                  xlim = c(-0.02, 0.08), ylim = c(0, 100),
#                  xlab = NA, ylab = "Density",
#                  main = paste("unscaled", use, "data"))
#             for (i in 1:ncol(eset)) {
#                 ## i <- 1
#                 lines(density(eset[, i]), col = "grey")
#             }
#             ## overlay density of all data
#             lines(density(eset), col = "black")
#             dev.off()
#         } else {
#             x11(xpos = 2000)
#             plot(x = NULL, type = "n",
#                  ## xlim = c(min(eset), max(eset)), ylim = c(0, 40),
#                  ## xlim = c(-15, 0), ylim = c(0, 0.50),
#                  xlim = c(-0.02, 0.08), ylim = c(0, 100),
#                  xlab = NA, ylab = "Density",
#                  main = paste("unscaled", use, "data"))
#             for (i in 1:ncol(eset)) {
#                 ## i <- 1
#                 lines(density(eset[, i]), col = "grey")
#             }
#             ## overlay density of all data
#             lines(density(eset), col = "black")
#             ## strongly right skewed
#         }
#         if (pdf == TRUE) {
#             ## density plot each sample, log2 transformed smoothed data # # 
#             pdf(file.path(root, "results", "plots", "PDFs",
#                           paste("log2", use, "density.alldata.pdf")),
#                 width = 8.27, height = 8.27)
#             plot(x = NULL, type = "n",
#                  ## xlim = c(min(eset), max(eset)), ylim = c(0, 40),
#                  xlim = c(-15, 0), ylim = c(0, 0.50),
#                  ## xlim = c(-0.02, 0.08), ylim = c(0, 100),
#                  xlab = NA, ylab = "Density",
#                  main = paste("log2 transformed unscaled", use, "data"))
#             for (i in 1:ncol(eset)) {
#                 ## i <- 1
#                 lines(density(log2(eset[, i])), col = "grey")
#             }
#             ## overlay density of all data
#             lines(density(log2(eset)), col = "black")
#             dev.off()
#         } else {
#             x11(xpos = 2000)
#             plot(x = NULL, type = "n",
#                  ## xlim = c(min(eset), max(eset)), ylim = c(0, 40),
#                  xlim = c(-15, 0), ylim = c(0, 0.50),
#                  ## xlim = c(-0.02, 0.08), ylim = c(0, 100),
#                  xlab = NA, ylab = "Density",
#                  main = paste("log2 transformed unscaled", use, "data"))
#             for (i in 1:ncol(eset)) {
#                 ## i <- 1
#                 lines(density(log2(eset[, i])), col = "grey")
#             }
#              ## overlay density of all data
#             lines(density(log2(eset)), col = "black")
#         }
#     }
# }
# 
# ## instpect
# boxplots(ese = ese, eset = eset)
# ## moderately left skewed, bi-modalish, hence need to transform
# dev.off()
# dev.off()
# dev.off()
# 
# ## save:
# ## boxplots(ese = ese, eset = eset, pdf = TRUE)
# 
# ## keep untransformed data for testing purposes:
# keep.unt <- eset
# ## eset <- keep.unt
# 
# if (use == "smoothed") {
#     eset <- log2(eset)
# }

###########################################################################
#                         scaling/normalisation                           #
###########################################################################
## choose a scaling method
## scal.meth <- "MAD"
## scal.meth <- "quantile"
## scal.meth <- "none"

###############################################
# MAD scaling (modified from Wijnen analysis) #
###############################################
if (scal.meth == "MAD") {
    ## function which subtracts the sample median, divides by the absolute
    ## sample median then adds back the gobal median to all results
    eset[1:5, 1:5]
    dim(eset)
    mad.scale <- function(x, method=match.fun(median), na.rm = TRUE, ...) {
        ## x <- eset           # used for debugging
        ## scl <- median       # used for debugging
        ## i <- 1              # used for debugging
        scl <- method
        ## global mean/median 
        glob.me <- scl(x, na.rm=na.rm, ...)
        print(glob.me)
        y <- apply(x, 2, function(x) {
            ## calculate residuals by subtracting sample median from
            ## each feature:
            med.res <- x - scl(x, na.rm=na.rm, ...)
            ## now get median of the residuals,
            ## i.e.,  median absolute deviation (MAD)
            mad <- scl(abs(med.res), na.rm=na.rm, ...)
            ## median distance away from the median (MAD) ...as
            ## opposed to the square-root of the average distance form
            ## the mean normalise the variance of each plate
            med.res.v <- med.res / mad
            ## med.res.v <- med.res/(sqrt(var(med.res)))
            ## add global median to all
            med.res.gl <- med.res.v + glob.me
            return(med.res.gl)
        }
                   )
        return(y)
    }

    eset.s <- mad.scale(eset, median)
    dim(eset.s)
    ## plot(density(eset.s[, 1]))

    ## density plot each sample, MAD scaled data # # # # # # # # # # # # # 
    x11(xpos= 2000)
    ## pdf(file.path(root, "results", "plots", "PDFs",
    ##               paste(use, "density.alldata.pdf")),
    ##     width = 8.27, height = 8.27)
    plot(x = NULL, type = "n",
         ## xlim = c(min(eset.s), max(eset.s)), ylim = c(0, 40),
         xlim = c(-7, 10), ylim = c(0, 0.365),
         xlab = NA, ylab = "Density",
         main = paste(scal.meth, "normalised", use, "data"))
    for (i in 1:ncol(eset.s)) {
        ## i <- 1
        lines(density(eset.s[, i]), col = "grey")
    }
    ## overlay density of all data
    lines(density(eset.s), col = "black")
    ## dev.off()

#     ## boxplot each sample, AMD scaled data# # # # # # # # # # # # # # # # 
#     x11(width = 19, height = 10)
#     ## pdf(file.path(root, "results", "plots", "PDFs",
#     ##               paste(use, "boxplot.alldata.pdf")),
#     ##     width = 11.69, height = 8.27)
#     ## png(file.path(root, "results", "plots", "PNGs",
#     ##               paste(use, "boxplot.alldata.png")),
#     ##     width=11.69, height=8.27, units="in", pointsize=12,
#     ##     bg="white",  res=600)
#     boxplot(eset.s,
#             main = paste(scal.meth, "normalised", use, "data"))
    ## to do: make sample names vertical
    ## dev.off()

###############################################
# quantiles normalisation                     #
###############################################
} else if (scal.meth == "quantile") {
    library(preprocessCore)
    any(is.na(eset))
    min(eset)
    if (use == "raw") {
        eset.s <- normalize.quantiles(x = eset, copy = TRUE)
    } else if (use == "smoothed") {
        eset.s <- normalize.quantiles.robust(
            x = (eset + abs(min(eset))), copy = TRUE, weights = NULL,
            remove.extreme = "both", ## c("variance","mean","both","none"),
            n.remove    = 2, use.median  = FALSE, use.log2 = FALSE) #TRUE)
    }
    ## quickly inspect
    eset.s[1:5, 1:5]
    eset[1:5, 1:5]
    class(eset.s)
    class(eset)

    ## add back column names to quant. noramlised matrix
    dimnames(eset.s) <- dimnames(eset)

    ## quick inspection of the density
    ## plot(density(eset.s[, 1],  na.rm = TRUE))

    ## density plot each sample, quantile normed data # # # # # # # # # # #
    x11(xpos= 2000)
    ## pdf(file.path(root, "results", "plots", "PDFs",
    ##               paste(scal.meth, "normed", use, "density.alldata.pdf",
    ##                     sep = ".")),
    ##     width = 8.27, height = 8.27)
#     if (use == "raw") {
#         plot(x = NULL, type = "n",
#              xlim = c(-0.07, 0.12), ylim = c(0, 28),
#              xlab = NA, ylab = "Density",
#              main = paste(scal.meth, "normalised", use, "data"))
#     } else if (use == "smoothed") {
#         plot(x = NULL, type = "n",
#              ## xlim = c((min(eset.s)-0.01), max(eset.s)), ylim = c(0, 90),
#              xlim = c(min(eset.s), max(eset.s)), ylim = c(0, 0.40),
#              xlab = NA, ylab = "Density",
#              main = paste(scal.meth, "normalised", "log2 transformed",
#                  use, "data"))
#     }
#     for (i in 1:ncol(eset.s)) {
#         ## i <- 1
#         lines(density(eset.s[, i]), col = "grey")
#     }
    ## overlay density of all data
#     lines(density(eset.s), col = "black")
    ## dev.off()

#     ## boxplot each sample, quantile scaled data # # # # # # # # # # # # # 
#     x11(width = 19, height = 10)
#     ## pdf(file.path(root, "results", "plots", "PDFs",
#     ##               paste(scal.meth, "normed", use, "boxplot.alldata.pdf",
#     ##                     sep = ".")),
#     ##     width = 11.69, height = 8.27)
#     ## png(file.path(root, "results", "plots", "PNGs",
#     ##               paste(use, "boxplot.alldata.png")),
#     ##     width=11.69, height=8.27, units="in", pointsize=12,
#     ##     bg="white",  res=600)
#     boxplot(eset.s,
#             main = paste(scal.meth, "normalised", "log2 transformed",
#                 use, "data"))
#     ## to do: make sample names vertical
#     ## dev.off()
}
# dev.off()
# dev.off()

###########################################################################
#                        feature selection using limma                    #
###########################################################################
## Use type scaled or unscaled data

eset[1:5, 1:5]
eset.s[1:5, 1:5]
if (scal.meth != "none") {
eset <- eset.s
}
all(eset == eset.s)

## export the entire data set, as used for limma, here:
dim(eset)
## write.table(eset, file = file.path(root, "results", "tables",
##               paste(scal.meth, "normalised", "log2 transformed", use,
##                 "all.data.b.txt", sep=".")), sep="\t", row.names = FALSE)

## load libraries
library(limma)

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

## Analysis 1 - all data at times 1 and 2 !deprecated 2104.07.21!
## first remove all the time points 3 and 4

## analyse all time points via:
timepoints <- "all"
## timepoints <- "only1and2"

if (timepoints == "only1and2") {
    rem <- which(targets$Time == "03" | targets$Time == "04")
    eset.1 <- eset[, -rem]
    targs.1 <- targets[-rem, ]
    level.labs <- c("healthy.01", "healthy.02",
                    "asthma.01", "asthma.02",
                    "CF.01", "CF.02")
    cat(paste("\nanalysing", ncol(eset.1), "samples\n"))
} else if (timepoints == "all") {
    targs.1 <- targets
    eset.1 <- eset
    level.labs <- c("healthy.01", "healthy.02",
                    "healthy.03", "healthy.04",
                    "asthma.01", "asthma.02",
                    "CF.01", "CF.02")
    cat(paste("\nanalysing", ncol(eset.1), "samples\n"))
}
dim(targs.1)
cbind(colnames(eset.1),  targs.1$Target)
unique(targs.1$Target)
    
## limma users guide example pg 50:
## Treat <- factor(paste(targets$Condition,targets$Tissue,sep="."))
## design <- model.matrix(~0+Treat)
## colnames(design) <- levels(Treat)
## Then we estimate the correlation between measurements made on the same subject:
## corfit <- duplicateCorrelation(eset,design,block=targets$Subject)
## corfit$consensus
## Then this inter-subject correlation is input into the linear model fit:
## fit <- lmFit(eset,design,block=targets$Subject,correlation=corfit$consensus)
## end example


# # # # # # #  set up design using Buffy's parameterization # # # # # # # #
# design.1 <- model.matrix(~ Treatment + Antagomir)
Target.fac.1 <- factor(targs.1$Target,
                       levels = level.labs)
design.1 <- model.matrix(~0 + Target.fac.1, data=targets)
as.data.frame( design.1)

## give nicer column names
colnames(design.1) <- levels(Target.fac.1)
design.1

## estimate the correlation between measurements made on the same subject:
## install.packages("statmod", lib = .libPaths()[1])
corfit <- duplicateCorrelation(eset.1, design.1, block = targs.1$Patient)
corfit$consensus                        # [1] 0.1012354 ## all timepoints
                                        # [1] 0.1612525 ## timepoints 1 & 2

## fit the primary model
fit.1 <- lmFit(eset.1, design.1, block = targs.1$Patient,
               correlation = corfit$consensus)
 
## ## ## ## ## ## ## get the standard errors for Nandini  ## ## ## ## ## ## 
## subgroup the data
healthy.01 <- eset.1[, targs.1$Target == "healthy.01"]
healthy.02 <- eset.1[, targs.1$Target == "healthy.02"]
if (timepoints == "all") {
    healthy.03 <- eset.1[, targs.1$Target == "healthy.03"]
    healthy.04 <- eset.1[, targs.1$Target == "healthy.04"]
}
asthma.01 <- eset.1[, targs.1$Target == "asthma.01"]
asthma.02 <- eset.1[, targs.1$Target == "asthma.02"]
CF.01 <- eset.1[, targs.1$Target == "CF.01"]
CF.02 <- eset.1[, targs.1$Target == "CF.02"]

## list it up and calculate  std. dev.
if (timepoints == "only1and2") {
    groups.list <- list(healthy.01, healthy.02, asthma.01,
                        asthma.02, CF.01, CF.02)
} else if (timepoints == "all") {
    groups.list <- list(healthy.01, healthy.02, healthy.03, healthy.04,
                        asthma.01, asthma.02, CF.01, CF.02)
}
sds <- lapply(groups.list, function(x) apply(x, 1, sd) )
names(sds) <- unique(targs.1$Target)

## calculate global std. dev.
g.sd <- apply(eset.1, 1, sd)

sds.df <- data.frame(do.call(cbind, sds), "all.sd" = g.sd)
head(sds.df)

## write.csv(sds.df, file = file.path(root, "results", "tables",
##             paste(use, scal.meth, timepoints,"t.points",
##                   "standard.devs.csv", sep=".")),
##           row.names = TRUE)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## define contrast(s) for priamry question(s) of interest
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

cont.mat.1 <- makeContrasts(
    "healthy.vs.asthma"        = (healthy.01 + healthy.02) -
                                 (asthma.01 + asthma.02),
    "healthy.01.vs.asthma.01"  = healthy.01 - asthma.01,
    ## "healthy.02.vs.asthma.02"  = healthy.02 - asthma.02,
    "healthy.vs.CF"            = (healthy.01 + healthy.02) -
                                 (CF.01 + CF.02),
    "healthy.01.vs.CF.01"      = healthy.01 - CF.01,
    ## "healthy.02.vs.CF.02"      = healthy.02 - CF.02,
    "asthma.vs.CF"             = (asthma.01 + asthma.02) -
                                 (CF.01 + CF.02),
    "asthma.01.vs.CF.01"       = asthma.01 - CF.01,
    levels = design.1
    )
cont.mat.1
fit.1.2 <- contrasts.fit(fit.1, cont.mat.1)
fit.1.3 <- eBayes(fit.1.2)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## define contrast(s) for secondary question(s) of interest
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## reproducibiltiy T1 vs T2 on all goups
if (timepoints == "all") {
    repro.cont.1.2 <- makeContrasts(
        ## reproducibiltiy T1 vs T2 on all goups
        "healthy.01.vs.healthy.02" = healthy.01 - healthy.02,
        "asthma.01.vs.asthma.02"   = asthma.01 - asthma.02,
        "CF.01.vs.CF.02"           = CF.01 - CF.02,
        "t01.vs.t02"                 = (healthy.01 + asthma.01 + CF.01) -
                                     (healthy.02 + asthma.02 + CF.02),
        levels = design.1
        )
    repro.cont.1.2
    fit.repro.1.2 <- contrasts.fit(fit.1, repro.cont.1.2)
    fit.repro.1.2.e <- eBayes(fit.repro.1.2)
## reproducivibiliy T1, T2, T3, vs T4 in healthy
    repro.cont.h1to4 <- makeContrasts(healthy.01 - healthy.02, 
                                      healthy.01 - healthy.03,
                                      healthy.01 - healthy.04,
                                      healthy.02 - healthy.03, 
                                      healthy.02 - healthy.04, 
                                      healthy.03 - healthy.04,
                                      levels = design.1
                                      ) 
    repro.cont.h1to4
    fit.repro.h1to4 <- contrasts.fit(fit.1, repro.cont.h1to4)
    fit.repro.h1to4.e <- eBayes(fit.repro.h1to4)
}



###########################################################################
#                               healthy.vs.asthma                         #
###########################################################################
tt.healthy.vs.asthma  <- topTable(fit.1.3, coef="healthy.vs.asthma",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.healthy.vs.asthma)
tt.healthy.vs.asthma[order(tt.healthy.vs.asthma$P.Value), ][1:20, ]
sum(tt.healthy.vs.asthma$adj.P.Val < 0.05)
# 1255     # NA   # 546 # 321 # 347 # 917* # 900** 
# pilot    # pil. # all # all # all # 
# smoothed # raw  # raw # MAD # quantile
# *all + smoothed + quantile normed
# **all + smoothed + quant.norm, log2 = T
# ***all 4 healthy time points included
# post meeting: 856** 762***
## hist(tt.healthy.vs.asthma$adj.P.Val)
## dev.off()
head(tt.healthy.vs.asthma)
write.table(tt.healthy.vs.asthma,
          file = file.path(root, "P-values of CO2H2O calibrated data",
              paste(use, scal.meth, "healthy.vs.asthma.alldata", datadate,
                    zeros, "csv", sep=".")),sep="\t",
          row.names = TRUE)

## "healthy.1.vs.asthma.1" = healthy.1 - asthma.1,
tt.healthy.1.vs.asthma.1  <- topTable(fit.1.3,
                                      coef="healthy.01.vs.asthma.01",
                                      number=Inf, p.value=1, sort.by="none")
dim(tt.healthy.1.vs.asthma.1)
tt.healthy.1.vs.asthma.1[order(tt.healthy.1.vs.asthma.1$P.Value), ][1:20, ]
sum(tt.healthy.1.vs.asthma.1$adj.P.Val < 0.05)
## 68 # NA # 56 # 64 # 71 # 205 # 199 # 113 # 195
## hist(tt.healthy.1.vs.asthma.1$adj.P.Val)
## dev.off()
## write.csv(tt.healthy.1.vs.asthma.1,
##           file = file.path(root, "results", "tables",
##               paste(use, scal.meth, "healthy.1.vs.asthma.1.alldata.csv",
##                     sep=".")), row.names = TRUE)

###########################################################################
#                               healthy.vs.CF                             #
###########################################################################
tt.healthy.vs.CF  <- topTable(fit.1.3, coef="healthy.vs.CF",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.healthy.vs.CF)
tt.healthy.vs.CF[order(tt.healthy.vs.CF$P.Value), ][1:20, ]
sum(tt.healthy.vs.CF$adj.P.Val < 0.05)
# 36 # NA # 353 # 387 # 483 # 1109 # 1097 # 1066 # 888
## hist(tt.healthy.vs.CF$adj.P.Val)
## dev.off()
write.table(tt.healthy.vs.CF,
          file = file.path(root, "P-values of CO2H2O calibrated data",
              paste(use, scal.meth, "healthy.vs.CF.alldata", datadate,
                    zeros, "csv", sep=".")),sep="\t",
          row.names = TRUE)

## ## healthy.1.vs.CF.1 ## ##
tt.healthy.1.vs.CF.1  <- topTable(fit.1.3, coef="healthy.01.vs.CF.01",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.healthy.1.vs.CF.1)
tt.healthy.1.vs.CF.1[order(tt.healthy.1.vs.CF.1$P.Value), ][1:20, ]
sum(tt.healthy.1.vs.CF.1$adj.P.Val < 0.05)
# [1] 0 #    # 59 # 58 # 119 # 560 # 559 # 493 # 498
sum(tt.healthy.1.vs.CF.1$adj.P.Val < 0.10)
# [1] 0                # 235 # 778 # 775 # 763 # 808
## hist(tt.healthy.1.vs.CF.1$adj.P.Val)
## dev.off()
## write.csv(tt.healthy.vs.CF,
##           file = file.path(root, "results", "tables",
##                paste(use, scal.meth, "healthy.1.vs.CF.1.alldata.csv",
##                      sep=".")), row.names = TRUE)

###########################################################################
#                              asthma.vs.CF                               #
###########################################################################
tt.asthma.vs.CF  <- topTable(fit.1.3, coef="asthma.vs.CF",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.asthma.vs.CF)
tt.asthma.vs.CF[order(tt.asthma.vs.CF$P.Value), ][1:20, ]
sum(tt.asthma.vs.CF$adj.P.Val < 0.05) # [1] 0 #  #  #  #  #  # 16 # 13 # 23
sum(tt.asthma.vs.CF$adj.P.Val < 0.10) # [1] 0 #  #  #  #  #  # 59 # 53 # 57
## hist(tt.asthma.vs.CF$adj.P.Val)
## hist(tt.asthma.vs.CF$P.Value)
## dev.off()
## write.csv(tt.asthma.vs.CF, file = file.path(root, "results", "tables",
##             paste(use, scal.meth, "asthma.vs.CF.alldata.csv", sep=".")),
##           row.names = TRUE)

tt.asthma.1.vs.CF.1  <- topTable(fit.1.3, coef="asthma.01.vs.CF.01",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.asthma.1.vs.CF.1)
tt.asthma.1.vs.CF.1[order(tt.asthma.1.vs.CF.1$P.Value), ][1:20, ]
sum(tt.asthma.1.vs.CF.1$adj.P.Val < 0.05)
                           # [1] 0 #  #  #  #  #  # 16 # 13 # 0 # 0
sum(tt.asthma.1.vs.CF.1$adj.P.Val < 0.10)
                           # [1] 0 #  #  #  #  #  # 59 # 53 # 0 # 4
## hist(tt.asthma.1.vs.CF.1$adj.P.Val)
## hist(tt.asthma.1.vs.CF.1$P.Value)
## dev.off()
## write.csv(tt.asthma.1.vs.CF.1,
##           file = file.path(root, "results", "tables", paste(use,
##               scal.meth, "asthma.1.vs.CF.1.alldata.csv", sep=".")),
##           row.names = TRUE)

###########################################################################
#                        reproducibility btwn times                       #
###########################################################################
## ## "healthy.1.vs.healthy.2" = healthy.1 - healthy.2,
tt.heal.1.vs.heal.2  <-
    topTable(fit.repro.1.2.e, coef="healthy.01.vs.healthy.02", number=Inf,
             p.value=1, sort.by="none")
dim(tt.heal.1.vs.heal.2)
tt.heal.1.vs.heal.2[order(tt.heal.1.vs.heal.2$P.Value), ][1:20, ]
sum(tt.heal.1.vs.heal.2$adj.P.Val < 0.05)
# [1] 0 #     # 2  # 2  # 3  # 0 # 0
sum(tt.heal.1.vs.heal.2$adj.P.Val < 0.10)
# [1] 0 #     # 3  # 27 # 21 # 0 # 0
sum(tt.heal.1.vs.heal.2$adj.P.Val < 0.20)
## hist(tt.heal.1.vs.heal.2$adj.P.Val)
## dev.off()

## ## "asthma.1.vs.asthma.2"   = asthma.1 - asthma.2,
tt.asthma.1.vs.asthma.2  <-
    topTable(fit.repro.1.2.e, coef="asthma.01.vs.asthma.02", number=Inf,
             p.value=1, sort.by="none")
dim(tt.asthma.1.vs.asthma.2)
tt.asthma.1.vs.asthma.2[order(tt.asthma.1.vs.asthma.2$P.Value), ][1:20, ]
sum(tt.asthma.1.vs.asthma.2$adj.P.Val < 0.05) # [1] 0
sum(tt.asthma.1.vs.asthma.2$adj.P.Val < 0.10) # [1] 0
sum(tt.asthma.1.vs.asthma.2$adj.P.Val < 0.20) # [1] 0
## hist(tt.asthma.1.vs.asthma.2$adj.P.Val)
## dev.off()

## ## "CF.1.vs.CF.2"      = CF.1 - CF.2,
tt.CF.1.vs.CF.2  <- topTable(fit.repro.1.2.e, coef="CF.01.vs.CF.02",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.CF.1.vs.CF.2)
tt.CF.1.vs.CF.2[order(tt.CF.1.vs.CF.2$P.Value), ][1:20, ]
sum(tt.CF.1.vs.CF.2$adj.P.Val < 0.05)
# [1] 0 #     # 1  # 0 # 0 # 0 # 0
sum(tt.CF.1.vs.CF.2$adj.P.Val < 0.10)
# [1] 0 #     # 1  # 1 # 1 # 0 # 0
sum(tt.CF.1.vs.CF.2$adj.P.Val < 0.20)
# 0
## hist(tt.CF.1.vs.CF.2$adj.P.Val)
## dev.off()

## ## "t1.vs.t2" = (healthy.01 + asthma.01 + CF.01) -
##                 (healthy.02 + asthma.02 + CF.02)
tt.t1.vs.t2  <- topTable(fit.repro.1.2.e, coef="t01.vs.t02",
                         number=Inf, p.value=1, sort.by="none")
dim(tt.t1.vs.t2)
tt.t1.vs.t2[order(tt.t1.vs.t2$P.Value), ][1:20, ]
sum(tt.t1.vs.t2$adj.P.Val < 0.05) # [1] 0                    # 0 # 0 # 0 
sum(tt.t1.vs.t2$adj.P.Val < 0.10) # [1] 0                    # 1 # 0 # 0
sum(tt.t1.vs.t2$adj.P.Val < 0.20) # [1] 0                    # 1 # 0 # 0
## hist(tt.t1.vs.t2$adj.P.Val)
## dev.off()

###########################################################################
# correlations btwn times                                                 #
###########################################################################
dim(eset.1)
colnames(eset.1)


###########################################################################
#                                PCA plots                                #
###########################################################################

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
heal.T1 <- which(targs.1$Target == "healthy.01")
heal.T2 <- which(targs.1$Target == "healthy.02")
heal.T3 <- which(targs.1$Target == "healthy.03")
heal.T4 <- which(targs.1$Target == "healthy.04")
health.samp.1.2 <- c(heal.T1, heal.T2)
health.samp.all <- c(heal.T1, heal.T2, heal.T3, heal.T4)

asthma.T1 <- which(targs.1$Target == "asthma.01")
asthma.T2 <- which(targs.1$Target == "asthma.02")
asthma.samp <- c(asthma.T1, asthma.T2)

CF.T1 <- which(targs.1$Target == "CF.01")
CF.T2 <- which(targs.1$Target == "CF.02")
CF.samp <- c(CF.T1, CF.T2)

## head(eset.1)
## plot(density(na.omit(eset.1[, 13])))

grp1.1.2 <- rep.nas(eset.1, i = health.samp.1.2)
grp1.all <- rep.nas(eset.1, i = health.samp.all)
head(grp1.1.2)
head(grp1.all)

grp2  <- rep.nas(eset.1, i = asthma.samp)
head(grp2)

grp3 <- rep.nas(eset.1, i = CF.samp)
head(grp3)
ncol(grp1) + ncol(grp2) + ncol(grp3) == ncol(eset.1)

## corelations
health.1.2.cor <- cor(grp1.1.2, method = "spearman")
health.all.cor <- cor(grp1.all, method = "spearman")
asthma.cor <- cor(grp2, method = "spearman")
CF.cor     <- cor(grp3, method = "spearman")

###########################################################################
# correlation heatmaps                                                    #
###########################################################################

## setup some color blind safe colors
library(gplots) # used for for heatmap.2()
library(RColorBrewer)
blu <- colorRampPalette(brewer.pal(9, "Blues"))
blues <- rev(blu(256))
blu.yel <- colorRampPalette(c("blue", "yellow"))
blues.yel <- blu.yel(256)

## and the classic red-blue
blu.white.red <- colorRampPalette(c("blue", "white", "red"))
blu.white.red <- blu.white.red(256)

## blue.white.red heatmaps
## ## ## ## ## ## ## ## ## healthy T1 and T2 ## ## ## ## ## ## ## ## ## ##
pdf(file = file.path(root, "results", "plots", "PDFs",
        "heathy.1.2.spearman.cor.blue.pdf"), width = 8.27, height = 8.27)
## x11(w = 9, h = 9, xpos = 1200)
heatmap.2(health.1.2.cor, col =  blu.white.red, trace = "none",
          main = "Spearman's rho btween time 1 & 2 healthy patients",
          density.info = "none", margins = c(5, 5))
dev.off()

## ## ## ## ## ## ## ## ## asthma T1 and T2 ## ## ## ## ## ## ## ## ## ## #
pdf(file = file.path(root, "results", "plots", "PDFs",
        "asthma.spearman.cor.blue.pdf"), width = 8.27, height = 8.27)
## x11(w = 9, h = 9, xpos = 1200)
heatmap.2(asthma.cor, col =  blu.white.red, trace = "none",
          main = "Spearman's rho btween all asthma patients",
          density.info = "none", margins = c(5, 5))
dev.off()

## ## ## ## ## ## ## ## ## CF T1 and T2 ## ## ## ## ## ## ## ## ## ## ## ##
pdf(file = file.path(root, "results", "plots", "PDFs",
        "CF.spearman.cor.blue.pdf"), width = 8.27, height = 8.27)
## x11(w = 9, h = 9, xpos = 1200)
heatmap.2(CF.cor, col =  blu.white.red, trace = "none",
          main = "Spearman's rho btween all CF patients",
          density.info = "none", margins = c(5, 5))
dev.off()



## PCA 3 main groups #####################################################
## bind all gropus back together
eset.1.nas <-na.omit(cbind(grp1, grp2, grp3))
nrow(eset.1.nas)
pca1 <- prcomp(x=t(eset.1.nas))

## pdf(file=file.path(root, "results", "plots", "PDFs",
##       paste("van_mastrigt.3grp.PCA.plot", use, ".pdf")))
x11(xpos=2000)
plot(pca1$x[,1:2], type="n")
points(pca1$x[,1:2], bg=c(rep("red", length(health.samp)),
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
dev.off()

## PCA 6 groups ###########################################################
grp1 <- rep.nas(eset.1, i = heal.T1)
grp2 <- rep.nas(eset.1, i = heal.T2)
grp3 <- rep.nas(eset.1, i = asthma.T1)
grp4 <- rep.nas(eset.1, i = asthma.T2)
grp5 <- rep.nas(eset.1, i = CF.T1)
grp6 <- rep.nas(eset.1, i = CF.T2)
## bind all gropus back together
eset.2.nas <-na.omit(cbind(grp1, grp2, grp3, grp4, grp5, grp6))
nrow(eset.2.nas)
pca2 <- prcomp(x=t(eset.2.nas))

## pdf(file=file.path(root, "results", "plots", "PDFs",
##       paste("van_mastrigt.6grp.PCA.plot", use, ".pdf")))
x11(xpos=2000)
plot(pca2$x[,1:2], type="n")
points(pca2$x[,1:2], bg=c(rep("red", length(heal.T1)),
                         rep("red4", length(heal.T2)),
                         rep("blue", length(asthma.T1)),
                         rep("lightblue", length(asthma.T2)),
                         rep("green", length(CF.T1)),
                         rep("lightgreen", length(CF.T1))),
       col=c(rep("red", length(heal.T1)),
             rep("red4", length(heal.T2)),
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

## PCA 3 time 1 groups ####################################################
grp1 <- rep.nas(eset.1, i = heal.T1)
## grp2 <- rep.nas(eset.1, i = heal.T2)
grp3 <- rep.nas(eset.1, i = asthma.T1)
## grp4 <- rep.nas(eset.1, i = asthma.T2)
grp5 <- rep.nas(eset.1, i = CF.T1)
## grp6 <- rep.nas(eset.1, i = CF.T2)
## bind all gropus back together
eset.3.nas <-na.omit(cbind(grp1, grp3, grp5))
nrow(eset.3.nas)
pca3 <- prcomp(x=t(eset.3.nas))

## pdf(file=file.path(root, "results", "plots", "PDFs",
##       paste("van_mastrigt.6grp.PCA.plot", use, ".pdf")))
x11(xpos=2000)
plot(pca3$x[,1:2], type="n")
points(pca3$x[,1:2], bg=c(rep("red", length(heal.T1)),
                         rep("blue", length(asthma.T1)),
                         rep("green", length(CF.T1))),
       col=c(rep("red", length(heal.T1)),
             rep("blue", length(asthma.T1)),
             rep("green", length(CF.T1))),
       pch=21)
## add a legend
legend("topleft",
       legend=c("healthy.T1",
                "asthma.T1",
                "CF.T1"),
       fill=c("red",
              "blue",
              "green"))
## add a title
title(paste("PCA for", use, "data at time 1", sep=" "))
dev.off()

###########################################################################
#                               3D PCA plot                               #
###########################################################################
pca1 <- prcomp(x=t(eset.1.nas))

## pdf(file=file.path(root, "results", "plots", "PDFs",
##       paste("van_mastrigt.3grp.PCA.plot", use, ".pdf")))
x11(xpos=2000)
plot(pca1$x[,1:2], type="n")
points(pca1$x[,1:2], bg=c(rep("red", length(health.samp)),
                         rep("blue", length(asthma.samp)),
                         rep("green", length(CF.samp))),
                        col=c(rep("red", length(health.samp)),
                              rep("blue", length(asthma.samp)),
                              rep("green", length(CF.samp))),
                            pch=21)
dev.off()
### PCA analysis
# on all 96 samples (non-means)
## pca3d <- prcomp(t(all.ord.DE), scale.=TRUE)
pca3d <- prcomp(x=t(eset.1.nas), scale.=TRUE) ## 3 main groups
str(eset.1.nas)

#Plot scatter plot using "scatterplot3d"
library(scatterplot3d)
scatterplot3d(pca3d$x[,1:3], highlight.3d=TRUE )

## Plot scatter plot using "plot3d"
#Add text labels just to the mean postions- ! RUN FIRST !
library(rgl)
## pca.means <- prcomp(t(all.DE.means), scale.=TRUE)
## text3d(pca.means$x[,1:3], 
##        texts=c("R.S", "R.E", "R.L", "C.S", "C.E", "C.L"),
##        adj=c(0,0))

### make rotatable 3d plot - ! RUN SECOND !
plot3d(pca3d$x[,1:3], 
       col=c(rep("blue", length(health.samp)),
             rep("red", length(asthma.samp)),
             rep("green", length(CF.samp))),
       #size=8,
       type="s",
       radius=2)#,
       ## add=TRUE)    #specifiy to add to prelabelled plot


### Add box, axis, acis lables- ! RUN THIRD !
axes3d(edges = "bbox")#, labels = TRUE, tick = TRUE, nticks = 5, ...)

### Add axis labels (could also do it using title3d()- ! RUN FOURTH !
## decorate3d(pca.means$x[,1:3], xlab="PC1", ylab="PC2", zlab="PC3", 
##            #main="PCA of 96 samples ofr DE genes (q<0.05)",
##            box=FALSE,     #Suppress redraweing box
##            axes=TRUE,
##            sub=NULL, 
##            aspect=TRUE) #must be true to correct non-cube ratio, ie, make it a cube!

### Add grid- ! RUN FITH !
grid3d(c("x", "y+", "z")) #at = NULL, col = "gray", lwd = 1, lty = 1, n = 5)
grid3d(c("x+", "y+", "z"), at = NULL, col = "gray", lwd = 1, lty = 1, n = 5)


### Add (large) title- ! RUN SIX !
#par3d(cex=1.5)
title3d(main="PCA of 96 samples ofr DE genes (q<0.05)", line=5, cex=1.2)

#    decorate3d(xlab = NULL, ylab = NULL, zlab = NULL,
#               box = FALSE, axes = FALSE,
#               main= NULL,
#               sub="PCA of 96 samples ofr DE genes (q<0.05)")#, 
#               #aspect=TRUE) #must be true to correct non-cube ratio, ie, make it a cube!

##saveit
# dir4 <- "C:/Documents and Settings/kbrand/work/F/Data/P.period/analyses&QC/Fig.s/PCA/" #home
 dir4 <- "F:/Data/P.period/analyses&QC/Fig.s/PCA/" #work
#  rgl.postscript(filename=paste(dir4, "r_PCA_v1_atempt2.pdf", sep=""),
#                 fmt="pdf", drawText=TRUE)  #ps, eps, tex, pdf, svg, pgf
#  rgl.postscript(filename=paste(dir4, "r_PCA_v1_atempt2.zoomedout.pdf", sep=""),
#                 fmt="pdf", drawText=TRUE)  #ps, eps, tex, pdf, svg, pgf
#  rgl.postscript(filename=paste(dir4, "r_PCA_v1.svg", sep=""),
#                 fmt="svg", drawText=TRUE)
#  rgl.postscript(filename=paste(dir4, "r_PCA_v1.eps", sep=""), 
#                 fmt="eps", drawText=TRUE)
#  rgl.postscript(filename=paste(dir4, "r_PCA_v1.ps", sep=""), 
#                 fmt="ps", drawText=TRUE)
#  rgl.postscript(filename=paste(dir4, "r_PCA_v1.tex", sep=""), 
#                 fmt="tex", drawText=TRUE)
rgl.snapshot(filename   = "3D.PCA.1.png", fmt="png", top=TRUE )



###########################################################################
#! we should verify thatthe lables are indeed in the correct palce 
#! relative to ploted points
###########################################################################

### PCA on means
pca.means <- prcomp(t(all.DE.means), scale.=TRUE)
scatterplot3d(pca.means$x[,1:3], 
              color=c("#FFEDA0", "#FEB24C", "#F03B20", 
                      "#DEEBF7", "#9ECAE1", "#3182BD"))

rgl.open()
rgl.clear()
#rgl.clear("lights")
plot3d(pca.means$x[,1:3], 
       col=c("#FFEDA0", "#FEB24C", "#F03B20", "#DEEBF7", "#9ECAE1", "#3182BD"),
       type="s") 
       #radius=8)
text3d(pca.means$x[,1:3], 
       texts=c("R.S", "R.E", "R.L", "C.S", "C.E", "C.L"),
       adj=c(1.5,1.5))
decorate3d(pca.means$x[,1:3], xlab="PC1", ylab="PC2", zlab="PC3", 
           #main = "PCA of 96 samples ofr DE genes (q<0.05)",
           aspect=TRUE) #must be true to correct non-cube ratio, ie, make it a cube!
       





     


 


###choosing colors
library(RColorBrewer)
library(gplots)

showpanel <- function(col)
      {
        image(z=matrix(1:100, ncol=1), col=col,
xaxt="n", yaxt="n" )
      }

display.brewer.all(3)

showpanel(rep(brewer.pal(3, "YlOrRd"), 16))
showpanel(rep(brewer.pal(3, "Blues"), 16))


gc()
## script finishes here ###################################################
## deprecated code: #######################################################

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


## ###############################################
## # MAD scaling (modified from Wijnen analysis) #
## ###############################################
## if (scal.meth == "MAD") {
##     ## function which subtracts the sample median, divides by the absolute
##     ## sample median then adds back the gobal median to all results
##     eset[1:5, 1:5]
##     dim(eset)
##     mad.scale <- function(x, method=match.fun(median), na.rm = TRUE, ...) {
##         ## x <- eset           # used for debugging
##         ## scl <- median       # used for debugging
##         ## i <- 1              # used for debugging
##         scl <- method
##         ## global mean/median 
##         glob.me <- scl(x, na.rm=na.rm, ...)
##         print(glob.me)
##         y <- apply(x, 2, function(x) {
##             ## calcualte residuals by subtracting sample med. from each feature
##             scal.pl <- x - scl(x, na.rm=na.rm, ...)
##             ## median distance away from the median (MAD) ...as
##             ## opposed to the square-root of the average distance form
##             ## the mean normalise the variance of each plate
##             scal.pl.v <- scal.pl/scl(abs(scal.pl), na.rm=na.rm, ...)
##             ## scal.pl.v <- scal.pl/(sqrt(var(scal.pl)))
##             ## add global median to all
##             scal.pl.gl <- scal.pl.v + glob.me
##             return(scal.pl.gl)
##         }
##                    )
##         return(y)
##     }

##     eset.s <- mad.scale(eset, median)
