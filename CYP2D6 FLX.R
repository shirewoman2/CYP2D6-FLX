# CYP2D6 FLX      
# 7/18/14 LS


# Housekeeping --------------------------------------------
library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
library(seqinr)
library(gridExtra)
library(plyr)
library(stringr)

# Order of datasets:
# 1. FLXEposP4
# 2. FLXEnegP5
# 3. FLXEposU2
# 4. FLXEnegU3

Dataset <- c("FLXEposP4", "FLXEnegP5", "FLXEposU2", "FLXEnegU3")

MainDir <- "G:/Data/Metabolomics/Jessica/FLXInhibition"
setwd(MainDir)

Data.orig <- list()

Directory <- c("G:/Data/Metabolomics/Jessica/FLXInhibition/20140402_FluoxetinePlasma/20140402_FluoxetinePlasmaESI+",
               "G:/Data/Metabolomics/Jessica/FLXInhibition/20140402_FluoxetinePlasma/20140402_FluoxetinePlasmaESIneg",
               "G:/Data/Metabolomics/Jessica/FLXInhibition/20130809_FLXInhibitionESI+",
               "G:/Data/Metabolomics/Jessica/FLXInhibition/20130809_FLXInhibitionESI-")

File <- c("FLXEposP4 peak table - RT above 2 min.csv",
          "FLXEnegP5 peak table - RT above 2 min.csv",
          "FLXEposU2 peak table - RT above 2 min.csv",
          "FLXEnegU3 peak table - RT above 2 min.csv")

for (j in 1:4){
      setwd(Directory[j])
      Data.orig[[j]] <- read.csv(File[j], na.strings=c("", "NA"), row.names=1)
}

setwd(MainDir)
Metadata <- read.csv("CYP2D6 fluoxetine metadata.csv", skip=1, na.strings=c("","NA", "#N/A"))
Metadata$EposFile2 <- make.names(str_trim(Metadata$EposFile))
Metadata$EnegFile2 <- make.names(str_trim(Metadata$EnegFile))


# All Sample names
Samp.P <- as.character(Metadata$Sample[Metadata$Matrix == "plasma"])
Samp.U <- as.character(Metadata$Sample[Metadata$Matrix == "urine"])


# Functions, etc. --------------------------------------

# Theme for graphs made using ggplot2
ThemeLaura <- function (base_size = 12, base_family = "") {
      theme_gray(base_size = base_size, base_family = base_family) %+replace%
            theme(
                  axis.text = element_text(colour = "black"),
                  axis.title.x = element_text(colour = "black"),
                  axis.title.y = element_text(colour = "black", angle=0),
                  panel.background = element_rect(fill="white", color=NA),
                  panel.grid.minor.y = element_line(color="white"),
                  panel.grid.minor.x = element_line(color="white"),
                  panel.grid.major = element_line(colour = "white"),
                  plot.background = element_rect(fill="white"),
                  panel.border = element_rect(color="black", fill=NA),
                  strip.background = element_rect(color=NA, fill="white"),
                  plot.background=element_rect(fill=NA, colour=NA),
                  legend.background = element_rect(color=NA, fill="white"),
                  legend.key = element_rect(color=NA, fill="white")
            )
}

theme_set(ThemeLaura())


# function for writing csv files with a note on the top line.
my.write <- function(x, file, header, f = write.csv, ...){
      # create and open the file connection
      datafile <- file(file, open = 'wt')
      # close on exit
      on.exit(close(datafile))
      # if a header is defined, write it to the file
      if(!missing(header)) writeLines(header,con=datafile)
      # write the file using the defined function and required addition arguments
      f(x, datafile,...)
}

FileNote <- paste("This file was created on", Sys.Date(), "using the script CYP2D6 metoprolol plasma.R.")

# Function for adding 95% CI ellipses to ggplots.
require(proto)

StatEllipse <- proto(ggplot2:::Stat,
{
      required_aes <- c("x", "y")
      default_geom <- function(.) GeomPath
      objname <- "ellipse"
      
      calculate_groups <- function(., data, scales, ...){
            .super$calculate_groups(., data, scales,...)
      }
      calculate <- function(., data, scales, level = 0.95, segments = 51,...){ # change the confidence level here. I don't know how to make it more flexible where you specify it everytime.
            dfn <- 2
            dfd <- length(data$x) - 1
            if (dfd < 3){
                  ellipse <- rbind(c(NA,NA))
            } else {
                  require(MASS)
                  v <- cov.trob(cbind(data$x, data$y))
                  shape <- v$cov
                  center <- v$center
                  radius <- sqrt(dfn * qf(level, dfn, dfd))
                  angles <- (0:segments) * 2 * pi/segments
                  unit.circle <- cbind(cos(angles), sin(angles))
                  ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
            }
            
            ellipse <- as.data.frame(ellipse)
            colnames(ellipse) <- c("x","y")
            return(ellipse)
      }
}
)

stat_ellipse <- function(mapping=NULL, data=NULL, geom="path", position="identity", ...) {
      StatEllipse$new(mapping=mapping, data=data, geom=geom, position=position, ...)
}


# extract legend from a ggplot2 plot and put it where you want
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)}


# Data tidying -------------------------------------
# Subsetting to only take column names with clinical Samples plus MF, mz, RT, count
FileNames <- list(
      as.character(Metadata$EposFile2[Metadata$Matrix == "plasma"]),
      as.character(Metadata$EnegFile2[Metadata$Matrix == "plasma"]),
      as.character(Metadata$EposFile2[Metadata$Matrix == "urine"]),
      as.character(Metadata$EnegFile2[Metadata$Matrix == "urine"]))


for (j in 1:4){
      ifelse(j < 3,
             Samp <- Samp.P,
             Samp <- Samp.U)
      
      Data.orig[[j]] <- Data.orig[[j]][, c("MassFeature", "mz", "RTmin",  
                                  FileNames[[j]], "Count")]
      names(Data.orig[[j]]) <- c("MassFeature", "mz", "RT", Samp, "Count")
      
}


# Replacing zeroes with 100s so that we don't take log(0) in future steps
for (j in 1:4){
      ifelse(j < 3,
             Samp <- Samp.P,
             Samp <- Samp.U)
      
      Mat <- Data.orig[[j]][, Samp]
      Impute <- Mat == 0
      Mat[Impute] <- 100
      Data.orig[[j]][, Samp] <- Mat
      rm(Mat)
}


# Preprocessing the data -----------------------------------------
# Objects I'll need to have made in advance
Data.filtered <- list()
Count.filtered <- list()
TCCsum <- list()
Data.TCCnorm <- list()
Data.log <- list()


for (j in 1:4){
      ifelse(j < 3,
             Samp <- Samp.P,
             Samp <- Samp.U)
      
      NumSamp <- length(Samp)
      
      # Filtering by frequency
      # Set the minimum number of Samples for which you want a mass feature to
      # have been detected originally to be included.
      Cutoff <- ceiling(NumSamp*0.25)
      Data.filtered[[j]] <- subset(Data.orig[[j]], Data.orig[[j]]$Count
                                   >= Cutoff)
      Count.filtered[[j]] <- Data.filtered[[j]]$Count
      
            
      # Calculate the TCC sum.
      TCCsum[[j]] <- apply(Data.filtered[[j]][, Samp], 2, sum, na.rm=TRUE)
      
      # Normalize each Sample by the TCC sum.
      Data.TCCnorm[[j]] <- Data.filtered[[j]]
      for (s in 1:NumSamp){
            Data.TCCnorm[[j]][,s+4] <- (Data.filtered[[j]][,s+4]/
                                              TCCsum[[j]][s])*1e6
            rm(s)
      }
      
      # log10 transforming
      Data.log[[j]] <- data.frame("MassFeature" = Data.TCCnorm[[j]]$MassFeature,
                                  "mz" = Data.TCCnorm[[j]]$mz,
                                  "RT" = Data.TCCnorm[[j]]$RT,
                                  log10(Data.TCCnorm[[j]][, Samp]))
      
      
      # Writing the final file.
      setwd(Directory[j])
      my.write(Data.log[[j]], paste(Dataset[j], "preprocessed.csv"),
               header=FileNote, row.names=F)
      
}

setwd(MainDir)
save(Data.orig, Metadata, Samp.P, Samp.U, Directory, File, 
     Data.log, MainDir, file="CYP2D6 FLX main data.RData")


# Checking for Jessica's 443 compound ---------------------------------
M1 <- list()

for (j in 1:4){
      ifelse(j %in% c(1,3),
             MZ <- 443.3028 + 1.0073,
             MZ <- 443.3028 - 1.0073)
      
      M1[[j]] <- Data.log[[j]][(Data.log[[j]]$mz < (MZ + 0.05) & Data.log[[j]]$mz > (MZ - 0.05)), ]
      
}

# Only detected in ESI+ urine samples.

setwd(MainDir)

# Making graph of M1 abundance based on inhibition status.
M1 <- melt(M1[[3]], id.vars=c("MassFeature", "mz", "RT"), 
             value.name="Abundance", variable.name="Sample")

M1 <- join(M1, Metadata, by = "Sample", type = "left")
M1 <- arrange(M1, Subject, StudyDay, Tx)

M1.means <- ddply(M1, c("Subject", "Tx"), summarize, "MeanAbund" = mean(Abundance))
M1.means.wide <- dcast(M1.means, Subject ~ Tx, value.var = "MeanAbund")

# Paired t test to compare mean abundance between treatments
with(M1.means.wide, t.test(Baseline, Fluoxetine, paired = T))

# Paired t-test
# 
# data:  Baseline and Fluoxetine
# t = 3.3046, df = 9, p-value = 0.009163
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#       0.2459354 1.3133421
# sample estimates:
#       mean of the differences 
# 0.7796387 

ggplot(M1, aes(x = Tx, y = Abundance, color = Subject, group = Subject)) +
      geom_point() + geom_line() + facet_wrap(~ UC)
ggsave("Line plot of M1 by subject and facetted by urine collection.png", 
       width=9, height=5, dpi=300)

ggplot(M1.means, aes(x = Tx, y = MeanAbund, color = Subject, group = Subject)) +
      geom_point() + geom_line()
ggsave("Line plot by subject of mean abundance of M1 over three urine collections.png",
       width=6, height=5, dpi=300)


