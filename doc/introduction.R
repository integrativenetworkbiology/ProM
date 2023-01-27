## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE------------------------------------------------
library(ProM)
library(Seurat)

## ------------------------------------------------------------------------
data(sc)
class(sc)
dim(sc)
data(bulk)
class(bulk)
dim(bulk)

## ------------------------------------------------------------------------
table(sc$majorlabel)
table(sc$minorlabel)
table(sc$sampleID)

## ------------------------------------------------------------------------
d<-prepareInput(bulk,sc,majorlabel = "majorlabel", minorlabel = "minorlabel", samplelabel = "sampleID")
summary(d)

## ------------------------------------------------------------------------
priorcellfreqM<-getcellfreqM(sc,clusterlabel="minorlabel")
dim(priorcellfreqM)
priorcellfreqM[,"SI_651"]

## ------------------------------------------------------------------------
ProMEstimate<-ProM(bulk=d$bulk,sc.basis=d$minorsc.basis,markers=d$markers,priorcellfreqM=priorcellfreqM)
dim(ProMEstimate)
ProMEstimate[,1:5]

## ------------------------------------------------------------------------
otherEstimate<-otherDeconvolution(bulk=d$bulk,sc.basis=d$minorsc.basis,markers=d$markers,sc=sc,method="MuSiC")
otherEstimate[,1:5]

