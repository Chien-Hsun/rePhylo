README
---
title: 'rePhylo: Re-investigating and improving phylogenies'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r cover, include = FALSE, echo = FALSE}
require(rePhylo)
```

```{r figure, echo = FALSE}
barc<-readRDS("data/barplotMulticf.rds")
plot(barc$groupings)
```



## Overview

This package is to improve phylogenies based on phylogenies. One of the main function is to select potential paralogs by given groupings / constraints, and the users can perform another turn of phylogenetic reconstructions for an improved phylogeny. 

## Description and notes
`cladeFilter` reports tips as potential paralogs that violate the groupings provided in `taxa` for each gene tree of `trees`. Users can provide `taxa` directly, or generate one by `concor.node` and `createGrouping`.
**IMPORTANT NOTE:** The grouping/constraints must be chosen very carefully.  Avoid any aribitrary selections of particular relatioships; suggested to use the well-confirmed relationship, unless for special purpose.

## Installation
```{r install}
#install.packages("rePhylo")

```


## Usage and examples
First get the datasets:

```{r prepare}
data("Brassidata")
trees <- Brassidata$trees
taxa <- Brassidata$taxaTable
ref <- Brassidata$ref
```


`trees` contains all gene trees to be investigated as a list.


`taxa` is a `data.frame` that gives grouping information. The first column must be tip names in any of the test trees (gene trees).


`ref` is a reference tree (when using **Method 1**)


#### Method 1

```{r method1}
#concorn <- concor.node(ref = ref, trees = trees, bp = c(0, 30, 50), getTreeNames = FALSE)
#groups <- createGrouping(concor = concorn, percent = 0.7)
#head(groups)
```

#### Method 2

```{r emthod2}
#taxa <- read.table("taxaTable.txt", as.is = TRUE, header = 1)
head(taxa)
```

## Details





