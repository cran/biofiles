## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"-----------------------
BiocStyle::latex(width = 80)

## ----mito-path----------------------------------------------------------------
mito.path <- system.file("extdata", "S_cerevisiae_mito.gb", package="biofiles")

## ----read-rds, eval=TRUE, echo=FALSE------------------------------------------
mito.rds <- system.file("extdata", "NC_001224.rds", package="biofiles")
mito <- biofiles::loadRecord(mito.rds)

## ----read-mito, eval=FALSE, echo=TRUE-----------------------------------------
#  mito <- biofiles::gbRecord(mito.path)
#  mito

## ----show.rds, eval=TRUE, echo=FALSE------------------------------------------
mito

## ----store-mito, eval=FALSE, echo=TRUE----------------------------------------
#  biofiles::saveRecord(mito)
#  rm(mito)
#  mito <- biofiles::loadRecord("NC_001224.rds")
#  biofiles::summary(mito, n = 3)

## ----display-store-mito, eval=TRUE, echo=FALSE--------------------------------
biofiles::summary(mito, n = 3)

## ----summarise----------------------------------------------------------------
biofiles::summary(mito)

## ----tabulate-----------------------------------------------------------------
biofiles::qualifTable(mito)
biofiles::featureTable(mito)

## ----extract-header-----------------------------------------------------------
biofiles::getAccession(mito)
biofiles::getDefinition(mito)
biofiles::getGeneID(mito)
biofiles::getOrganism(mito)
biofiles::getLength(mito)
biofiles::getComment(mito)

## ----extract-sequence---------------------------------------------------------
biofiles::getSequence(mito)

## ----extract-ft---------------------------------------------------------------
biofiles::ft(mito)

## ----extract-key--------------------------------------------------------------
cds <- biofiles::filter(mito, key = "CDS")
biofiles::summary(cds[1:2])

## ----extact-key2--------------------------------------------------------------
cds <- mito["CDS"]
biofiles::summary(cds[3:4])

## ----extract-range------------------------------------------------------------
f10000 <- biofiles::filter(mito, range = "..10000")
biofiles::summary(f10000)

## ----extract-product----------------------------------------------------------
cytb <- biofiles::filter(mito, key = "CDS", product = "^cytochrome b$")
cytb

## ----access-feature-----------------------------------------------------------
biofiles::start(cds[1:3])
biofiles::end(cds[1:3])
biofiles::span(cds[1:3])
biofiles::strand(cds[1:3])
biofiles::locusTag(cds[1:3])
biofiles::dbxref(cds[1:3])
biofiles::product(cds[1:3])
biofiles::translation(cds[1:3])

## ----access-sequence----------------------------------------------------------
biofiles::getSequence(cds[1:6])

## ----access-qualif------------------------------------------------------------
biofiles::qualif(cds[1:3])
biofiles::qualif(cds[1:3], which = c("gene", "locus_tag", "EC_number", "product", "db_xref.GeneID"))

## ----access-select------------------------------------------------------------
cols <- c("key", "gene", "locus_tag", "product")
biofiles::select(cds[1:4], .cols = cols)

## ----granges------------------------------------------------------------------
biofiles::ranges(cds)

## ----granges-more-------------------------------------------------------------
biofiles::ranges(cds, join = TRUE, include = c("gene", "product", "db_xref"))

## ----sessinfo-----------------------------------------------------------------
utils::sessionInfo()

