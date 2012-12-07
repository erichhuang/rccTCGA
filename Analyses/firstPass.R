## firstPass.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## REQUIRE
require(synapseClient)
require(annotate)
require(org.Hs.eg.db)
require(Biobase)

## LOAD DATA
rccEnt <- loadEntity('syn345400')
rccEset <- rccEnt$objects$eset
rccRnaSeq <- exprs(rccEset)