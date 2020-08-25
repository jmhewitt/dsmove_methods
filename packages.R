library(conflicted)
library(dotenv)
library(drake)

devtools::document('packages/dsmovetools/')
library(ggplot2)
library(ggthemes)
library(igraph)
library(metR)
library(nimble)
library(sp)
library(spdep)
library(viridis)

# compile likelihood function
cctds_nbhd_ll = compileNimble(ctds_nbhd_ll)