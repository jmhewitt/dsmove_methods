library(conflicted)
library(dotenv)
library(drake)

devtools::document('packages/dsmovetools/')
library(ctmcmove)
library(fda)
library(ggplot2)
library(ggthemes)
library(igraph)
library(metR)
library(nimble)
library(sp)
library(spdep)
library(viridis)

# compile likelihood function if needed
if(!exists('cctds_nbhd_ll')) {
  cctds_nbhd_ll = compileNimble(ctds_nbhd_ll)
}
