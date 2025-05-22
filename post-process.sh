#!/bin/bash

#clusterfile="/groups/guttman/.aviti/AV233703/20240212_AV233703_MARKO1/demultiplex/Samples/MARKO/pipeline/workup/clusters/Aliquot2.clusters"
#clusterfile="/groups/guttman/.aviti/AV233703/20240212_AV233703_MARKO1/demultiplex/Samples/MARKO/pipeline/workup/clusters/Aliquot10.clusters"
#clusterfile="/groups/guttman/.aviti/AV233703/20240212_AV233703_MARKO1/demultiplex/Samples/MARKO/pipeline/workup/clusters/Aliquot9.clusters"
#clusterfile="/groups/guttman/.aviti/AV233703/20240212_AV233703_MARKO1/demultiplex/Samples/MARKO/pipeline/workup/clusters/Aliquot11.clusters"
clusterfile="/groups/guttman/.aviti/AV233703/20240212_AV233703_MARKO1/demultiplex/Samples/MARKO/pipeline/workup/clusters/Aliquot8.clusters"

#split cluster file by condition
awk 'BEGIN { FS = "." } ; $6 ~ /CNTRL/ {print $0 }' $clusterfile > $clusterfile".control.clusters"
awk 'BEGIN { FS = "." } ; $6 ~ /TORIN/ {print $0 }' $clusterfile > $clusterfile".torin.clusters"
