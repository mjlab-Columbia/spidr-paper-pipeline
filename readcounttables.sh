#!/bin/bash
clusterfile=$1
echo "Working on Cluster file: $clusterfile"

#Filtering clusters for condition based on first barcode:
awk 'BEGIN { FS = "." } ; $7 ~ /HEALTHY/ {print $0}' $clusterfile > $clusterfile".cntrl.clusters"

awk 'BEGIN { FS = "." } ; $7 ~ /FUS/ {print $0 }' $clusterfile > $clusterfile".fus.clusters"

#Read and Bead Counts:
python /groups/guttman/projects/spidr/scripts/plot_assignment_distribution_RPM.py -c $clusterfile".cntrl.clusters"  --min_oligos 1 --proportion 0.8 --max_size 100 --outprefix $clusterfile".cntrldistribution_1_0.8_100"

python /groups/guttman/projects/spidr/scripts/plot_assignment_distribution_RPM.py -c $clusterfile".fus.clusters"  --min_oligos 1 --proportion 0.8 --max_size 100 --outprefix $clusterfile".fusdistribution_1_0.8_100"

grep -v Cluster $clusterfile".cntrldistribution_1_0.8_100.counts" | sort -k1,1  > $clusterfile".controldistribution_1_0.8_100.sorted.counts"

grep -v Cluster $clusterfile".fusdistribution_1_0.8_100.counts" | sort -k1,1  > $clusterfile".fusdistribution_1_0.8_100.sorted.counts"
