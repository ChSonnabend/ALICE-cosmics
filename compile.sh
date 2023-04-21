singularity shell -B /lustre /lustre/alice/ctf/container/singularity_o2deploy.sif<<\EOF
alienv -w /opt/alibuild/sw enter O2Physics::latest

root -l
.L dumpClusters.C++

EOF
