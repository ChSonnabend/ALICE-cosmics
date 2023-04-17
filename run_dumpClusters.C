

R__LOAD_LIBRARY(dumpClusters_cxx.so);

void run_dumpClusters()
{
    printf("run_dumpClusters started \n");
    gSystem ->Load("dumpClusters_cxx.so");
    dumpClustersFromTracksITS("/Users/aschmah/alice/TPC_calibration/reco/tpctracksC.root","/Users/aschmah/alice/TPC_calibration/reco/tpc-native-clustersC.root","/Users/aschmah/alice/TPC_calibration/reco/o2clus_itsC.root");

}