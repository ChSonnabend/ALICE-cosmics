R__LOAD_LIBRARY(/lustre/alice/ctf/scripts/ALICE-cosmics/dumpClusters_C.so);

void run_dumpClusters()
{
    printf("run_dumpClusters started \n");
    gSystem ->Load("/lustre/alice/ctf/scripts/ALICE-cosmics/dumpClusters_C.so");
    const char* directory = gSystem->Getenv("DIRECTORY_FILES");
    std::string tracks=directory, native=directory, clus_its=directory;
    tracks+="/tpctracks.root"; native+="/tpc-native-clusters.root"; clus_its+="/o2clus_its.root";
    dumpClustersFromTracksITS(tracks,native,clus_its);

}
