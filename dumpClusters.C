#define GPUCA_TPC_GEOMETRY_O2
#include <fmt/format.h>
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsTPC/ClusterNative.h"
#include "DataFormatsTPC/ClusterNativeHelper.h"
#include "DataFormatsTPC/Defs.h"
#include "DataFormatsTPC/Constants.h"
#include "TPCBase/Mapper.h"
#include "TPCCalibration/TrackDump.h"
#include "GPU/GPUO2Interface.h"
#include "GPU/GPUDefOpenCL12Templates.h"
#include "GPU/GPUDefConstantsAndSettings.h"
#include "GPU/GPUTPCGeometry.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "functions.h"

#include "AlTrackHitEvent.h"
#include "AlTrackHitLinkDef.h"

#include <iostream>
#include <array>
#include <algorithm>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGFrame.h>
#include <TGTab.h>
#include <TGLCameraOverlay.h>
#include <TEveFrameBox.h>
#include <TEveQuadSet.h>
#include <TEveTrans.h>
#include <TEvePointSet.h>
#include <TEveTrackPropagator.h>
#include <TEveTrack.h>
#include <TEveEventManager.h>
#include <TEveScene.h>

#include "EventVisualisationView/MultiView.h"

#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "ITSMFTReconstruction/DigitPixelReader.h"
#include "ITSMFTReconstruction/RawPixelReader.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "DataFormatsITSMFT/Digit.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITS/TrackITS.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"

#include "CommonDataFormat/TFIDInfo.h"
#include "CommonDataFormat/InteractionRecord.h"

#include <vector>

ClassImp(AlTPCCluster)
ClassImp(AlITSHit)
ClassImp(AlTrack)
ClassImp(AlTrackHitEvent)

using namespace o2::itsmft;

extern TEveManager* gEve;
static TEveScene* chipScene;

static TGNumberEntry* gEntry;
static TGNumberEntry* gChipID;

o2::itsmft::TopologyDictionary dict;


// o2::tpc::TrackDump::ClusterNativeAdd::loadCorrMaps("./spline/TPC/Calib/CorrectionMap/TPCFastTransform_VoxRes_0_2cm_in_XYZ.root")

using namespace std;
using namespace o2::tpc;
using namespace o2::tpc::constants;
using ClusterRefVec = std::vector<o2::tpc::TPCClRefElem>;
using ClExcludes = std::vector<int>[MAXSECTOR][MAXGLOBALPADROW];

struct ClInfo : public TrackDump::ClusterNativeAdd {
  ClInfo() = default;
  ClInfo(const ClInfo&) = default;
  ClInfo(const ClusterNative& cl) : ClusterNativeAdd(cl){};
  int trackID = -1;

  ClassDefNV(ClInfo, 1);
};

#pragma link C++ class ClInfo + ;
#pragma link C++ class std::vector < ClInfo> + ;
#pragma link C++ class std::vector < std::vector < std::vector < Int_t> > > + ;

void fillClInfo(ClusterNativeAccess const& clusterIndex, std::vector<ClInfo>& clInfos, ClExcludes* excludes = nullptr);

// .L dumpClusters.C++
// dumpClusters("tpc-native-clusters.root")
// dumpClusters("tpc-native-clusters_ITS_low_pT_full_ext0_z.root")
// dumpClustersFromTracks("tpctracks.root","tpc-native-clusters.root")
// dumpClustersFromTracks("tpctracks_ITS_low_pT_full_ext0_z_dxpos.root","tpc-native-clusters_ITS_low_pT_full_ext0_z_dxpos.root")
// dumpClustersFromTracks("tpctracks_no_corr.root","tpc-native-clusters_no_corr.root")
// dumpClustersFromTracks("tpctracks_ITS_low_pt_V16_dydz.root","tpc-native-clusters_ITS_low_pt_V16_dydz.root")
// dumpClustersFromTracks("tpctracks_vD_only.root","tpc-native-clusters_vD_only.root")

// dumpClustersFromTracksITS("/Users/aschmah/alice/TPC_calibration/reco/tpctracksC.root","/Users/aschmah/alice/TPC_calibration/reco/tpc-native-clustersC.root","/Users/aschmah/alice/TPC_calibration/reco/o2clus_itsC.root")



void dumpClusters(std::string_view file = "tpc-native-clusters.root")
{
  std::vector<ClInfo> clInfos;

  o2::tpc::ClusterNativeHelper::Reader tpcClusterReader;
  tpcClusterReader.init(file.data());

  o2::tpc::ClusterNativeAccess clusterIndex;
  std::unique_ptr<o2::tpc::ClusterNative[]> clusterBuffer;
  // o2::tpc::MCLabelContainer clusterMCBuffer;
  o2::tpc::ClusterNativeHelper::ConstMCLabelContainerViewWithBuffer clusterMCBuffer;

  memset(&clusterIndex, 0, sizeof(clusterIndex));

  o2::utils::TreeStreamRedirector pcstream("ClusterPos.root", "recreate");

  TH2D* h2D_y_vs_x_cluster = new TH2D("h2D_y_vs_x_cluster","h2D_y_vs_x_cluster",2000,-50,50,300,85,255);
  TH1D* h_cluster_radius   = new TH1D("h_cluster_radius","h_cluster_radius",300,80,260);

  printf("Entries: %zu \n",tpcClusterReader.getTreeSize());

  for(size_t i = 0; i < tpcClusterReader.getTreeSize(); ++i)
  //for (size_t i = 0; i < 1; ++i)
  {
    clInfos.clear();
    tpcClusterReader.read(i);
    tpcClusterReader.fillIndex(clusterIndex, clusterBuffer, clusterMCBuffer);
    fillClInfo(clusterIndex, clInfos);

    // AliceO2/Detectors/TPC/calibration/include/TPCCalibration/TrackDump.h
    // struct ClusterNativeAdd : public ClusterNative

    Int_t size_cls = (Int_t)clInfos.size();
    //printf("size: %d \n",size_cls);
    for(Int_t i_cls = 0; i_cls < size_cls; i_cls++)
    {
        //cout << clInfos[i_cls].timeFlagsPacked << endl;
        //printf("gxy: {%4.3f, %4.3f}, gxyc: {%4.3f, %4.3f}, lxy: {%4.3f, %4.3f}, lxyc: {%4.3f, %4.3f} \n",clInfos[i_cls].gx(),clInfos[i_cls].gy(),clInfos[i_cls].gxc(),clInfos[i_cls].gyc(),clInfos[i_cls].lx(),clInfos[i_cls].ly(),clInfos[i_cls].lxc(),clInfos[i_cls].lyc());
        Double_t radius = clInfos[i_cls].lx();
        h2D_y_vs_x_cluster ->Fill(clInfos[i_cls].ly(),clInfos[i_cls].lx());
        h_cluster_radius   ->Fill(radius);
    }

    fmt::print("Read event {} with {} converted clusters\n", i, clInfos.size());

    pcstream << "cl"
             << "cls=" << clInfos
             << "\n";
  }

  //pcstream.Close();


  //-------------------------------------------------------
  h2D_y_vs_x_cluster ->GetXaxis()->SetTitle("x (cm)");
  h2D_y_vs_x_cluster ->GetYaxis()->SetTitle("y (cm)");
  h2D_y_vs_x_cluster ->GetZaxis()->SetTitle("entries");
  TCanvas* can_h2D_y_vs_x_cluster = Draw_2D_histo_and_canvas((TH2D*)h2D_y_vs_x_cluster,"can_h2D_y_vs_x_cluster",1100,600,0.0,0.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
  can_h2D_y_vs_x_cluster->SetLogz(1);
  //-------------------------------------------------------



  //-------------------------------------------------------
  h_cluster_radius ->SetLineColor(kBlack);
  h_cluster_radius ->GetYaxis()->SetTitleOffset(1.0);
  h_cluster_radius ->GetXaxis()->SetTitle("radius (cm)");
  h_cluster_radius ->GetYaxis()->SetTitle("entries");
  //h_cluster_radius ->GetYaxis()->SetRangeUser(-4.5,4.5);
  TCanvas* can_h_cluster_radius  = Draw_1D_histo_and_canvas((TH1D*)h_cluster_radius,"can_h_cluster_radius",800,500,0,0,"hist"); //TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option)
  //-------------------------------------------------------

}

void fillClInfo(ClusterNativeAccess const& clusterIndex, std::vector<ClInfo>& clInfos, ClExcludes* excludes)
{
    GPUCA_NAMESPACE::gpu::GPUTPCGeometry gpuGeom;

    for(int sector = 0; sector < MAXSECTOR; ++sector)
    {
        for(int padrow = 0; padrow < MAXGLOBALPADROW; ++padrow)
        {
            for(size_t icl = 0; icl < clusterIndex.nClusters[sector][padrow]; ++icl)
            {
                if(excludes)
                {
                    const auto& exRow = (*excludes)[sector][padrow];
                    //printf("begin: %d, end: %d, icl: %zu \n",exRow.begin(), exRow.end(), icl);
                    if(std::find(exRow.begin(), exRow.end(), icl) != exRow.end())
                    {
                        //printf("icl: %zu \n",icl);
                        continue;
                    }
                }
                const auto& cl = clusterIndex.clusters[sector][padrow][icl];
                auto& clInfo = clInfos.emplace_back(cl);
                clInfo.sector = sector;
                clInfo.padrow = padrow;
            }
        }
    }
}

void dumpClustersFromTracks(std::string_view trackFile = "tpctracks.root", std::string_view clusterFile = "tpc-native-clusters.root")
{
  GPUCA_NAMESPACE::gpu::GPUTPCGeometry gpuGeom;

  auto file = TFile::Open(trackFile.data());
  auto tree = (TTree*)file->Get("tpcrec");
  if(tree == nullptr)
  {
    std::cout << "Error getting tree\n";
    return;
  }

  // ---| branch setup |--------------------------------------------------------
  std::vector<o2::tpc::TrackTPC>* tpcTracks = nullptr;
  tree->SetBranchAddress("TPCTracks", &tpcTracks);
  ClusterRefVec* clusRefPtr = nullptr;
  tree->SetBranchAddress("ClusRefs", &clusRefPtr);

  std::vector<ClInfo> clInfos;

  o2::tpc::ClusterNativeHelper::Reader tpcClusterReader;
  tpcClusterReader.init(clusterFile.data());

  o2::tpc::ClusterNativeAccess clusterIndex;
  std::unique_ptr<o2::tpc::ClusterNative[]> clusterBuffer;
  // o2::tpc::MCLabelContainer clusterMCBuffer;
  o2::tpc::ClusterNativeHelper::ConstMCLabelContainerViewWithBuffer clusterMCBuffer;

  memset(&clusterIndex, 0, sizeof(clusterIndex));

  //o2::utils::TreeStreamRedirector pcstream("ClusterPosTracks.root", "recreate");

  TH1D* h_clusters_vs_padrow = new TH1D("h_clusters_vs_padrow","h_clusters_vs_padrow",155,0,154);
  TH1D* h_clusters_vs_padrow_all = new TH1D("h_clusters_vs_padrow_all","h_clusters_vs_padrow_all",155,0,154);
  TProfile* TP_clusters_vs_padrow = new TProfile("TP_clusters_vs_padrow","TP_clusters_vs_padrow",155,0,154);
  TH1D* h_time = new TH1D("h_time","h_time",1000,0,100000);
  TH2D* h_cls_y_vs_x = new TH2D("h_cls_y_vs_x","h_cls_y_vs_x",200,-250,250,200,-250,250);




  //------------------------------------------------
  TFile* outputfile = new TFile("Cluster_points.root","RECREATE");
  AlTPCCluster    *TPCCluster;
  AlITSHit        *ITSHit;
  AlTrack         *TPCTrack;
  AlTrackHitEvent *TrackHitEvent;
  TTree           *Tree_TrackHitEvent;

  TrackHitEvent       = new AlTrackHitEvent();
  TPCTrack            = new AlTrack();
  TPCCluster          = new AlTPCCluster();
  ITSHit              = new AlITSHit();
  Tree_TrackHitEvent  = NULL;
  Tree_TrackHitEvent  = new TTree("TrackHitEvent" , "TrackHitEvent" );
  Tree_TrackHitEvent  ->Branch("Events"  , "TrackHitEvent", TrackHitEvent);
  //------------------------------------------------




  TTree output_tree("cls_tree","a simple Tree with simple variables");
  Float_t id, x, y, time, row, sec;
  output_tree.Branch("id",&id);
  output_tree.Branch("x",&x);
  output_tree.Branch("y",&y);
  output_tree.Branch("time",&time);
  output_tree.Branch("row",&row);
  output_tree.Branch("sec",&sec);
  // SetAlias("clZ","(250-(clusters.getTime() - mTime0)*0.2*2.58)*(1-2*((clusters.sector%36)>17))");

  TTree output_tree_tracks("tracks_tree","TPC tracks");
  Float_t trkid, trkX, alpha, par0, par1, par2, par3, par4, time0;
  output_tree_tracks.Branch("trkid",&trkid);
  output_tree_tracks.Branch("trkX",&trkX);
  output_tree_tracks.Branch("alpha",&alpha);
  output_tree_tracks.Branch("par0",&par0);
  output_tree_tracks.Branch("par1",&par1);
  output_tree_tracks.Branch("par2",&par2);
  output_tree_tracks.Branch("par3",&par3);
  output_tree_tracks.Branch("par4",&par4);
  output_tree_tracks.Branch("time0",&time0);


  // === event loop |===========================================================
  int nEntries = tree->GetEntriesFast();
  int sum_clusters = 0;
  for(int i = 0; i < nEntries; ++i)
  {
      TrackHitEvent  ->clearTPCTrackList();
      TrackHitEvent  ->clearTPCClusterList();

      clInfos.clear();
      tree->GetEntry(i);
      tpcClusterReader.read(i);
      tpcClusterReader.fillIndex(clusterIndex, clusterBuffer, clusterMCBuffer);

#if 1
      // XALEX
      fillClInfo(clusterIndex, clInfos);
      Int_t size_clsAll = (Int_t)clInfos.size();
      //printf("N clusters: %d \n",size_clsAll);
      for(Int_t i_cls = 0; i_cls < size_clsAll; i_cls++)
      {
          Double_t radius = clInfos[i_cls].lx();
          Double_t cls_time = clInfos[i_cls].getTime();
          h_cls_y_vs_x ->Fill(clInfos[i_cls].gx(),clInfos[i_cls].gy());

          x = clInfos[i_cls].gx();
          y = clInfos[i_cls].gy();
          //float zc = clInfos[i_cls].zc();

          time = clInfos[i_cls].getTime();
          id = -1;
          sec = (Float_t)clInfos[i_cls].sector;
          row = clInfos[i_cls].padrow;
          //printf("row: %4.3f, sector: %d \n",row,sec);
          //output_tree.Fill();

          TPCCluster = TrackHitEvent->createTPCCluster();
          TPCCluster ->set_cluster_pos_time(x,y,time);
          TPCCluster ->set_track_id(id);
          TPCCluster ->set_sector(sec);
          TPCCluster ->set_row(row);
      }
#endif

      clInfos.clear();

      const size_t nTracks = tpcTracks->size();

      fmt::print("Processing event {} with {} tracks ({} MC indices)\n", i, nTracks, clusterMCBuffer.first.getIndexedSize());

      ClExcludes excludes;
      // ---| track loop |---
      for(size_t k = 0; k < nTracks; k++)
      {
          // AliceO2/DataFormats/Reconstruction/include/ReconstructionDataFormats/TrackParametrization.h
          const auto& track = (*tpcTracks)[k];
          const int nCl = track.getNClusterReferences();
          sum_clusters += nCl;
          float pt_track  = track.getPt();
          float eta_track = track.getEta();

          trkid = k;
          trkX  = track.getX();
          alpha = track.getAlpha();
          par0  = track.getParam(0);
          par1  = track.getParam(1);
          par2  = track.getParam(2);
          par3  = track.getParam(3);
          par4  = track.getParam(4);
          time0 = track.getTime0();
          //output_tree_tracks.Fill();

          TPCTrack = TrackHitEvent->createTPCTrack();
          TPCTrack ->set_track_id(k);
          TPCTrack ->set_X(track.getX());
          TPCTrack ->set_alpha(track.getAlpha());
          TPCTrack ->set_par(track.getParam(0),track.getParam(1),track.getParam(2),track.getParam(3),track.getParam(4));
          TPCTrack ->set_time0(track.getTime0());
          TPCTrack ->clearTPCClusterList();


          int padrow_cls[152] = {0};
          for(int j = nCl - 1; j >= 0; j--)
          {
              uint8_t sector, padrow;
              uint32_t clusterIndexInRow;
              track.getClusterReference(*clusRefPtr, j, sector, padrow, clusterIndexInRow);
              const auto& cl = clusterIndex.clusters[sector][padrow][clusterIndexInRow];
              excludes[sector][padrow].emplace_back(clusterIndexInRow);
              // const auto& cl = track.getCluster(*clusRefPtr, j, clusterIndex, sector, padrow);


              auto& clInfo = clInfos.emplace_back(cl);

              if(pt_track > 0.6)
              {
                  h_clusters_vs_padrow ->Fill(padrow);
                  padrow_cls[padrow] = 1;
              }

              //printf("track: %zu, cluster: %d, padrow: %d \n",k,j,padrow);


              //printf("pad: %4.3f, eta: %4.3f \n",cl.getPad(),eta_track);
              //printf("track: %d, cluser: %d, lx: %4.3f \n",k,j,clInfos[(Int_t)clInfos.size()-1].lx());
              //printf("track: %zu, cluster: %d, lx: %4.3f, gx: %4.3f, gy: %4.3f, padrow: %d, time: %4.3f \n",k,j,clInfo.lx(),clInfo.gx(),clInfo.gy(),padrow,clInfo.getTime());
              //printf("cl: %4.3f \n",cl.lx());

              clInfo.sector = sector;
              clInfo.padrow = padrow;

              clInfo.trackID = k;

              TPCCluster = TPCTrack->createTPCCluster();
              TPCCluster ->set_cluster_pos_time(clInfo.gx(),clInfo.gy(),clInfo.getTime());
              TPCCluster ->set_track_id(k);
              TPCCluster ->set_sector(sector);
              TPCCluster ->set_row(padrow);
          }

          for(int ipad = 0; ipad < 152; ipad++)
          {
              TP_clusters_vs_padrow ->Fill(ipad,padrow_cls[ipad]);
          }
      }

      fillClInfo(clusterIndex, clInfos, &excludes);

      Int_t size_cls = (Int_t)clInfos.size();
      //printf("size: %d \n",size_cls);
      for(Int_t i_cls = 0; i_cls < size_cls; i_cls++)
      {
          //cout << clInfos[i_cls].timeFlagsPacked << endl;
          //printf("gxy: {%4.3f, %4.3f}, gxyc: {%4.3f, %4.3f}, lxy: {%4.3f, %4.3f}, lxyc: {%4.3f, %4.3f}, time: %4.3f \n",clInfos[i_cls].gx(),clInfos[i_cls].gy(),clInfos[i_cls].gxc(),clInfos[i_cls].gyc(),clInfos[i_cls].lx(),clInfos[i_cls].ly(),clInfos[i_cls].lxc(),clInfos[i_cls].lyc(),clInfos[i_cls].getTime());
          Double_t radius = clInfos[i_cls].lx();
          h_time ->Fill(clInfos[i_cls].getTime());
          h_clusters_vs_padrow_all ->Fill(clInfos[i_cls].padrow);

          x = clInfos[i_cls].gx();
          y = clInfos[i_cls].gy();
          time = clInfos[i_cls].getTime();
          id = clInfos[i_cls].trackID;
          row = clInfos[i_cls].padrow;
          sec = (Float_t)clInfos[i_cls].sector;
          //output_tree.Fill();

          //TPCCluster = TrackHitEvent->createTPCCluster();
          //TPCCluster ->set_cluster_pos_time(x,y,time);
          //TPCCluster ->set_track_id(id);
          //TPCCluster ->set_sector(sec);
          //TPCCluster ->set_row(row);
      }

      //pcstream << "cl"
      //         << "cls=" << clInfos
      //         << "\n";

      Tree_TrackHitEvent ->Fill();
  }

  printf("sum_clusters: %d \n",sum_clusters);

  //-------------------------------------------------------
  h_clusters_vs_padrow ->SetLineColor(kBlack);
  h_clusters_vs_padrow ->GetYaxis()->SetTitleOffset(1.0);
  h_clusters_vs_padrow ->GetXaxis()->SetTitle("pad row");
  h_clusters_vs_padrow ->GetYaxis()->SetTitle("entries");
  //h_clusters_vs_padrow ->GetYaxis()->SetRangeUser(-4.5,4.5);
  TCanvas* can_h_clusters_vs_padrow  = Draw_1D_histo_and_canvas((TH1D*)h_clusters_vs_padrow,"can_h_clusters_vs_padrow",800,500,0,0,"hist"); //TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option)
  h_clusters_vs_padrow_all ->SetLineColor(kRed);
  h_clusters_vs_padrow_all ->DrawCopy("same hist");
  //-------------------------------------------------------


  //-------------------------------------------------------
  TP_clusters_vs_padrow ->SetLineColor(kBlack);
  TP_clusters_vs_padrow ->GetYaxis()->SetTitleOffset(1.0);
  TP_clusters_vs_padrow ->GetXaxis()->SetTitle("pad row");
  TP_clusters_vs_padrow ->GetYaxis()->SetTitle("<cluster>");
  TP_clusters_vs_padrow ->GetYaxis()->SetRangeUser(0.0,0.18);
  TCanvas* can_TP_clusters_vs_padrow  = Draw_1D_histo_and_canvas((TH1D*)TP_clusters_vs_padrow,"can_TP_clusters_vs_padrow",800,500,0,0,"hist"); //TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option)
  //-------------------------------------------------------


  //-------------------------------------------------------
  h_time ->SetLineColor(kBlack);
  h_time ->GetYaxis()->SetTitleOffset(1.0);
  h_time ->GetXaxis()->SetTitle("time");
  h_time ->GetYaxis()->SetTitle("entries");
  //h_time ->GetYaxis()->SetRangeUser(-4.5,4.5);
  TCanvas* can_h_time  = Draw_1D_histo_and_canvas((TH1D*)h_time,"can_h_time",800,500,0,0,"hist"); //TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option)
  //-------------------------------------------------------



  //-------------------------------------------------------
  h_cls_y_vs_x ->GetXaxis()->SetTitle("x (cm)");
  h_cls_y_vs_x ->GetYaxis()->SetTitle("y (cm)");
  h_cls_y_vs_x ->GetZaxis()->SetTitle("entries");
  TCanvas* can_h_cls_y_vs_x = Draw_2D_histo_and_canvas((TH2D*)h_cls_y_vs_x,"can_h_cls_y_vs_x",1100,800,0.0,0.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
  can_h_cls_y_vs_x->SetLogz(0);
  //-------------------------------------------------------



  //-------------------------------------------------------
  outputfile  ->cd();
  //output_tree.Write();
  //output_tree_tracks.Write();
  Tree_TrackHitEvent ->Write();
  outputfile  ->Close();
  //-------------------------------------------------------



  //pcstream.Close();
}



void dumpClustersFromTracksITS(std::string_view trackFile = "tpctracks.root", std::string_view clusterFile = "tpc-native-clusters.root", std::string_view ITSclusterFile = "o2clus_its.root")
{

    //------------------------------------------------------------
    TFile* file_ITS_noisy_pixels;
    Float_t NITS_id, NITS_row, NITS_col;

    TFile* inputfile = TFile::Open("./file_ITS_noisy_pixels.root");
    TTree *input_tree = (TTree*)inputfile->Get("ITS_noisy");
    Float_t NITS_in_id, NITS_in_row, NITS_in_col;
    input_tree->SetBranchAddress("id",&NITS_in_id);
    input_tree->SetBranchAddress("row",&NITS_in_row);
    input_tree->SetBranchAddress("col",&NITS_in_col);
    //------------------------------------------------------------



    //----------------------------------------------------
    std::vector< std::vector< std::vector<Int_t> > > vec_ITS_id_row_column_noisy;
    vec_ITS_id_row_column_noisy.resize(24120);
    for(Int_t i_id = 0; i_id < 24120; i_id++)
    {
        vec_ITS_id_row_column_noisy[i_id].resize(512);
    }
    Int_t nentries_noisy_ITS = (Int_t)input_tree->GetEntries();
    printf("Noisy ITS pixels: %d \n",nentries_noisy_ITS);     // 21483
    for(Int_t ientry = 0; ientry < nentries_noisy_ITS; ientry++)
    {
        input_tree->GetEntry(ientry);
        if((Int_t)NITS_in_id < 0 || (Int_t)NITS_in_id >= 24120) printf("WARNING: NITS_in_id out of range! \n");
        if((Int_t)NITS_in_row < 0 || (Int_t)NITS_in_row >= 512) printf("WARNING: NITS_in_row out of range! \n");
        vec_ITS_id_row_column_noisy[(Int_t)NITS_in_id][(Int_t)NITS_in_row].push_back((Int_t)NITS_in_col);

        //printf("i_entry: %d, id: %d, row: %d, col: %d \n",ientry,(Int_t)NITS_in_id,(Int_t)NITS_in_row,(Int_t)NITS_in_col);
    }

    //printf("col noisy: %d \n",vec_ITS_id_row_column_noisy[286][150][0]);
    //i_hit: 7006, id: 286, row: 150, col: 206
    //----------------------------------------------------



    //----------------------------------------------------
    // ITS related
    printf("Load ITS data and geometry \n");
    //auto fileITS = TFile::Open("o2clus_its.root");
    auto fileITS = TFile::Open(ITSclusterFile.data());
    TTree* ITSclstree = (TTree*)fileITS->Get("o2sim");
    std::vector<CompClusterExt> *mClusterBuffer = nullptr;
    std::vector<o2::itsmft::ROFRecord> *mClustersROF = nullptr;
    ITSclstree->SetBranchAddress("ITSClusterComp", &mClusterBuffer);
    ITSclstree->SetBranchAddress("ITSClustersROF", &mClustersROF);

    gsl::span<CompClusterExt> mClusters;

    // Geometry
    o2::base::GeometryManager::loadGeometry("o2sim_geometry-aligned.root");
    auto gman = o2::its::GeometryTGeo::Instance();
    gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot,
                                                   o2::math_utils::TransformType::L2G));

    // https://github.com/AliceO2Group/AliceO2/blob/8397b143cb59d8df5c95ed70ea76b7f87f132d2b/EventVisualisation/Workflow/README.md
    // https://alice.its.cern.ch/jira/browse/O2-2288  -> ITSdictionary.bin
    dict.readFromFile(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(o2::detectors::DetID::ITS));

    Long64_t entries = ITSclstree ->GetEntries();

    std::vector<o2::dataformats::TFIDInfo>* tfids = nullptr;
    TFile* fInTFID = TFile::Open("o2_tfidinfo.root");
    if (fInTFID) {
        // for a simulation this file is not available
        tfids = (std::vector<o2::dataformats::TFIDInfo>*) fInTFID->Get("tfidinfo");
    }
    printf("ITS data and geometry loaded \n");
    //----------------------------------------------------



    GPUCA_NAMESPACE::gpu::GPUTPCGeometry gpuGeom;

    auto file = TFile::Open(trackFile.data());
    auto tree = (TTree*)file->Get("tpcrec");
    if(tree == nullptr)
    {
        std::cout << "Error getting tree\n";
        return;
    }

    // ---| branch setup |--------------------------------------------------------
    std::vector<o2::tpc::TrackTPC>* tpcTracks = nullptr;
    tree->SetBranchAddress("TPCTracks", &tpcTracks);
    ClusterRefVec* clusRefPtr = nullptr;
    tree->SetBranchAddress("ClusRefs", &clusRefPtr);

    std::vector<ClInfo> clInfos;

    o2::tpc::ClusterNativeHelper::Reader tpcClusterReader;
    tpcClusterReader.init(clusterFile.data());

    o2::tpc::ClusterNativeAccess clusterIndex;
    std::unique_ptr<o2::tpc::ClusterNative[]> clusterBuffer;
    // o2::tpc::MCLabelContainer clusterMCBuffer;
    o2::tpc::ClusterNativeHelper::ConstMCLabelContainerViewWithBuffer clusterMCBuffer;

    memset(&clusterIndex, 0, sizeof(clusterIndex));

    //o2::utils::TreeStreamRedirector pcstream("ClusterPosTracks.root", "recreate");

    TH1D* h_clusters_vs_padrow = new TH1D("h_clusters_vs_padrow","h_clusters_vs_padrow",155,0,154);
    TH1D* h_clusters_vs_padrow_all = new TH1D("h_clusters_vs_padrow_all","h_clusters_vs_padrow_all",155,0,154);
    TProfile* TP_clusters_vs_padrow = new TProfile("TP_clusters_vs_padrow","TP_clusters_vs_padrow",155,0,154);
    TH1D* h_time = new TH1D("h_time","h_time",1000,0,100000);
    TH2D* h_cls_y_vs_x = new TH2D("h_cls_y_vs_x","h_cls_y_vs_x",200,-250,250,200,-250,250);




    //------------------------------------------------
    TFile* outputfile = new TFile("Tree_Cluster_points_cosmics_V5b.root","RECREATE");
    AlTPCCluster    *TPCCluster;
    AlITSHit        *ITSHit;
    AlTrack         *TPCTrack;
    AlTrackHitEvent *TrackHitEvent;
    TTree           *Tree_TrackHitEvent;

    TrackHitEvent       = new AlTrackHitEvent();
    TPCTrack            = new AlTrack();
    TPCCluster          = new AlTPCCluster();
    ITSHit              = new AlITSHit();
    Tree_TrackHitEvent  = NULL;
    Tree_TrackHitEvent  = new TTree("TrackHitEvent" , "TrackHitEvent" );
    Tree_TrackHitEvent  ->Branch("Events"  , "TrackHitEvent", TrackHitEvent);
    //------------------------------------------------




    TTree output_tree("cls_tree","a simple Tree with simple variables");
    Float_t id, x, y, time, row, sec;
    output_tree.Branch("id",&id);
    output_tree.Branch("x",&x);
    output_tree.Branch("y",&y);
    output_tree.Branch("time",&time);
    output_tree.Branch("row",&row);
    output_tree.Branch("sec",&sec);
    // SetAlias("clZ","(250-(clusters.getTime() - mTime0)*0.2*2.58)*(1-2*((clusters.sector%36)>17))");

    TTree output_tree_tracks("tracks_tree","TPC tracks");
    Float_t trkid, trkX, alpha, par0, par1, par2, par3, par4, time0;
    output_tree_tracks.Branch("trkid",&trkid);
    output_tree_tracks.Branch("trkX",&trkX);
    output_tree_tracks.Branch("alpha",&alpha);
    output_tree_tracks.Branch("par0",&par0);
    output_tree_tracks.Branch("par1",&par1);
    output_tree_tracks.Branch("par2",&par2);
    output_tree_tracks.Branch("par3",&par3);
    output_tree_tracks.Branch("par4",&par4);
    output_tree_tracks.Branch("time0",&time0);


    // === event loop |===========================================================
    int nEntries = tree->GetEntriesFast();
    printf("Number of time frames: %d \n",nEntries);

    vector< vector<Float_t> > vec_TPC_clusters;
    vector<Float_t> vec_TPC_cluster_info;
    vec_TPC_cluster_info.resize(6); // x,y,time,id,sec,row

    vector< vector<Float_t> > vec_ITS_hits;
    vector<Float_t> vec_ITS_hit_info;
    vec_ITS_hit_info.resize(7); // x,y,z,time_in_BC,row,col,id


    int sum_clusters = 0;
    Int_t N_accepted_events = 0;
    for(int i = 0; i < nEntries; ++i)  // time frames
    {
        if(i%20 == 0) printf("TF: %d out of %d \n",i,nEntries);


        TrackHitEvent  ->clearTPCTrackList();
        TrackHitEvent  ->clearTPCClusterList();
        TrackHitEvent  ->clearITSHitList();





        clInfos.clear();
        tree->GetEntry(i);
        const size_t nTracks = tpcTracks->size();
        //printf("TF: %d out of %d, nTracks TPC : %zu \n",i,nEntries,nTracks);
        if(nTracks <= 0) continue;
        tpcClusterReader.read(i);
        tpcClusterReader.fillIndex(clusterIndex, clusterBuffer, clusterMCBuffer);

        // XALEX
        fillClInfo(clusterIndex, clInfos);
        Int_t size_clsAll = (Int_t)clInfos.size();
        //printf("N clusters: %d \n",size_clsAll);
        for(Int_t i_cls = 0; i_cls < size_clsAll; i_cls++)
        {
            Double_t radius = clInfos[i_cls].lx();
            Double_t cls_time = clInfos[i_cls].getTime();
            h_cls_y_vs_x ->Fill(clInfos[i_cls].gx(),clInfos[i_cls].gy());

            x = clInfos[i_cls].gx();
            y = clInfos[i_cls].gy();
            //float zc = clInfos[i_cls].zc();

            time = clInfos[i_cls].getTime();
            id = -1;
            sec = (Float_t)clInfos[i_cls].sector;
            row = clInfos[i_cls].padrow;
            //printf("row: %4.3f, sector: %d \n",row,sec);
            //output_tree.Fill();

            vec_TPC_cluster_info[0] = x;
            vec_TPC_cluster_info[1] = y;
            vec_TPC_cluster_info[2] = time;
            vec_TPC_cluster_info[3] = id;
            vec_TPC_cluster_info[4] = sec;
            vec_TPC_cluster_info[5] = row;
            vec_TPC_clusters.push_back(vec_TPC_cluster_info);

            //TPCCluster = TrackHitEvent->createTPCCluster();
            //TPCCluster ->set_cluster_pos_time(x,y,time);
            //TPCCluster ->set_track_id(id);
            //TPCCluster ->set_sector(sec);
            //TPCCluster ->set_row(row);
        }

        clInfos.clear();

        //fmt::print("Processing event {} with {} tracks ({} MC indices)\n", i, nTracks, clusterMCBuffer.first.getIndexedSize());

        ClExcludes excludes;
        // ---| track loop |---
        for(size_t k = 0; k < nTracks; k++)
        {
            // AliceO2/DataFormats/Reconstruction/include/ReconstructionDataFormats/TrackParametrization.h
            const auto& track = (*tpcTracks)[k];
            const int nCl = track.getNClusterReferences();
            sum_clusters += nCl;
            float pt_track  = track.getPt();
            float eta_track = track.getEta();

            trkid = k;
            trkX  = track.getX();
            alpha = track.getAlpha();
            par0  = track.getParam(0);
            par1  = track.getParam(1);
            par2  = track.getParam(2);
            par3  = track.getParam(3);
            par4  = track.getParam(4);
            time0 = track.getTime0();
            //output_tree_tracks.Fill();

            TPCTrack = TrackHitEvent->createTPCTrack();
            TPCTrack ->set_track_id(k);
            TPCTrack ->set_X(track.getX());
            TPCTrack ->set_alpha(track.getAlpha());
            TPCTrack ->set_par(track.getParam(0),track.getParam(1),track.getParam(2),track.getParam(3),track.getParam(4));
            TPCTrack ->set_time0(track.getTime0());
            TPCTrack ->clearTPCClusterList();


            int padrow_cls[152] = {0};
            for(int j = nCl - 1; j >= 0; j--)
            {
                uint8_t sector, padrow;
                uint32_t clusterIndexInRow;
                track.getClusterReference(*clusRefPtr, j, sector, padrow, clusterIndexInRow);
                const auto& cl = clusterIndex.clusters[sector][padrow][clusterIndexInRow];
                excludes[sector][padrow].emplace_back(clusterIndexInRow);
                // const auto& cl = track.getCluster(*clusRefPtr, j, clusterIndex, sector, padrow);


                auto& clInfo = clInfos.emplace_back(cl);

                if(pt_track > 0.6)
                {
                    h_clusters_vs_padrow ->Fill(padrow);
                    padrow_cls[padrow] = 1;
                }

                //printf("track: %zu, cluster: %d, padrow: %d \n",k,j,padrow);


                //printf("pad: %4.3f, eta: %4.3f \n",cl.getPad(),eta_track);
                //printf("track: %d, cluser: %d, lx: %4.3f \n",k,j,clInfos[(Int_t)clInfos.size()-1].lx());
                //printf("track: %zu, cluster: %d, lx: %4.3f, gx: %4.3f, gy: %4.3f, padrow: %d, time: %4.3f \n",k,j,clInfo.lx(),clInfo.gx(),clInfo.gy(),padrow,clInfo.getTime());
                //printf("cl: %4.3f \n",cl.lx());

                clInfo.sector = sector;
                clInfo.padrow = padrow;

                clInfo.trackID = k;

                TPCCluster = TPCTrack->createTPCCluster();
                TPCCluster ->set_cluster_pos_time(clInfo.gx(),clInfo.gy(),clInfo.getTime());
                TPCCluster ->set_track_id(k);
                TPCCluster ->set_sector(sector);
                TPCCluster ->set_row(padrow);
            }

            for(int ipad = 0; ipad < 152; ipad++)
            {
                TP_clusters_vs_padrow ->Fill(ipad,padrow_cls[ipad]);
            }
        }

        fillClInfo(clusterIndex, clInfos, &excludes);

#if 0
        Int_t size_cls = (Int_t)clInfos.size();
        //printf("size: %d \n",size_cls);
        for(Int_t i_cls = 0; i_cls < size_cls; i_cls++)
        {
            //cout << clInfos[i_cls].timeFlagsPacked << endl;
            //printf("gxy: {%4.3f, %4.3f}, gxyc: {%4.3f, %4.3f}, lxy: {%4.3f, %4.3f}, lxyc: {%4.3f, %4.3f}, time: %4.3f \n",clInfos[i_cls].gx(),clInfos[i_cls].gy(),clInfos[i_cls].gxc(),clInfos[i_cls].gyc(),clInfos[i_cls].lx(),clInfos[i_cls].ly(),clInfos[i_cls].lxc(),clInfos[i_cls].lyc(),clInfos[i_cls].getTime());
            Double_t radius = clInfos[i_cls].lx();
            h_time ->Fill(clInfos[i_cls].getTime());
            h_clusters_vs_padrow_all ->Fill(clInfos[i_cls].padrow);

            x = clInfos[i_cls].gx();
            y = clInfos[i_cls].gy();
            time = clInfos[i_cls].getTime();
            id = clInfos[i_cls].trackID;
            row = clInfos[i_cls].padrow;
            sec = (Float_t)clInfos[i_cls].sector;
            //output_tree.Fill();

            //TPCCluster = TrackHitEvent->createTPCCluster();
            //TPCCluster ->set_cluster_pos_time(x,y,time);
            //TPCCluster ->set_track_id(id);
            //TPCCluster ->set_sector(sec);
            //TPCCluster ->set_row(row);
        }

        //pcstream << "cl"
        //         << "cls=" << clInfos
        //         << "\n";
#endif

        //------------------------------------------------------------
        ITSclstree ->GetEntry(i);
        Int_t cls_size = mClusterBuffer ->size();
        Int_t ROF_size = mClustersROF   ->size();
        //printf("cls_size: %d \n",cls_size);
        for(Int_t i_ROF = 0; i_ROF < ROF_size; i_ROF++) // ITS read out frames
        {
            auto rof = (*mClustersROF)[i_ROF];
            auto ir = rof.getBCData();
            int64_t time_in_BC = (*tfids)[i].isDummy() ? ir.bc2ns() * 1e-3 : ir.differenceInBC(o2::InteractionRecord{0, (*tfids)[i].firstTForbit});
            //printf("time_in_BC: %lld \n",time_in_BC);

            //std::cout << "Orbit: " << ir.orbit << "  BC: " << ir.bc << '\n';
            Int_t first = rof.getFirstEntry();
            Int_t last = first + rof.getNEntries();
            //printf("first: %d, last: %d \n",first,last);

            mClusters = gsl::make_span(&(*mClusterBuffer)[first], last - first);
            //std::cout << "Number of ITSClusters: " << mClusters.size() << '\n';
            for(Int_t i_cls = 0; i_cls < (Int_t)mClusters.size(); i_cls++)
            {
                //const CompClusterExt& mClusters = (*mClusterBuffer)[i_cls];
                UShort_t row = mClusters[i_cls].getRow();
                UShort_t col = mClusters[i_cls].getCol();
                UShort_t cid = mClusters[i_cls].getChipID();
                //printf("row: %d \n",row);
                auto locC = dict.getClusterCoordinates(mClusters[i_cls]);
                //printf("loc: {%4.3f, %4.3f, %4.3f} \n",locC.X(),locC.Y(),locC.Z());
                auto id = mClusters[i_cls].getSensorID();
                //printf("id: %d \n",id);

                // https://github.com/AliceO2Group/AliceO2/blob/3d08e57ccfdbf823046e65470e26fb121d571926/Detectors/ITSMFT/ITS/macros/EVE/DisplayEventsComp.C
                const auto gloC = gman->getMatrixL2G(id) * locC;
                //printf("row: %d, point: {%4.3f, %4.3f, %4.3f} \n",row,gloC.X(), gloC.Y(), gloC.Z());
                //tpoints->SetNextPoint(gloC.X(), gloC.Y(), gloC.Z());

                vec_ITS_hit_info[0] = gloC.X();
                vec_ITS_hit_info[1] = gloC.Y();
                vec_ITS_hit_info[2] = gloC.Z();
                vec_ITS_hit_info[3] = time_in_BC;
                vec_ITS_hit_info[4] = row;
                vec_ITS_hit_info[5] = col;
                vec_ITS_hit_info[6] = id;
                vec_ITS_hits.push_back(vec_ITS_hit_info);

                //ITSHit = TrackHitEvent->createITSHit();
                //ITSHit ->set_cluster_pos(gloC.X(), gloC.Y(), gloC.Z());
                //ITSHit ->set_BC(time_in_BC);
                //ITSHit ->set_row_col_id(row,col,id);
            }
        }
        //------------------------------------------------------------



        //------------------------------------------------------------
        UShort_t NTPCTracks  = TrackHitEvent ->getNumTPCTrack();
        Int_t    NTPCCluster = TrackHitEvent ->getNumTPCCluster();
        //Int_t    NITSHit     = TrackHitEvent ->getNumITSHit();
        Int_t NITSHit = (Int_t)vec_ITS_hits.size();
        //printf("NTPCTracks: %d, NTPCCluster: %d, NITSHit: %d \n",NTPCTracks,NTPCCluster,NITSHit);


        // TPC track loop
        vector<TVector3> vec_TV3_dir;
        vector<TVector3> vec_TV3_base;
        TVector3 TV3_cluster;
        vector<Int_t>    vec_N_TPC_clusters;
        Float_t  track_p[5];
        vector<Float_t> vec_helix_par;
        vec_helix_par.resize(6);
        vector< vector<Float_t> > vec_cls_points;
        vector<Float_t> vec_cls_point;
        vec_cls_point.resize(3);
        B_field = -0.001; // -5.0 kG
        Int_t max_ITS_chip_id = 0;
        Int_t max_ITS_row     = 0;
        Int_t max_ITS_col     = 0;

        for(Int_t i_track = 0; i_track < NTPCTracks; i_track++)
        {
            TPCTrack = TrackHitEvent->getTPCTrack(i_track);
            Int_t    track_id    = TPCTrack->get_track_id();
            Float_t  track_X     = TPCTrack->get_X();
            Float_t  track_alpha = TPCTrack->get_alpha();
            for(Int_t i_par = 0; i_par < 5; i_par++)
            {
                track_p[i_par] = TPCTrack->get_par(i_par);
            }
            Float_t track_time0 = TPCTrack->get_time0();
            set_helix(track_X,track_alpha,track_p,B_field);
            for(Int_t i_par = 0; i_par < 6; i_par++)
            {
                vec_helix_par[i_par] = fHelix[i_par];
                //printf("i_par: %d, par: %4.3f \n",i_par,fHelix[i_par]);
            }

            Double_t pT_track_charge = 0.0;
            TLorentzVector TLV_helix_prim = get_TLV_helix(B_field,pT_track_charge);
            Double_t pt_track  = TLV_helix_prim.Pt();
            Double_t eta_track = TLV_helix_prim.Eta();
            Double_t phi_track = TLV_helix_prim.Phi();
            Double_t px_track  = TLV_helix_prim.Px();
            Double_t py_track  = TLV_helix_prim.Py();
            Double_t pz_track  = TLV_helix_prim.Pz();

#if 0
            Float_t phi_track_at_inner_wall = -999.0;
            evaluate_helix(0.0,track_pos);
            Double_t radius_track = TMath::Sqrt(track_pos[0]*track_pos[0] + track_pos[1]*track_pos[1]);
            Float_t track_path_add = +1.0;
            if(radius_track > 87.225) track_path_add = -1.0; // go inwards
            for(Float_t track_path = 0.0; fabs(track_path) < 350.0; track_path += track_path_add)
            {
                evaluate_helix(track_path,track_pos);
                radius_track = TMath::Sqrt(track_pos[0]*track_pos[0] + track_pos[1]*track_pos[1]);
                if(radius_track < (85.225+1.0) && radius_track > (85.225 - 1.0))
                {
                    phi_track_at_inner_wall = TMath::ATan2(track_pos[1],track_pos[0]);
                    break;
                }
            }

            Float_t phi_track_at_inner_wall_sector = -999.0;
            if(phi_track_at_inner_wall > -999.0)
            {
                phi_track_at_inner_wall *= TMath::RadToDeg(); // -180..180
                phi_track_at_inner_wall += 180.0; // 0..360
                Int_t sector = (Int_t)(phi_track_at_inner_wall / 20.0);
                phi_track_at_inner_wall_sector = phi_track_at_inner_wall - sector*20.0;
                //printf("phi: %4.3f, sector: %d, phi_sec: %4.3f \n",phi_track_at_inner_wall,sector,phi_track_at_inner_wall_sector);
            }
#endif

            Int_t    NTPCCluster_track =  TPCTrack->getNumTPCCluster();


            // Estimate straight lines for TPC tracks
            TVector3 TV3_beam_center;
            TV3_beam_center.SetXYZ(0.0,0.0,0.0);
            TPCCluster = TPCTrack->getTPCCluster(0);
            Float_t x_cls_A        = TPCCluster ->get_cluster_x();
            Float_t y_cls_A        = TPCCluster ->get_cluster_y();
            Float_t t_cls_A        = TPCCluster ->get_cluster_time();
            TVector3 TV3_base_track;
            TV3_base_track.SetXYZ(x_cls_A,y_cls_A,t_cls_A);
            TVector3 TV3_dir_trackA;
            TV3_dir_trackA.SetXYZ(x_cls_A,y_cls_A,t_cls_A);
            TPCCluster = TPCTrack->getTPCCluster(NTPCCluster_track-1);
            Float_t x_cls_B        = TPCCluster ->get_cluster_x();
            Float_t y_cls_B        = TPCCluster ->get_cluster_y();
            Float_t t_cls_B        = TPCCluster ->get_cluster_time();
            TVector3 TV3_dir_trackB;
            TV3_dir_trackB.SetXYZ(x_cls_B,y_cls_B,t_cls_B);
            TVector3 TV3_dir_track = TV3_dir_trackB - TV3_dir_trackA;

            vec_TV3_dir.push_back(TV3_dir_track);
            vec_TV3_base.push_back(TV3_base_track);
            vec_N_TPC_clusters.push_back(NTPCCluster_track);
            TVector3 TV3_dca_center = calculateDCA_vec_StraightToPoint(TV3_base_track,TV3_dir_track,TV3_beam_center);
            Double_t radius_xy = TMath::Sqrt(TV3_dca_center.X()*TV3_dca_center.X() + TV3_dca_center.Y()*TV3_dca_center.Y());

            //printf("i_track: %d, par: {%4.3f, %4.3f, %4.3f, %4.3f, %4.3f} \n",i_track,vec_helix_par[0],vec_helix_par[1],vec_helix_par[2],vec_helix_par[3],vec_helix_par[4]);
        } // end of TPC track loop


        // Check for TPC to TPC track matches
        Int_t flag_good_TPC_match = 0;
        for(Int_t i_trackA = 0; i_trackA < (Int_t)vec_TV3_dir.size(); i_trackA++)
        {
            for(Int_t i_trackB = (i_trackA+1); i_trackB < (Int_t)vec_TV3_dir.size(); i_trackB++)
            {
                if(i_trackA == i_trackB) continue;
                TVector3 TV3_dca_TPC_hit = calculateDCA_vec_StraightToPoint(vec_TV3_base[i_trackA],vec_TV3_dir[i_trackA],vec_TV3_base[i_trackB]);
                if(fabs(TV3_dca_TPC_hit.Z()) < 6.0 && TV3_dca_TPC_hit.Perp() < 3.0)
                {
                    flag_good_TPC_match = 1;
                }
            }
        }



        //printf("Track loop done \n");

        // ITS loop
        Int_t N_good_ITS_hits   = 0;
        Int_t N_bad_ITS_hits    = 0;
        Int_t i_point_ITS_match = 0;
        vector<TVector3> vec_TV3_ITS_hits;
        TVector3 TV3_ITS_hit;
        std::vector<Int_t>::iterator itterator;
        Int_t flag_good_ITS_match = 0;
        Int_t N_matched_ITS_hits = 0;
        for(Int_t i_hit = 0; i_hit < NITSHit; i_hit++)
        {
            //ITSHit = TrackHitEvent->getITSHit(i_hit);
            //Float_t  x_cls        = ITSHit ->get_cluster_x();
            //Float_t  y_cls        = ITSHit ->get_cluster_y();
            //Float_t  z_cls        = ITSHit ->get_cluster_z();
            //int64_t  BC_cls       = ITSHit ->get_BC();
            //UShort_t row          = ITSHit ->get_row();
            //UShort_t col          = ITSHit ->get_col();
            //UShort_t id           = ITSHit ->get_id();

            Float_t  x_cls        =  (Float_t)vec_ITS_hits[i_hit][0];
            Float_t  y_cls        =  (Float_t)vec_ITS_hits[i_hit][1];
            Float_t  z_cls        =  (Float_t)vec_ITS_hits[i_hit][2];
            int64_t  BC_cls       =  (int64_t)vec_ITS_hits[i_hit][3];
            UShort_t row          =  (UShort_t)vec_ITS_hits[i_hit][4];
            UShort_t col          =  (UShort_t)vec_ITS_hits[i_hit][5];
            UShort_t id           =  (UShort_t)vec_ITS_hits[i_hit][6];

            Float_t  time_bin     = 400.0 + ((Float_t)BC_cls)/8.0;
            if(id > max_ITS_chip_id) max_ITS_chip_id = id;
            if(row > max_ITS_row) max_ITS_row = row;
            if(col > max_ITS_col) max_ITS_col = col;

           // printf("size: %d \n",(Int_t)vec_ITS_id_row_column_noisy[id][row].size());

#if 0
            if(id == 770 && row == 119 && col == 366)
            //if((Int_t)vec_ITS_id_row_column_noisy[id][row].size() == 0)
            {
                printf("Match!, size: %d \n",(Int_t)vec_ITS_id_row_column_noisy[id][row].size());
                for(Int_t i_size = 0; i_size < (Int_t)vec_ITS_id_row_column_noisy[id][row].size(); i_size++)
                {
                    printf(" col: %d \n, ",vec_ITS_id_row_column_noisy[id][row][i_size]);
                }
            }
#endif
            Int_t flag_good_ITS_hit = 0;
            if((Int_t)vec_ITS_id_row_column_noisy[id][row].size() > 0)
            {
                //printf("sizeA: %d \n",sizeA);
                itterator = std::find(vec_ITS_id_row_column_noisy[id][row].begin(), vec_ITS_id_row_column_noisy[id][row].end(), col);
                //if(std::find(vec_ITS_id_row_column_noisy[id][row].begin(), vec_ITS_id_row_column_noisy[id][row].end(), (Int_t)col) != vec_ITS_id_row_column_noisy[id][row].end() )
                if( itterator != vec_ITS_id_row_column_noisy[id][row].end() )
                {
                    if(id == 770 && row == 119 && col == 366)
                    {
                        printf(" \n");
                        printf("i_hit: %d, id: %d, row: %d, col: %d, size: %d \n",i_hit,id,row,col,(Int_t)vec_ITS_id_row_column_noisy[id][row].size());
                        for(Int_t i_size = 0; i_size < (Int_t)vec_ITS_id_row_column_noisy[id][row].size(); i_size++)
                        {
                            printf(" col: %d, ",vec_ITS_id_row_column_noisy[id][row][i_size]);
                        }
                    }
                    // pixel flagged as noisy
                    N_bad_ITS_hits++;
                    continue;
                }
                else
                {
                    N_good_ITS_hits++;
                    flag_good_ITS_hit = 1;
                }
            }
            else
            {
                //printf("Good hit \n");
                N_good_ITS_hits++;
                flag_good_ITS_hit = 1;
            }

            vec_cls_point[0] = x_cls;
            vec_cls_point[1] = y_cls;
            vec_cls_point[2] = time_bin;

            TVector3 TV3_ITS_hit;
            TV3_ITS_hit.SetXYZ(x_cls,y_cls,time_bin);
            vec_TV3_ITS_hits.push_back(TV3_ITS_hit);

            // Check for TPC track to ITS hit matches
            if(flag_good_ITS_hit)
            {
                for(Int_t i_track = 0; i_track < (Int_t)vec_TV3_dir.size(); i_track++)
                {
                    TVector3 TV3_dca_ITS_hit = calculateDCA_vec_StraightToPoint(vec_TV3_base[i_track],vec_TV3_dir[i_track],TV3_ITS_hit);
                    if(fabs(TV3_dca_ITS_hit.Z()) < 6.0 && TV3_dca_ITS_hit.Perp() < 3.0)
                    {
                        if(vec_N_TPC_clusters[i_track] > 30)
                        {
                            flag_good_ITS_match = 1;
                            printf("Match at TF: %d, track: %d, ITS hit: {%4.3f, %4.3f, %4.3f} \n",i,i_track,x_cls,y_cls,time_bin);

                            ITSHit = TrackHitEvent->createITSHit();
                            ITSHit ->set_cluster_pos(x_cls,y_cls,z_cls);
                            ITSHit ->set_BC(BC_cls);
                            ITSHit ->set_row_col_id(row,col,id);
                            N_matched_ITS_hits++;
                        }
                    }
                }
            }
            //vec_cls_points.push_back(vec_cls_point);
        } // end of ITS hit loop
        //if(N_good_ITS_hits > 0) printf("N_good_ITS_hits: %d, N_bad_ITS_hits: %d \n",N_good_ITS_hits,N_bad_ITS_hits);
        //--------------------


        // Store only those TPC clusters which are close to a TPC track
        for(Int_t i_cls = 0; i_cls < (Int_t)vec_TPC_clusters.size(); i_cls++)
        {
            Float_t x    = vec_TPC_clusters[i_cls][0];
            Float_t y    = vec_TPC_clusters[i_cls][1];
            Float_t time = vec_TPC_clusters[i_cls][2];
            Int_t id     = (Int_t)vec_TPC_clusters[i_cls][3];
            Int_t sec    = (Int_t)vec_TPC_clusters[i_cls][4];
            Int_t row    = (Int_t)vec_TPC_clusters[i_cls][5];
            TV3_cluster.SetXYZ(x,y,time);

            for(Int_t i_track = 0; i_track < (Int_t)vec_TV3_dir.size(); i_track++)
            {
                TVector3 TV3_dca_cls_hit = calculateDCA_vec_StraightToPoint(vec_TV3_base[i_track],vec_TV3_dir[i_track],TV3_cluster);
                if(fabs(TV3_dca_cls_hit.Z()) < 6.0 && TV3_dca_cls_hit.Perp() < 3.0)
                {
                    TPCCluster = TrackHitEvent->createTPCCluster();
                    TPCCluster ->set_cluster_pos_time(x,y,time);
                    TPCCluster ->set_track_id(id);
                    TPCCluster ->set_sector(sec);
                    TPCCluster ->set_row(row);
                }
            }
        }
        //------------------------------------------------------------


        if(flag_good_ITS_match || flag_good_TPC_match)
        {
            Tree_TrackHitEvent ->Fill();
            printf("accepted event: %d with %d attached ITS hits \n",N_accepted_events,N_matched_ITS_hits);
            N_accepted_events++;
        }
        vec_TPC_clusters.clear();
        vec_ITS_hits.clear();
    } // end of time frame loop

    printf("sum_clusters: %d \n",sum_clusters);

    //-------------------------------------------------------
    h_clusters_vs_padrow ->SetLineColor(kBlack);
    h_clusters_vs_padrow ->GetYaxis()->SetTitleOffset(1.0);
    h_clusters_vs_padrow ->GetXaxis()->SetTitle("pad row");
    h_clusters_vs_padrow ->GetYaxis()->SetTitle("entries");
    //h_clusters_vs_padrow ->GetYaxis()->SetRangeUser(-4.5,4.5);
    TCanvas* can_h_clusters_vs_padrow  = Draw_1D_histo_and_canvas((TH1D*)h_clusters_vs_padrow,"can_h_clusters_vs_padrow",800,500,0,0,"hist"); //TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option)
    h_clusters_vs_padrow_all ->SetLineColor(kRed);
    h_clusters_vs_padrow_all ->DrawCopy("same hist");
    //-------------------------------------------------------


    //-------------------------------------------------------
    TP_clusters_vs_padrow ->SetLineColor(kBlack);
    TP_clusters_vs_padrow ->GetYaxis()->SetTitleOffset(1.0);
    TP_clusters_vs_padrow ->GetXaxis()->SetTitle("pad row");
    TP_clusters_vs_padrow ->GetYaxis()->SetTitle("<cluster>");
    TP_clusters_vs_padrow ->GetYaxis()->SetRangeUser(0.0,0.18);
    TCanvas* can_TP_clusters_vs_padrow  = Draw_1D_histo_and_canvas((TH1D*)TP_clusters_vs_padrow,"can_TP_clusters_vs_padrow",800,500,0,0,"hist"); //TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option)
    //-------------------------------------------------------


    //-------------------------------------------------------
    h_time ->SetLineColor(kBlack);
    h_time ->GetYaxis()->SetTitleOffset(1.0);
    h_time ->GetXaxis()->SetTitle("time");
    h_time ->GetYaxis()->SetTitle("entries");
    //h_time ->GetYaxis()->SetRangeUser(-4.5,4.5);
    TCanvas* can_h_time  = Draw_1D_histo_and_canvas((TH1D*)h_time,"can_h_time",800,500,0,0,"hist"); //TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option)
    //-------------------------------------------------------



    //-------------------------------------------------------
    h_cls_y_vs_x ->GetXaxis()->SetTitle("x (cm)");
    h_cls_y_vs_x ->GetYaxis()->SetTitle("y (cm)");
    h_cls_y_vs_x ->GetZaxis()->SetTitle("entries");
    TCanvas* can_h_cls_y_vs_x = Draw_2D_histo_and_canvas((TH2D*)h_cls_y_vs_x,"can_h_cls_y_vs_x",1100,800,0.0,0.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_h_cls_y_vs_x->SetLogz(0);
    //-------------------------------------------------------



    //-------------------------------------------------------
    outputfile  ->cd();
    //output_tree.Write();
    //output_tree_tracks.Write();
    Tree_TrackHitEvent ->Write();
    outputfile  ->Close();
    //-------------------------------------------------------



    //pcstream.Close();
}
