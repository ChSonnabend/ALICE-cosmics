
// static const Double_t kAlmost1=1. - Double_t(FLT_EPSILON);
// static const Double_t kAlmost0=Double_t(FLT_MIN);
// static const Double_t kB2C=-0.299792458e-3;
// static Float_t fHelix[9];
// static Float_t B_field = -0.001; // -5.0 kG
// static Float_t track_pos[3];

#include "functions.h"
#include "AlTrackHitEvent.h"
#include "AlTrackHitLinkDef.h"
#include "TPrincipal.h"

#include <algorithm>
#include <numeric>
#include <vector>

ClassImp(AlTPCCluster)
    ClassImp(AlITSHit)
        ClassImp(AlTrack)
            ClassImp(AlTrackHitEvent)

                template <typename T>
                vector<float> sort_indices(const vector<T> &v)
{

    vector<float> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2)
                { return v[i1] < v[i2]; });

    return idx;
}

void AnaTrackHitEvent(Long64_t N_events = -1, Int_t event_plot = 0, Int_t flag_ITS_noisy = 0, bool saveToTree = true)
{
    // To analyze cosmics events, output from dumpClusters.C
    //.L AnaTrackHitEvent.cc++
    // AnaTrackHitEvent(-2,0,0)

    // AnaTrackHitEvent(-2,1,0) -> central membrane, misalignment?
    // AnaTrackHitEvent(-2,23,0) -> central membrane, misalignment?
    // AnaTrackHitEvent(-2,26,0) -> central membrane, misalignment?

    saveToTree = saveToTree && (N_events == -1);

    printf("AnaTrackHitEvent started \n");

    // N_events:
    // -1 -> loop over all events
    // -2 -> loop only over event_plot
    // any number -> loop over this number of events

    // AnaTrackHitEvent(-2,686,0)

    // 519 -> multiple cosmics, two machtes
    // 1156 -> one match
    // 759 -> one match
    // 686 -> one match + 8 ITS hits

    //------------------------------------------------------------
    TPrincipal *PCA = new TPrincipal(3, "");
    int fill_index = 0;
    bool anomalous_sector = false;                                                                                                   // Changes in iterations, leave default as false
    vector<int> tpc_ncls, event_num, num_ITS_hits, anomalous_sectors;                                                                //{4, 8, 17, 26, 35};                                              // sector 4 is only for reference
    vector<float> chi2_pca, chi2_pca_xy, chi2_pca_time, sigma_values_tpc = {0.22, 0.22, 0.7}, sigma_values_its = {0.01, 0.01, 0.01}; /// sigma_x = sigma_y = 1mm, sigma_time=2
    vector<Int_t> crosses_anomalous_sectors, crosses_anomalous_ITS;
    vector<vector<float>> residuals_in_anomalous_sectors;
    vector<TVector3> direction_PCA, base_PCA;
    TChain *input_chain = NULL;
    input_chain = new TChain("TrackHitEvent", "TrackHitEvent");
    // TString addfile = "./Data/Tree_Cluster_points_cosmics_V3.root";
    // TString addfile = "./Data/Tree_Cluster_points_cosmics_V5c.root";
    // TString addfile = "./Data/Tree_Cluster_points_cosmics_V5e.root";
    TString addfile = "./Data/IchFindAlexToll1100.root";
    input_chain->AddFile(addfile.Data(), -1, "TrackHitEvent");
    Long64_t file_entries = input_chain->GetEntries();
    printf("entries: %lld \n", file_entries);

    anomalous_sectors.resize(36);
    std::iota(anomalous_sectors.begin(), anomalous_sectors.end(), 0);

    AlTrackHitEvent *TrackHitEvent = new AlTrackHitEvent();
    input_chain->SetBranchAddress("Events", &TrackHitEvent);

    AlTrack *TPCTrack = new AlTrack();
    AlTPCCluster *TPCCluster = new AlTPCCluster();
    AlITSHit *ITSHit = new AlITSHit();
    Float_t track_p[5];
    vector<vector<Float_t>> vec_track_par_helix;
    vector<Float_t> vec_helix_par;
    vec_helix_par.resize(6);
    //------------------------------------------------------------

    //------------------------------------------------------------
    TFile *file_ITS_noisy_pixels = NULL;
    Float_t NITS_id, NITS_row, NITS_col;
    TTree output_tree("ITS_noisy", "a simple Tree with simple variables");
    if (flag_ITS_noisy)
    {
        file_ITS_noisy_pixels = new TFile("file_ITS_noisy_pixels_V2.root", "RECREATE");
        output_tree.Branch("id", &NITS_id);
        output_tree.Branch("row", &NITS_row);
        output_tree.Branch("col", &NITS_col);
    }

    TFile *inputfile = TFile::Open("file_ITS_noisy_pixels.root");
    TTree *input_tree = (TTree *)inputfile->Get("ITS_noisy");
    Float_t NITS_in_id, NITS_in_row, NITS_in_col;
    input_tree->SetBranchAddress("id", &NITS_in_id);
    input_tree->SetBranchAddress("row", &NITS_in_row);
    input_tree->SetBranchAddress("col", &NITS_in_col);
    //------------------------------------------------------------

    //------------------------------------------------------------

    TString open_mode = "OPEN";
    if (saveToTree)
        open_mode = "RECREATE";
    TFile *qafileITSTPC = new TFile("QA_tree.root", open_mode, "Tree with QA histograms and data");
    TTree *data_tree_TPC = new TTree("TPC", "TPC alignment QA");
    TTree *data_tree_TPC_raw_woCrit = new TTree("TPC_raw_woCrit", "TPC raw information with DCA's");
    TTree *data_tree_TPC_raw_wCrit = new TTree("TPC_raw_wCrit", "TPC raw information with DCA's");
    TTree *data_tree_TPC_raw_all = new TTree("TPC_raw_all", "TPC raw information with DCA's");
    TTree *data_tree_ITS = new TTree("ITS", "ITS alignment QA");

    TH2D *h_cls_y_vs_x = new TH2D("h_cls_y_vs_x", "h_cls_y_vs_x", 1000, -250, 250, 1000, -250, 250);
    TH2D *h_cls_y_vs_x_coarse = new TH2D("h_cls_y_vs_x_coarse", "h_cls_y_vs_x_coarse", 300, -250, 250, 300, -250, 250);

    // ITS / TPC misalignment
    TH1D *angle_dist_PCA = new TH1D("angle_dist_PCA", "angle_dist_PCA", 10000, 0, TMath::Pi());
    TH1D *angle_dist = new TH1D("angle_dist", "angle_dist", 10000, 0, TMath::Pi());
    TH1D *outlier_events = new TH1D("outlier_events", "outlier_events", 50000, 0, 50000);
    TH1D *chi2_TPC_tracks = new TH1D("chi2_TPC_tracks", "chi2_TPC_tracks", 500, 0, 50);
    TH1D *chi2_TPC_tracks_xy = new TH1D("chi2_TPC_tracks_xy", "chi2_TPC_tracks_xy", 500, 0, 50);
    TH1D *chi2_TPC_tracks_time = new TH1D("chi2_TPC_tracks_time", "chi2_TPC_tracks_time", 500, 0, 50);
    TH1D *residuals_2Point_ITS_histo = new TH1D("residuals_2Point_ITS_histo", "residuals_2Point_ITS_histo", 5000, 0, 10);
    TH1D *residuals_PCAfit_ITS_histo = new TH1D("residuals_PCAfit_ITS_histo", "residuals_PCAfit_ITS_histo", 5000, 0, 10);
    TH1D *chi2_PCAfit_ITS_histo = new TH1D("chi2_PCAfit_ITS_histo", "chi2_PCAfit_ITS_histo", 100, 0, 20);
    TNtuple *residuals_ITS_dx = new TNtuple("residuals_ITS_dx", "residuals_ITS_dx", "x:y:dcax:color");
    TNtuple *residuals_ITS_dy = new TNtuple("residuals_ITS_dy", "residuals_ITS_dy", "x:y:dcay:color");
    TNtuple *residuals_ITS_dz = new TNtuple("residuals_ITS_dz", "residuals_ITS_dz", "x:y:dcaz:color");
    TNtuple *residuals_ITS_dr = new TNtuple("residuals_ITS_dr", "residuals_ITS_dr", "x:y:dcar:color");
    TNtuple *residuals_ITS_full = new TNtuple("residuals_ITS_full", "residuals_ITS_full", "x:y:z:dca_x:dca_y:dca_z:chi2:n_points");

    // TNtuple* tpc_raw_clusters = new TNtuple("tpc_raw_clusters", "tpc_raw_clusters", "x:y:time:sector:row");
    TNtuple *tpc_raw_dca_wCrit = new TNtuple("tpc_raw_dca_wCrit", "tpc_raw_dca_wCrit", "x:y:time:sector:row:dca_x:dca_y:dca_time");
    TNtuple *tpc_raw_dca_woCrit = new TNtuple("tpc_raw_dca_woCrit", "tpc_raw_dca_woCrit", "x:y:time:sector:row:dca_x:dca_y:dca_time");

    // TPC SC anomaly, sector 8, 17, 26, 35
    std::vector<TH1D *> anomalous_TPC_dca3D_histos;
    std::vector<TH1D *> anomalous_TPC_dcaRad_histos;
    TString histo_name;
    TH1D *anomalous_TPC_sector_events = new TH1D("anomalous_TPC_sector_events", "anomalous_TPC_sector_events", 50000, 0, 50000);
    for (int i = 0; i < anomalous_sectors.size(); i++)
    {
        histo_name.Form("residual_dca3D_sector%i_histo", anomalous_sectors[i]);
        anomalous_TPC_dca3D_histos.push_back(new TH1D(histo_name, histo_name, 1000, 0, 5));
        histo_name.Form("residual_dcaRad_sector%i_histo", anomalous_sectors[i]);
        anomalous_TPC_dcaRad_histos.push_back(new TH1D(histo_name, histo_name, 2000, -5, 5));
    }

    if (saveToTree)
    {
        for (int i = 0; i < anomalous_TPC_dca3D_histos.size(); i++)
        {
            histo_name.Form("residual_dca3D_sector%i_histo", anomalous_sectors[i]);
            data_tree_TPC->Branch(histo_name, "TH1D", &anomalous_TPC_dca3D_histos[i], 128000, 0);
            histo_name.Form("residual_dcaRad_sector%i_histo", anomalous_sectors[i]);
            data_tree_TPC->Branch(histo_name, "TH1D", &anomalous_TPC_dcaRad_histos[i], 128000, 0);
        }

        data_tree_TPC->Branch("outlier_events", "TH1D", &outlier_events, 128000, 0);
        data_tree_TPC->Branch("angle_dist_PCA", "TH1D", &angle_dist_PCA, 128000, 0);
        data_tree_TPC->Branch("angle_dist", "TH1D", &angle_dist, 128000, 0);
        data_tree_TPC->Branch("chi2_TPC_tracks", "TH1D", &chi2_TPC_tracks, 128000, 0);
        data_tree_TPC->Branch("chi2_TPC_tracks_xy", "TH1D", &chi2_TPC_tracks_xy, 128000, 0);
        data_tree_TPC->Branch("chi2_TPC_tracks_time", "TH1D", &chi2_TPC_tracks_time, 128000, 0);

        // data_tree_TPC_raw_all->Branch("tpc_raw_clusters", "TNtuple", &tpc_raw_clusters, 128000, 0);
        data_tree_TPC_raw_wCrit->Branch("tpc_raw_clusters_wCrit", "TNtuple", &tpc_raw_dca_wCrit, 128000, 0);
        data_tree_TPC_raw_woCrit->Branch("tpc_raw_clusters_woCrit", "TNtuple", &tpc_raw_dca_woCrit, 128000, 0);

        data_tree_ITS->Branch("residuals_2Point_ITS_histo", "TH1D", &residuals_2Point_ITS_histo, 128000, 0);
        data_tree_ITS->Branch("residuals_PCAfit_ITS_histo", "TH1D", &residuals_PCAfit_ITS_histo, 128000, 0);
        data_tree_ITS->Branch("chi2_PCAfit_ITS_histo", "TH1D", &chi2_PCAfit_ITS_histo, 128000, 0);
        data_tree_ITS->Branch("residuals_ITS_full", "TNtuple", &residuals_ITS_full, 128000, 0);
    }

    TGraph *tg_cls_y_vs_x = new TGraph();
    TGraph *tg_cls_x_vs_time = new TGraph();
    TGraph *tg_cls_y_vs_time = new TGraph();
    TGraph *tg_ITS_hit_x_vs_y = new TGraph();
    TGraph *tg_ITS_hit_x_vs_time = new TGraph();
    TGraph *tg_ITS_hit_y_vs_time = new TGraph();
    vector<TGraph *> vec_tg_cls_y_vs_x;
    vector<TGraph *> vec_tg_cls_x_vs_time;
    vector<TGraph *> vec_tg_cls_y_vs_time;
    //------------------------------------------------------------

    //------------------------------------------------------------
    vec_track_par_helix.clear();
    Long64_t start_events = 0;
    if (N_events == -1)
        N_events = file_entries;
    if (N_events == -2)
    {
        start_events = event_plot;
        N_events = event_plot + 1;
    }
    Float_t time_min = +9999999.0;
    Float_t time_max = -9999999.0;
    vector<vector<Float_t>> vec_cls_points;
    vector<Float_t> vec_cls_point;
    vec_cls_point.resize(3);
    Int_t max_ITS_chip_id = 0;
    Int_t max_ITS_row = 0;
    Int_t max_ITS_col = 0;

    printf("Define ITS pixels vectors \n");
    // 108 144 180 2688 3360 8232 9408 chips for each ITS layer, total 24120
    // 7 layers
    // 24120 chips, 512 rows, 1024 columns

    vector<vector<vector<Int_t>>> vec_ITS_id_row_column;
    vector<vector<vector<Int_t>>> vec_ITS_id_row_column_noisy;
    vector<vector<vector<Int_t>>> vec_ITS_id_row_counter;
    std::vector<Int_t>::iterator itterator;
    vec_ITS_id_row_column.resize(24120);
    vec_ITS_id_row_column_noisy.resize(24120);
    vec_ITS_id_row_counter.resize(24120);
    for (Int_t i_id = 0; i_id < 24120; i_id++)
    {
        vec_ITS_id_row_column[i_id].resize(512);
        vec_ITS_id_row_column_noisy[i_id].resize(512);
        vec_ITS_id_row_counter[i_id].resize(512);
    }

    Long64_t N_pixels_hit = 0;

    //------------------------------------------------------------------------------
    Int_t N_noisy_ITS_pixels = 0;
    if (flag_ITS_noisy)
    {
        printf("Start time frame loop -> determine noisy ITS pixels \n");
        for (Long64_t counter = 0; counter < N_events; counter++)
        {
            if (counter != 0 && counter % 10 == 0)
                cout << "." << flush;
            if (counter != 0 && counter % 100 == 0)
            {
                if ((file_entries - 0) > 0)
                {
                    Double_t event_percent = 100.0 * ((Double_t)(counter - 0)) / ((Double_t)(file_entries - 0));
                    cout << " " << counter << " (" << event_percent << "%) "
                         << "\n"
                         << "==> Processing data " << flush;
                }
            }

            if (!input_chain->GetEntry(counter)) // take the event -> information is stored in event
                break;                           // end of data chunk

            //--------------------
            // Event information
            UShort_t NTPCTracks = TrackHitEvent->getNumTPCTrack();
            Int_t NTPCCluster = TrackHitEvent->getNumTPCCluster();
            Int_t NITSHit = TrackHitEvent->getNumITSHit();
            // printf("NTPCTracks: %d, NTPCCluster: %d, NITSHit: %d \n",NTPCTracks,NTPCCluster,NITSHit);

            // ITS loop
            for (Int_t i_hit = 0; i_hit < NITSHit; i_hit++)
            {
                ITSHit = TrackHitEvent->getITSHit(i_hit);
                Float_t x_cls = ITSHit->get_cluster_x();
                Float_t y_cls = ITSHit->get_cluster_y();
                Float_t z_cls = ITSHit->get_cluster_z();
                int64_t BC_cls = ITSHit->get_BC();
                UShort_t row = ITSHit->get_row();
                UShort_t col = ITSHit->get_col();
                UShort_t id = ITSHit->get_id();
                Float_t time_bin = ((Float_t)BC_cls) / 8.0;
                if (id > max_ITS_chip_id)
                    max_ITS_chip_id = id;
                if (row > max_ITS_row)
                    max_ITS_row = row;
                if (col > max_ITS_col)
                    max_ITS_col = col;

                itterator = std::find(vec_ITS_id_row_column[id][row].begin(), vec_ITS_id_row_column[id][row].end(), col);
                if (itterator != vec_ITS_id_row_column[id][row].end())
                {
                    // pixel already in list
                    Int_t index = itterator - vec_ITS_id_row_column[id][row].begin();
                    vec_ITS_id_row_counter[id][row][index]++;
                    // printf("pixel at id: %d, row: %d, col: %d, has %d hits \n",id,row,col,vec_ITS_id_row_counter[id][row][index]);
                }
                else
                {
                    // pixel not in list yet
                    vec_ITS_id_row_column[id][row].push_back(col);
                    vec_ITS_id_row_counter[id][row].push_back(1);
                    N_pixels_hit++;
                }
            }
            //--------------------
        } // end if time frame loop

        for (Int_t i_id = 0; i_id < 24120; i_id++) // total number of ITS-2 chips
        {
            for (Int_t i_row = 0; i_row < 512; i_row++)
            {
                Int_t N_columns_hit = (Int_t)vec_ITS_id_row_column[i_id][i_row].size();
                for (Int_t i_hit = 0; i_hit < N_columns_hit; i_hit++)
                {
                    Int_t i_col = vec_ITS_id_row_column[i_id][i_row][i_hit];
                    Int_t N_hits = vec_ITS_id_row_counter[i_id][i_row][i_hit];
                    // printf("id: %d, row: %d, hit: %d, col: %d, N_hits: %d \n",i_id,i_row,i_hit,i_col,N_hits);
                    if (N_hits > 20)
                    {
                        N_noisy_ITS_pixels++;
                        vec_ITS_id_row_column_noisy[i_id][i_row].push_back(i_col);

                        NITS_id = i_id;
                        NITS_row = i_row;
                        NITS_col = i_col;
                        output_tree.Fill();
                    }
                }
            }
        }
        printf("Noisy ITS pixels: %d \n", N_noisy_ITS_pixels);
    }

    Int_t nentries_noisy_ITS = (Int_t)input_tree->GetEntries();
    printf("Noisy ITS pixels: %d \n", nentries_noisy_ITS);
    for (Int_t ientry = 0; ientry < nentries_noisy_ITS; ientry++)
    {
        input_tree->GetEntry(ientry);
        vec_ITS_id_row_column_noisy[(Int_t)NITS_in_id][(Int_t)NITS_in_row].push_back((Int_t)NITS_in_col);
    }
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    printf("Start time frame loop \n");
    for (Long64_t counter = start_events; counter < N_events; counter++)
    {
        if (counter != 0 && counter % 10 == 0)
            cout << "." << flush;
        if (counter != 0 && counter % 100 == 0)
        {
            if ((file_entries - 0) > 0)
            {
                Double_t event_percent = 100.0 * ((Double_t)(counter - 0)) / ((Double_t)(file_entries - 0));
                cout << " " << counter << " (" << event_percent << "%) "
                     << "\n"
                     << "==> Processing data " << flush;
            }
        }

        if (!input_chain->GetEntry(counter)) // take the event -> information is stored in event
            break;                           // end of data chunk

        //--------------------
        // Event information
        UShort_t NTPCTracks = TrackHitEvent->getNumTPCTrack();
        Int_t NTPCCluster = TrackHitEvent->getNumTPCCluster();
        Int_t NITSHit = TrackHitEvent->getNumITSHit();
        // printf("NTPCTracks: %d, NTPCCluster: %d, NITSHit: %d \n",NTPCTracks,NTPCCluster,NITSHit);

        // TPC track loop
        vector<TVector3> vec_TV3_dir;
        vector<TVector3> vec_TV3_base;
        vector<Int_t> vec_N_TPC_clusters;
        for (Int_t i_track = 0; i_track < NTPCTracks; i_track++)
        {
            int current_elem = counter * NTPCTracks + i_track;
            TPCTrack = TrackHitEvent->getTPCTrack(i_track);
            Int_t track_id = TPCTrack->get_track_id();
            Float_t track_X = TPCTrack->get_X();
            Float_t track_alpha = TPCTrack->get_alpha();
            for (Int_t i_par = 0; i_par < 5; i_par++)
            {
                track_p[i_par] = TPCTrack->get_par(i_par);
            }
            Float_t track_time0 = TPCTrack->get_time0();
            set_helix(track_X, track_alpha, track_p, B_field);
            for (Int_t i_par = 0; i_par < 6; i_par++)
            {
                vec_helix_par[i_par] = fHelix[i_par];
                // printf("i_par: %d, par: %4.3f \n",i_par,fHelix[i_par]);
            }

            Double_t pT_track_charge = 0.0;
            TLorentzVector TLV_helix_prim = get_TLV_helix(B_field, pT_track_charge);
            Double_t pt_track = TLV_helix_prim.Pt();
            Double_t eta_track = TLV_helix_prim.Eta();
            Double_t phi_track = TLV_helix_prim.Phi();
            Double_t px_track = TLV_helix_prim.Px();
            Double_t py_track = TLV_helix_prim.Py();
            Double_t pz_track = TLV_helix_prim.Pz();

            Float_t phi_track_at_inner_wall = -999.0;
            evaluate_helix(0.0, track_pos);
            Double_t radius_track = TMath::Sqrt(track_pos[0] * track_pos[0] + track_pos[1] * track_pos[1]);
            Float_t track_path_add = +1.0;
            if (radius_track > 87.225)
                track_path_add = -1.0; // go inwards
            for (Float_t track_path = 0.0; fabs(track_path) < 350.0; track_path += track_path_add)
            {
                evaluate_helix(track_path, track_pos);
                radius_track = TMath::Sqrt(track_pos[0] * track_pos[0] + track_pos[1] * track_pos[1]);
                if (radius_track < (85.225 + 1.0) && radius_track > (85.225 - 1.0))
                {
                    phi_track_at_inner_wall = TMath::ATan2(track_pos[1], track_pos[0]);
                    break;
                }
            }

            Float_t phi_track_at_inner_wall_sector = -999.0;
            if (phi_track_at_inner_wall > -999.0)
            {
                phi_track_at_inner_wall *= TMath::RadToDeg(); // -180..180
                phi_track_at_inner_wall += 180.0;             // 0..360
                Int_t sector = (Int_t)(phi_track_at_inner_wall / 20.0);
                phi_track_at_inner_wall_sector = phi_track_at_inner_wall - sector * 20.0;
                // printf("phi: %4.3f, sector: %d, phi_sec: %4.3f \n",phi_track_at_inner_wall,sector,phi_track_at_inner_wall_sector);
            }

            Int_t NTPCCluster_track = TPCTrack->getNumTPCCluster();
            tpc_ncls.push_back(NTPCCluster_track);
            event_num.push_back(counter);

            TVector3 TV3_beam_center;
            TV3_beam_center.SetXYZ(0.0, 0.0, 0.0);
            TPCCluster = TPCTrack->getTPCCluster(0);
            Float_t x_cls_A = TPCCluster->get_cluster_x();
            Float_t y_cls_A = TPCCluster->get_cluster_y();
            Float_t t_cls_A = TPCCluster->get_cluster_time();
            TVector3 TV3_base_track;
            TV3_base_track.SetXYZ(x_cls_A, y_cls_A, t_cls_A);
            TVector3 TV3_dir_trackA;
            TV3_dir_trackA.SetXYZ(x_cls_A, y_cls_A, t_cls_A);
            TPCCluster = TPCTrack->getTPCCluster(NTPCCluster_track - 1);
            Float_t x_cls_B = TPCCluster->get_cluster_x();
            Float_t y_cls_B = TPCCluster->get_cluster_y();
            Float_t t_cls_B = TPCCluster->get_cluster_time();
            TVector3 TV3_dir_trackB;
            TV3_dir_trackB.SetXYZ(x_cls_B, y_cls_B, t_cls_B);
            TVector3 TV3_dir_track = TV3_dir_trackB - TV3_dir_trackA;

            // printf("t_cls first: %4.3f, t_cls last: %4.3f \n",t_cls_A,t_cls_B);

            // vec_TV3_dir.push_back(TV3_dir_track);
            // vec_TV3_base.push_back(TV3_base_track);

            // PCA fit for straight line
            crosses_anomalous_sectors.push_back(-1);
            int count_clusters_used = 0;
            TVector3 current_point, proj_sector;

            // ---------------------
            // Dumping information on for DCA vector, fitted to all clusters
            for (int i = 0; i < NTPCCluster_track; i++)
            {
                TPCCluster = TPCTrack->getTPCCluster(i);
                Double_t cluster[3] = {static_cast<Double_t>(TPCCluster->get_cluster_x()), static_cast<Double_t>(TPCCluster->get_cluster_y()), static_cast<Double_t>(TPCCluster->get_cluster_time())};
                PCA->AddRow(static_cast<Double_t *>(cluster));
            }
            PCA->MakePrincipals();
            TMatrixD eigen = *(PCA->GetEigenVectors());
            TVectorD mean = *(PCA->GetMeanValues());
            TVector3 base;
            TVector3 direction;
            direction.SetXYZ(eigen[0][0], eigen[1][0], eigen[2][0]);
            direction = direction.Unit();
            base.SetXYZ(mean[0], mean[1], mean[2]);
            PCA->Clear();
            float proj_sector_angle;

            for (int i = 0; i < NTPCCluster_track; i++)
            {
                // Filling raw information for wCrit
                TPCCluster = TPCTrack->getTPCCluster(i);
                current_point.SetXYZ(TPCCluster->get_cluster_x(), TPCCluster->get_cluster_y(), TPCCluster->get_cluster_time());
                current_point = calculateDCA_vec_StraightToPoint(base, direction, current_point);
                tpc_raw_dca_wCrit->Fill(TPCCluster->get_cluster_x(), TPCCluster->get_cluster_y(), TPCCluster->get_cluster_time(), TPCCluster->get_cluster_sector(), TPCCluster->get_cluster_row(), current_point[0], current_point[1], current_point[2]);
            }
            // ---------------------

            // ---------------------
            // Anomalous sectors - Per sector check
            for (auto sec : anomalous_sectors)
            {
                for (int i = 0; i < NTPCCluster_track; i++)
                {
                    TPCCluster = TPCTrack->getTPCCluster(i);
                    if (sec == TPCCluster->get_cluster_sector())
                    {
                        anomalous_sector = true;
                        break;
                    }
                }
                proj_sector_angle = ((sec % 18) * 20. + 10.) * 2. * TMath::Pi() / 360.;
                proj_sector.SetXYZ(TMath::Cos(proj_sector_angle), TMath::Sin(proj_sector_angle), 0);
                if (anomalous_sector)
                {
                    int count_fitted_points = 0;
                    for (int i = 0; i < NTPCCluster_track; i++)
                    {
                        TPCCluster = TPCTrack->getTPCCluster(i);
                        if (TPCCluster->get_cluster_sector() != sec)
                        {
                            Double_t cluster[3] = {static_cast<Double_t>(TPCCluster->get_cluster_x()), static_cast<Double_t>(TPCCluster->get_cluster_y()), static_cast<Double_t>(TPCCluster->get_cluster_time())};
                            PCA->AddRow(static_cast<Double_t *>(cluster));
                            count_fitted_points++;
                        }
                    }
                    if(count_fitted_points >= 0){
                        PCA->MakePrincipals();
                        eigen = *(PCA->GetEigenVectors());
                        mean = *(PCA->GetMeanValues());
                        direction.SetXYZ(eigen[0][0], eigen[1][0], eigen[2][0]);
                        direction = direction.Unit();
                        base.SetXYZ(mean[0], mean[1], mean[2]);
                        for (int i = 0; i < NTPCCluster_track; i++)
                        {
                            // Filling raw information for wCrit
                            TPCCluster = TPCTrack->getTPCCluster(i);
                            if (TPCCluster->get_cluster_sector() == sec)
                            {
                                current_point.SetXYZ(TPCCluster->get_cluster_x(), TPCCluster->get_cluster_y(), TPCCluster->get_cluster_time());
                                current_point = calculateDCA_vec_StraightToPoint(base, direction, current_point);
                                tpc_raw_dca_woCrit->Fill(TPCCluster->get_cluster_x(), TPCCluster->get_cluster_y(), TPCCluster->get_cluster_time(), TPCCluster->get_cluster_sector(), TPCCluster->get_cluster_row(), current_point[0], current_point[1], current_point[2]);
                                anomalous_TPC_dca3D_histos[sec]->Fill(std::sqrt(std::pow(current_point[0], 2) + std::pow(current_point[1], 2)));
                                anomalous_TPC_dcaRad_histos[sec]->Fill(proj_sector * current_point);
                            }
                        }
                    }
                }
                PCA->Clear();
                anomalous_sector = false;
            }
            // ---------------------


            for (int i = 0; i < NTPCCluster_track; i++)
            {
                TPCCluster = TPCTrack->getTPCCluster(i);
                Double_t cluster[3] = {static_cast<Double_t>(TPCCluster->get_cluster_x()), static_cast<Double_t>(TPCCluster->get_cluster_y()), static_cast<Double_t>(TPCCluster->get_cluster_time())};
                PCA->AddRow(static_cast<Double_t *>(cluster));
                count_clusters_used++;
            }

            PCA->MakePrincipals();
            eigen = *(PCA->GetEigenVectors());
            mean = *(PCA->GetMeanValues());
            direction.SetXYZ(eigen[0][0], eigen[1][0], eigen[2][0]);
            direction = direction.Unit();
            base.SetXYZ(mean[0], mean[1], mean[2]);
            PCA->Clear();

            if (counter == event_plot)
            {
                direction_PCA.push_back(direction);
                base_PCA.push_back(base);
            }

            angle_dist_PCA->Fill(direction.Angle(TV3_dir_track));

            crosses_anomalous_ITS.push_back(straight_intersect_circle(direction, base, 85.225));
            residuals_in_anomalous_sectors.push_back(vector<float>{});

            if (!crosses_anomalous_ITS.back() && crosses_anomalous_sectors.back() > -1 && count_clusters_used > 10)
            {
                anomalous_TPC_sector_events->Fill(counter);
                if (start_events > 0)
                {
                    cout << "----------------" << endl;
                    cout << "Event: " << counter << endl;
                    cout << "Direction vector: [" << direction[0] << ", " << direction[1] << ", " << direction[2] << "]" << endl;
                    cout << "Base vector: [" << base[0] << ", " << base[1] << ", " << base[2] << "]" << endl;
                    cout << "----------------" << endl;
                }
            }

            chi2_pca.push_back(0);
            chi2_pca_xy.push_back(0);
            chi2_pca_time.push_back(0);

            for (int i = 0; i < NTPCCluster_track; i++)
            {
                TPCCluster = TPCTrack->getTPCCluster(i);
                current_point.SetXYZ(TPCCluster->get_cluster_x(), TPCCluster->get_cluster_y(), TPCCluster->get_cluster_time());
                current_point = calculateDCA_vec_StraightToPoint(base, direction, current_point); // Reusing the vector for DCA
                current_point[0] /= sigma_values_tpc[0];
                current_point[1] /= sigma_values_tpc[1];
                current_point[2] /= sigma_values_tpc[2];
                chi2_pca[fill_index] += (TMath::Power(current_point[0], 2) + TMath::Power(current_point[1], 2) + TMath::Power(current_point[2], 2)) / (NTPCCluster_track - 3);
                chi2_pca_xy[fill_index] += (TMath::Power(current_point[0], 2) + TMath::Power(current_point[1], 2)) / (NTPCCluster_track - 2);
                chi2_pca_time[fill_index] += (TMath::Power(current_point[2], 2)) / (NTPCCluster_track - 1);

                // Filling raw information for woCrit
                current_point.SetXYZ(TPCCluster->get_cluster_x(), TPCCluster->get_cluster_y(), TPCCluster->get_cluster_time());
                current_point = calculateDCA_vec_StraightToPoint(base, direction, current_point);
                // tpc_raw_dca_woCrit->Fill(TPCCluster->get_cluster_x(), TPCCluster->get_cluster_y(), TPCCluster->get_cluster_time(), TPCCluster->get_cluster_sector(), TPCCluster->get_cluster_row(), current_point[0], current_point[1], current_point[2]);

                // if (!crosses_anomalous_ITS.back() && crosses_anomalous_sectors.back() > -1 && count_clusters_used > 10) // && crosses_sector8_highRad>10)
                // {
                //     current_point.SetXYZ(TPCCluster->get_cluster_x(), TPCCluster->get_cluster_y(), TPCCluster->get_cluster_time());
                //     // current_point[0] /= sigma_values_tpc[0];
                //     // current_point[1] /= sigma_values_tpc[1];
                //     // current_point[2] /= sigma_values_tpc[2];
                //     current_point = calculateDCA_vec_StraightToPoint(base, direction, current_point); // Reusing the vector for DCA
                //     residuals_in_anomalous_sectors.back().push_back(current_point.Mag());
                //     proj_sector_angle = ((crosses_anomalous_sectors.back() % 18) * 20. - 10.) * 2. * TMath::Pi() / 360.;
                //     proj_sector.SetXYZ(TMath::Cos(proj_sector_angle), TMath::Sin(proj_sector_angle), 0);
                //     for (int an_sec = 0; an_sec < anomalous_sectors.size(); an_sec++)
                //     {
                //         if (crosses_anomalous_sectors.back() == anomalous_sectors[an_sec])
                //         {
                //             anomalous_TPC_dca3D_histos[an_sec]->Fill(std::sqrt(std::pow(current_point[0], 2) + std::pow(current_point[1], 2)));
                //             anomalous_TPC_dcaRad_histos[an_sec]->Fill(proj_sector * current_point);
                //         }
                //     }
                // }
            }
            chi2_TPC_tracks->Fill(chi2_pca[fill_index]);
            chi2_TPC_tracks_xy->Fill(chi2_pca_xy[fill_index]);
            chi2_TPC_tracks_time->Fill(chi2_pca_time[fill_index]);
            fill_index++;
            count_clusters_used = 0;

            // if(saveToTree) data_tree_TPC->Fill();

            // cout << "--------------------------\n";
            // cout << "Number of TPC clusters for this track: " << NTPCCluster_track << "\n" << endl;
            // PCA->Print("MSE");
            // cout << "First PCA: " << eigen[0][0] << " " << eigen[1][0] << " " << eigen[2][0] << endl;
            // cout << "Direction vector: " << TV3_dir_track[0] << " " << TV3_dir_track[1] << " " << TV3_dir_track[2] << endl;
            // cout << "Base vector: " << TV3_base_track[0] << " " << TV3_base_track[1] << " " << TV3_base_track[2] << endl;
            // cout << "--------------------------\n";

            vec_TV3_dir.push_back(direction);
            vec_TV3_base.push_back(base);
            TV3_dir_track = direction;
            TV3_base_track = base;

            vec_N_TPC_clusters.push_back(NTPCCluster_track);
            TVector3 TV3_dca_center = calculateDCA_vec_StraightToPoint(TV3_base_track, TV3_dir_track, TV3_beam_center);
            Double_t radius_xy = TMath::Sqrt(TV3_dca_center.X() * TV3_dca_center.X() + TV3_dca_center.Y() * TV3_dca_center.Y());

            if (radius_xy < 20.0 && fabs(phi_track_at_inner_wall_sector - 10.0) < 7.0 && NTPCCluster_track > 50)
            {
                // printf("TF: %lld, track: %d, phi_track_at_inner_wall: %4.3f, phi_track_at_inner_wall_sector: %4.3f, TV3_dca_center: {%4.3f, %4.3f, %4.3f}, radius_xy: %4.3f, NTPCCluster_track: %d \n",counter,i_track,phi_track_at_inner_wall,phi_track_at_inner_wall_sector,TV3_dca_center.X(),TV3_dca_center.Y(),TV3_dca_center.Z(),radius_xy,NTPCCluster_track);
            }
            if (NTPCCluster_track > 100)
            {
                // printf("TF: %lld, track: %d, phi_track_at_inner_wall: %4.3f, phi_track_at_inner_wall_sector: %4.3f, TV3_dca_center: {%4.3f, %4.3f, %4.3f}, radius_xy: %4.3f, NTPCCluster_track: %d \n",counter,i_track,phi_track_at_inner_wall,phi_track_at_inner_wall_sector,TV3_dca_center.X(),TV3_dca_center.Y(),TV3_dca_center.Z(),radius_xy,NTPCCluster_track);
            }

            // printf("i_track: %d, par: {%4.3f, %4.3f, %4.3f, %4.3f, %4.3f} \n",i_track,vec_helix_par[0],vec_helix_par[1],vec_helix_par[2],vec_helix_par[3],vec_helix_par[4]);
            if (counter == event_plot)
            {
                vec_track_par_helix.push_back(vec_helix_par);

                // clusters attached to track
                Int_t i_point = 0;
                for (Int_t i_cls = 0; i_cls < NTPCCluster_track; i_cls++)
                {
                    TPCCluster = TPCTrack->getTPCCluster(i_cls);
                    Float_t x_cls = TPCCluster->get_cluster_x();
                    Float_t y_cls = TPCCluster->get_cluster_y();
                    Float_t t_cls = TPCCluster->get_cluster_time();
                    Int_t row_cls = TPCCluster->get_cluster_row();
                    Int_t sec_cls = TPCCluster->get_cluster_sector();
                    Int_t track_id_cls = TPCCluster->get_cluster_track_id();
                    tg_cls_y_vs_x->SetPoint(i_point, x_cls, y_cls);
                    tg_cls_x_vs_time->SetPoint(i_point, t_cls, x_cls);
                    tg_cls_y_vs_time->SetPoint(i_point, t_cls, y_cls);
                    i_point++;

                    if (t_cls < time_min)
                        time_min = t_cls;
                    if (t_cls > time_max)
                        time_max = t_cls;
                }
                vec_tg_cls_y_vs_x.push_back((TGraph *)tg_cls_y_vs_x->Clone());
                vec_tg_cls_x_vs_time.push_back((TGraph *)tg_cls_x_vs_time->Clone());
                vec_tg_cls_y_vs_time.push_back((TGraph *)tg_cls_y_vs_time->Clone());
            }
        }

        // TPC cluster loop -> not attached to tracks
        for (Int_t i_cls = 0; i_cls < NTPCCluster; i_cls++)
        {
            TPCCluster = TrackHitEvent->getTPCCluster(i_cls);
            Float_t x_cls = TPCCluster->get_cluster_x();
            Float_t y_cls = TPCCluster->get_cluster_y();
            Float_t t_cls = TPCCluster->get_cluster_time();
            Int_t row_cls = TPCCluster->get_cluster_row();
            Int_t sec_cls = TPCCluster->get_cluster_sector();
            Int_t track_id_cls = TPCCluster->get_cluster_track_id();
            vec_cls_point[0] = x_cls;
            vec_cls_point[1] = y_cls;
            vec_cls_point[2] = t_cls;
            if (counter == event_plot)
            {
                h_cls_y_vs_x->Fill(x_cls, y_cls);
                h_cls_y_vs_x_coarse->Fill(x_cls, y_cls);
                vec_cls_points.push_back(vec_cls_point);
                // printf("TPC i_cls: %d, time_bin: %4.3f \n",i_cls,t_cls);
            }
            // tpc_raw_clusters->Fill(TPCCluster->get_cluster_x(), TPCCluster->get_cluster_y(), TPCCluster->get_cluster_time(), TPCCluster->get_cluster_sector(), TPCCluster->get_cluster_row());
        }

        // ITS loop
        Int_t N_good_ITS_hits = 0;
        Int_t i_point_ITS_match = 0;
        vector<TVector3> vec_TV3_ITS_hits, vec_TV3_ITS_hits_raw_multi, vec_TV3_ITS_hits_raw_unique;
        TVector3 TV3_ITS_hit, ITS_raw_hit;
        for (Int_t i_hit = 0; i_hit < NITSHit; i_hit++)
        {
            ITSHit = TrackHitEvent->getITSHit(i_hit);
            Float_t x_cls = ITSHit->get_cluster_x();
            Float_t y_cls = ITSHit->get_cluster_y();
            Float_t z_cls = ITSHit->get_cluster_z();
            int64_t BC_cls = ITSHit->get_BC();
            UShort_t row = ITSHit->get_row();
            UShort_t col = ITSHit->get_col();
            UShort_t id = ITSHit->get_id();
            Float_t time_bin = 400.0 + ((Float_t)BC_cls) / 8.0; // 400.0 ???
            if (id > max_ITS_chip_id)
                max_ITS_chip_id = id;
            if (row > max_ITS_row)
                max_ITS_row = row;
            if (col > max_ITS_col)
                max_ITS_col = col;

            itterator = std::find(vec_ITS_id_row_column_noisy[id][row].begin(), vec_ITS_id_row_column_noisy[id][row].end(), col);
            if (itterator != vec_ITS_id_row_column_noisy[id][row].end())
            {
                // pixel flagged as noisy
                continue;
            }
            else
            {
                N_good_ITS_hits++;
            }

            vec_cls_point[0] = x_cls;
            vec_cls_point[1] = y_cls;
            vec_cls_point[2] = time_bin;

            TVector3 TV3_ITS_hit;
            TV3_ITS_hit.SetXYZ(x_cls, y_cls, time_bin);
            ITS_raw_hit.SetXYZ(x_cls, y_cls, z_cls);
            vec_TV3_ITS_hits.push_back(TV3_ITS_hit);
            vec_TV3_ITS_hits_raw_multi.push_back(ITS_raw_hit);
            for (Int_t i_track = 0; i_track < (Int_t)vec_TV3_dir.size(); i_track++)
            {
                TVector3 TV3_dca_ITS_hit = calculateDCA_vec_StraightToPoint(vec_TV3_base[i_track], vec_TV3_dir[i_track], TV3_ITS_hit);
                if (fabs(TV3_dca_ITS_hit.Z()) < 100.0 && TV3_dca_ITS_hit.Perp() < 3.0)
                {
                    if (vec_N_TPC_clusters[i_track] > 30)
                    {
                        printf("Match at TF: %lld, track: %d, ITS hit: {%4.3f, %4.3f, %4.3f}, id,row,col: {%d, %d, %d} \n", counter, i_track, x_cls, y_cls, time_bin, id, row, col);
                    }
                    if (counter == event_plot)
                    {
                        tg_ITS_hit_x_vs_y->SetPoint(i_point_ITS_match, x_cls, y_cls);
                        // if(i_track == 5) vec_TV3_ITS_hits.push_back(TV3_ITS_hit);
                        i_point_ITS_match++;
                    }
                }
            }

            if (counter == event_plot)
            {
                h_cls_y_vs_x->Fill(x_cls, y_cls);
                h_cls_y_vs_x_coarse->Fill(x_cls, y_cls);
                vec_cls_points.push_back(vec_cls_point);
                tg_ITS_hit_x_vs_time->SetPoint(N_good_ITS_hits - 1, time_bin, x_cls);
                tg_ITS_hit_y_vs_time->SetPoint(N_good_ITS_hits - 1, time_bin, y_cls);
                // printf("ITS i_hit: %d, time_bin: %4.3f, row/col/chipid: {%d, %d, %d} \n",i_hit,time_bin,row,col,id);
            }
        }

        num_ITS_hits.push_back(N_good_ITS_hits);

        /// Making the ITS hits unique
        vec_TV3_ITS_hits_raw_unique = unique_ITS_hit(vec_TV3_ITS_hits_raw_multi);

        /// Defining some variables
        TVector3 current_its_cluster;
        vector<TVector3> temp_its_clusters;
        vector<float> residuals_2Point_ITS, residual_DCA_XYZ_ITS, residuals_outliers_ITS;
        float max_distance = 0, its_residual_dca = 0;
        int id_hit1 = 0, id_hit2 = 0, counter_hit1 = 0, counter_hit2 = 0;
        bool print_event = true;
        Int_t found_outlier = 0;
        TVector3 ITS_outer_base, ITS_outer_dir;

        /// CHEKING ITS HITS ///
        if ((Int_t)vec_TV3_ITS_hits_raw_unique.size() > 4)
        {

            // Calculate mutual distance and find two most distant points
            for (auto hit1 : vec_TV3_ITS_hits_raw_unique)
            {
                for (auto hit2 : vec_TV3_ITS_hits_raw_unique)
                {
                    if (max_distance < (hit1 - hit2).Mag())
                    {
                        max_distance = (hit1 - hit2).Mag();
                        id_hit1 = counter_hit1;
                        id_hit2 = counter_hit2;
                    }
                    counter_hit2++;
                }
                counter_hit2 = 0;
                counter_hit1++;
            }
            ITS_outer_dir.SetXYZ(vec_TV3_ITS_hits_raw_unique[id_hit1][0] - vec_TV3_ITS_hits_raw_unique[id_hit2][0],
                                 vec_TV3_ITS_hits_raw_unique[id_hit1][1] - vec_TV3_ITS_hits_raw_unique[id_hit2][1],
                                 vec_TV3_ITS_hits_raw_unique[id_hit1][2] - vec_TV3_ITS_hits_raw_unique[id_hit2][2]);
            ITS_outer_base.SetXYZ(vec_TV3_ITS_hits_raw_unique[id_hit1][0],
                                  vec_TV3_ITS_hits_raw_unique[id_hit1][1],
                                  vec_TV3_ITS_hits_raw_unique[id_hit1][2]);

            vector<float> event_dcas, chi2_values;
            for (Int_t hits_ITS = 0; hits_ITS < (Int_t)vec_TV3_ITS_hits_raw_unique.size(); hits_ITS++)
            {

                /// Calculate straight line between outer most points and then DCA's ///

                current_its_cluster = vec_TV3_ITS_hits_raw_unique[hits_ITS];
                temp_its_clusters.push_back(current_its_cluster);
                its_residual_dca = calculateMinimumDistanceStraightToPoint(ITS_outer_base, ITS_outer_dir, current_its_cluster);
                residual_DCA_XYZ_ITS.push_back(TMath::Abs(calculateDCA_vec_StraightToPoint(ITS_outer_base, ITS_outer_dir, current_its_cluster)[0]));
                residual_DCA_XYZ_ITS.push_back(TMath::Abs(calculateDCA_vec_StraightToPoint(ITS_outer_base, ITS_outer_dir, current_its_cluster)[1]));
                residual_DCA_XYZ_ITS.push_back(TMath::Abs(calculateDCA_vec_StraightToPoint(ITS_outer_base, ITS_outer_dir, current_its_cluster)[2]));
                residuals_2Point_ITS.push_back(its_residual_dca);
                residuals_2Point_ITS_histo->Fill(its_residual_dca);

                //////////////////////////////////////////////////////////////////////////////////////////

                /// For n ITS hits calculate PCA between n-1 hits and DCA to last point ///

                for (Int_t hits_ITS_except = 0; hits_ITS_except < (Int_t)vec_TV3_ITS_hits_raw_unique.size(); hits_ITS_except++)
                {
                    if (hits_ITS_except != hits_ITS)
                    {
                        Double_t cluster[3] = {static_cast<Double_t>(vec_TV3_ITS_hits_raw_unique[hits_ITS_except][0]), static_cast<Double_t>(vec_TV3_ITS_hits_raw_unique[hits_ITS_except][1]), static_cast<Double_t>(vec_TV3_ITS_hits_raw_unique[hits_ITS_except][2])};
                        PCA->AddRow(static_cast<Double_t *>(cluster));
                        // cout << "Filling matrix with: " << vec_TV3_ITS_hits_raw_unique[hits_ITS_except][0] << " " << vec_TV3_ITS_hits_raw_unique[hits_ITS_except][1] << " " << vec_TV3_ITS_hits_raw_unique[hits_ITS_except][2] << endl;
                    }
                }

                float dca_x, dca_y, dca_z;
                TVector3 base, direction, current_dca;

                PCA->MakePrincipals();
                TMatrixD eigen = *(PCA->GetEigenVectors());
                TVectorD mean = *(PCA->GetMeanValues());
                direction.SetXYZ(eigen[0][0], eigen[1][0], eigen[2][0]);
                base.SetXYZ(mean[0], mean[1], mean[2]);
                PCA->Clear();

                current_dca = calculateDCA_vec_StraightToPoint(base, direction, vec_TV3_ITS_hits_raw_unique[hits_ITS]);
                dca_x = TMath::Abs(current_dca[0]);
                dca_y = TMath::Abs(current_dca[1]);
                dca_z = TMath::Abs(current_dca[2]);
                residuals_PCAfit_ITS_histo->Fill(current_dca.Mag());

                if (dca_x < 0.1)
                    residuals_ITS_dx->Fill(vec_TV3_ITS_hits_raw_unique[hits_ITS][0], vec_TV3_ITS_hits_raw_unique[hits_ITS][1], dca_x, dca_x);
                if (dca_y < 0.1)
                    residuals_ITS_dy->Fill(vec_TV3_ITS_hits_raw_unique[hits_ITS][0], vec_TV3_ITS_hits_raw_unique[hits_ITS][1], dca_y, dca_y);
                if (dca_z < 0.1)
                    residuals_ITS_dz->Fill(vec_TV3_ITS_hits_raw_unique[hits_ITS][0], vec_TV3_ITS_hits_raw_unique[hits_ITS][1], dca_z, dca_z);
                // if(dca_r < 0.1) residuals_ITS_dr->Fill(vec_TV3_ITS_hits_raw_unique[hits_ITS][0], vec_TV3_ITS_hits_raw_unique[hits_ITS][1], dca_r, dca_r);

                event_dcas.push_back(current_dca.Mag());
                event_dcas.push_back(dca_x);
                event_dcas.push_back(dca_y);
                event_dcas.push_back(dca_z);

                residuals_PCAfit_ITS_histo->Fill(current_dca.Mag());

                float chi2_pca = 0;
                for (Int_t hits_ITS_except = 0; hits_ITS_except < (Int_t)vec_TV3_ITS_hits_raw_unique.size(); hits_ITS_except++)
                {
                    if (hits_ITS_except != hits_ITS)
                    {
                        TVector3 dca_current_point = calculateDCA_vec_StraightToPoint(base, direction, vec_TV3_ITS_hits_raw_unique[hits_ITS_except]);
                        dca_current_point[0] /= sigma_values_its[0];
                        dca_current_point[1] /= sigma_values_its[1];
                        dca_current_point[2] /= sigma_values_its[2];
                        chi2_pca += TMath::Power(dca_current_point[0], 2) + TMath::Power(dca_current_point[1], 2) + TMath::Power(dca_current_point[2], 2);
                    }
                }
                chi2_PCAfit_ITS_histo->Fill(chi2_pca / ((Int_t)vec_TV3_ITS_hits_raw_unique.size() - 1 - 3));
                chi2_values.push_back(chi2_pca / ((Int_t)vec_TV3_ITS_hits_raw_unique.size() - 1 - 3));
                // if(saveToTree) data_tree_ITS->Fill();

                residuals_ITS_full->Fill(vec_TV3_ITS_hits_raw_unique[hits_ITS][0], vec_TV3_ITS_hits_raw_unique[hits_ITS][1], vec_TV3_ITS_hits_raw_unique[hits_ITS][2], current_dca[0], current_dca[1], current_dca[2], chi2_pca, (Int_t)vec_TV3_ITS_hits_raw_unique.size());

                // if(current_dca.Mag()>0.03 && current_dca.Mag()<0.05){
                if (true)
                {
                    found_outlier = 1;
                }
            }

            if (print_event)
            {
                cout << "Straight-line fit to outer-most points and DCA's to all hits:" << endl;
                Int_t count_loop = 0;
                for (Int_t hits_ITS = 0; hits_ITS < (Int_t)vec_TV3_ITS_hits_raw_unique.size(); hits_ITS++)
                {
                    cout << "DCA_total = " << residuals_2Point_ITS[count_loop] << "; DCA_X = " << residual_DCA_XYZ_ITS[3 * count_loop + 0] << "; DCA_Y = " << residual_DCA_XYZ_ITS[3 * count_loop + 1] << "; DCA_Z = " << residual_DCA_XYZ_ITS[3 * count_loop + 2] << "; X = " << temp_its_clusters[count_loop][0] << "; Y = " << temp_its_clusters[count_loop][1] << "; Z = " << temp_its_clusters[count_loop][2] << endl;
                    count_loop++;
                }
                cout << "Total clusters ITS: " << (Int_t)vec_TV3_ITS_hits_raw_unique.size() << endl;

                cout << "PCA fit to n-1 points and DCA to left-out point" << endl;
                for (int lc = 0; lc < (Int_t)vec_TV3_ITS_hits_raw_unique.size(); lc++)
                {
                    cout << "Chi2 = " << chi2_values[lc] << "; DCA_total = " << event_dcas[4 * lc + 0] << "; DCA_X = " << event_dcas[4 * lc + 1] << "; DCA_Y = " << event_dcas[4 * lc + 2] << "; DCA_Z = " << event_dcas[4 * lc + 3] << "; X = " << vec_TV3_ITS_hits_raw_unique[lc][0] << "; Y = " << vec_TV3_ITS_hits_raw_unique[lc][1] << "; Z = " << vec_TV3_ITS_hits_raw_unique[lc][2] << endl;
                }
                cout << "Total clusters ITS: " << (Int_t)vec_TV3_ITS_hits_raw_unique.size() << endl;
            }

            //////////////////////////////////////////////////////////////////////////////////////////
        }
        else
        {
            residuals_2Point_ITS.push_back(-1);
        }
        residuals_outliers_ITS.push_back(found_outlier);

#if 0
        TVector3 TV3_base = vec_TV3_ITS_hits[6];
        TVector3 TV3_dir  = vec_TV3_ITS_hits[6] - vec_TV3_ITS_hits[7];
        for(Int_t i_hit = 0; i_hit < (Int_t)vec_TV3_ITS_hits.size(); i_hit++)
        {
            TVector3 TV3_dca_ITS_hit = calculateDCA_vec_StraightToPoint(TV3_base,TV3_dir,vec_TV3_ITS_hits[i_hit]);
            TV3_dca_ITS_hit.Print();
        }
#endif

        // To be continued for ITS residuals
        vector<Int_t> vec_ITS_hits_used;
        vec_ITS_hits_used.resize((Int_t)vec_TV3_ITS_hits.size());
        for (Int_t i_hit = 0; i_hit < (Int_t)vec_TV3_ITS_hits.size(); i_hit++)
        {
            vec_ITS_hits_used[i_hit] = 0;
        }
        for (Int_t i_hit = 0; i_hit < (Int_t)vec_TV3_ITS_hits.size() - 1; i_hit++)
        {
            if (vec_ITS_hits_used[i_hit])
                continue;
            vector<TVector3> vec_TV3_ITS_hits_same_time;
            Double_t time_hit = vec_TV3_ITS_hits[i_hit].Z();
            vec_TV3_ITS_hits_same_time.push_back(vec_TV3_ITS_hits[i_hit]);
            vec_ITS_hits_used[i_hit] = 1;
            for (Int_t i_hitB = (i_hit + 1); i_hitB < (Int_t)vec_TV3_ITS_hits.size(); i_hitB++)
            {
                if (vec_ITS_hits_used[i_hitB])
                    continue;
                Double_t time_hitB = vec_TV3_ITS_hits[i_hitB].Z();
                if (time_hit != time_hitB)
                    continue;
                vec_TV3_ITS_hits_same_time.push_back(vec_TV3_ITS_hits[i_hitB]);
                vec_ITS_hits_used[i_hitB] = 1;
            }

            Int_t N_hits_same_time = (Int_t)vec_TV3_ITS_hits_same_time.size();
            Int_t N_hits_same_time_close = 0;
            if (N_hits_same_time > 2)
            {
                TVector3 TV3_base = vec_TV3_ITS_hits_same_time[0];
                TVector3 TV3_dir = vec_TV3_ITS_hits_same_time[0] - vec_TV3_ITS_hits_same_time[1];
                for (Int_t i_hitC = 2; i_hitC < N_hits_same_time; i_hitC++)
                {
                    TVector3 TV3_dca_ITS_hit = calculateDCA_vec_StraightToPoint(TV3_base, TV3_dir, vec_TV3_ITS_hits_same_time[i_hitC]);
                    if (TV3_dca_ITS_hit.Perp() < 0.3)
                    {
                        N_hits_same_time_close++;
                        // printf("TF: %lld \n",counter);
                        // TV3_dca_ITS_hit.Print();
                    }
                }
            }
            if (N_hits_same_time_close > 3)
                printf("TF: %lld, N_hits_same_time_close: %d \n", counter, N_hits_same_time_close);
        }
        //--------------------

        //--------------------
        // Check for TPC to TPC track matches
        Int_t flag_good_TPC_match = 0;
        for (Int_t i_trackA = 0; i_trackA < (Int_t)vec_TV3_dir.size(); i_trackA++)
        {
            for (Int_t i_trackB = (i_trackA + 1); i_trackB < (Int_t)vec_TV3_dir.size(); i_trackB++)
            {
                if (i_trackA == i_trackB)
                    continue;
                TVector3 TV3_dca_TPC_hit = calculateDCA_vec_StraightToPoint(vec_TV3_base[i_trackA], vec_TV3_dir[i_trackA], vec_TV3_base[i_trackB]);
                Float_t dca_3D_lines = calculateDCA_vec_StraightToStraight(vec_TV3_base[i_trackA], vec_TV3_dir[i_trackA], vec_TV3_base[i_trackB], vec_TV3_dir[i_trackB]);
                if (start_events > 0)
                    printf("tracks: {%d, %d}, perp dist: %4.3f, 3D dca: %4.3f \n", i_trackA, i_trackB, TV3_dca_TPC_hit.Perp(), dca_3D_lines);
                // if(fabs(TV3_dca_TPC_hit.Z()) < 6.0 && TV3_dca_TPC_hit.Perp() < 3.0)
                // if(dca_3D_lines < 5.0 && tpc_ncls[i_trackA]>40 && tpc_ncls[i_trackB]>40 && (chi2_pca[i_trackA] < 5 || chi2_pca[i_trackB] < 5) && (num_ITS_hits[i_trackA] >= 5 || num_ITS_hits[i_trackB] >= 5))
                // if((residuals_outliers_ITS[i_trackA] == 0 || residuals_outliers_ITS[i_trackB] == 0))
                if (dca_3D_lines < 5.0 && tpc_ncls[i_trackA] > 40 && tpc_ncls[i_trackB] > 40 && (chi2_pca[i_trackA] < 5 || chi2_pca[i_trackB] < 5) && N_good_ITS_hits > 5)
                {
                    flag_good_TPC_match = 1;
                    printf("Good TPC-to-TPC match found! tracks: {%d, %d} \n", i_trackA, i_trackB);
                    Double_t angle = TMath::ACos((vec_TV3_dir[i_trackA].Dot(vec_TV3_dir[i_trackB])) / (vec_TV3_dir[i_trackA].Mag() * vec_TV3_dir[i_trackB].Mag()));
                    cout << "Angle: " << angle << endl;
                    // if(angle < 0.1 && angle > 0.){
                    //     outlier_events->Fill(counter);
                    // }
                    if (angle > 0.03 && angle < 0.2)
                    {
                        outlier_events->Fill(counter);
                    }
                    angle_dist->Fill(angle);
                }
            }
        }
        //--------------------

        // printf("Time frame: %lld, N_good_ITS_hits: %d \n",counter,N_good_ITS_hits);
    } // end if time frame loop
    //------------------------------------------------------------------------------

    /// Beautifications

    // angular distribution of TPC tracks
    angle_dist->SetTitle("angle distribution between TPC tracks");
    angle_dist->GetXaxis()->SetTitle("angle distribution [rad]");
    angle_dist->GetYaxis()->SetTitle("counts");

    // All outlier events passing the selection criteria
    outlier_events->SetTitle("Outlier events passing the selection criteria");
    outlier_events->GetXaxis()->SetTitle("Event [#]");
    outlier_events->GetYaxis()->SetTitle("counts");

    // Chi^2 of the TPC track fits (using base and dir vectors; either PCA or 2-point fit)
    chi2_TPC_tracks->SetTitle("#chi^{2}_{red} of the TPC track fits");
    chi2_TPC_tracks->GetXaxis()->SetTitle("#chi^{2}_{red}");
    chi2_TPC_tracks->GetYaxis()->SetTitle("counts");

    // Chi^2 of the TPC track fits (using base and dir vectors; either PCA or 2-point fit), only for xy coordinates
    chi2_TPC_tracks_xy->SetTitle("#chi^{2}_{red} of the TPC track fits; only xy coordinates");
    chi2_TPC_tracks_xy->GetXaxis()->SetTitle("#chi^{2}_{red} (xy)");
    chi2_TPC_tracks_xy->GetYaxis()->SetTitle("counts");

    // Chi^2 of the TPC track fits (using base and dir vectors; either PCA or 2-point fit), only for time coordinate
    chi2_TPC_tracks_time->SetTitle("#chi^{2}_{red} of the TPC track fits; only time coordinate");
    chi2_TPC_tracks_time->GetXaxis()->SetTitle("#chi^{2}_{red} (time)");
    chi2_TPC_tracks_time->GetYaxis()->SetTitle("counts");

    // Residual DCA - Track to ITS cluster, 2-point fit
    residuals_2Point_ITS_histo->SetTitle("Residual DCA - Track to ITS cluster, 2 point");
    residuals_2Point_ITS_histo->GetXaxis()->SetTitle("DCA [cm]");
    residuals_2Point_ITS_histo->GetYaxis()->SetTitle("counts");

    // Residual DCA - Track to ITS cluster, PCA
    residuals_PCAfit_ITS_histo->SetTitle("Residual DCA - Track to ITS cluster, PCA");
    residuals_PCAfit_ITS_histo->GetXaxis()->SetTitle("DCA [cm]");
    residuals_PCAfit_ITS_histo->GetYaxis()->SetTitle("counts");

    // Chi2 - Track to ITS cluster
    chi2_PCAfit_ITS_histo->SetTitle("#chi^{2}_{red} - Track to ITS cluster, PCA");
    chi2_PCAfit_ITS_histo->GetXaxis()->SetTitle("#chi^{2}_{red}");
    chi2_PCAfit_ITS_histo->GetYaxis()->SetTitle("counts");

    printf("time: {%4.3f, %4.3f}, max_ITS_chip_id: %d, max_ITS_row: %d, max_ITS_col: %d, N_pixels_hit: %lld, N_noisy_ITS_pixels: %d \n", time_min, time_max, max_ITS_chip_id, max_ITS_row, max_ITS_col, N_pixels_hit, N_noisy_ITS_pixels);
    TH2D *h_cls_x_vs_time = new TH2D("h_cls_x_vs_time", "h_cls_x_vs_time", 1000, time_min, time_max, 1000, -250, 250);
    TH2D *h_cls_y_vs_time = new TH2D("h_cls_y_vs_time", "h_cls_y_vs_time", 1000, time_min, time_max, 1000, -250, 250);

    for (Int_t i_cls = 0; i_cls < (Int_t)vec_cls_points.size(); i_cls++)
    {
        h_cls_x_vs_time->Fill(vec_cls_points[i_cls][2], vec_cls_points[i_cls][0]);
        h_cls_y_vs_time->Fill(vec_cls_points[i_cls][2], vec_cls_points[i_cls][1]);
    }
    //------------------------------------------------------------

    //------------------------------------------------------------
    Draw_1D_histo_and_canvas(angle_dist_PCA, "angular_distribution_PCA", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(angle_dist, "angular_distribution", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(outlier_events, "outlier_events", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(chi2_TPC_tracks, "chi2_TPC_tracks", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(chi2_TPC_tracks_xy, "chi2_TPC_tracks_xy", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(chi2_TPC_tracks_time, "chi2_TPC_tracks_time", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(residuals_2Point_ITS_histo, "residuals_2Point_ITS_histo", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(residuals_PCAfit_ITS_histo, "residuals_PCAfit_ITS_histo", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(chi2_PCAfit_ITS_histo, "chi2_PCAfit_ITS_histo", 720, 720, 720, 720, "");

    Draw_2D_ntuple_and_canvas(residuals_ITS_dx, "residuals_ITS_dx", 1000, 1000, "x:y:dcax");
    Draw_2D_ntuple_and_canvas(residuals_ITS_dy, "residuals_ITS_dy", 1000, 1000, "x:y:dcay");
    Draw_2D_ntuple_and_canvas(residuals_ITS_dz, "residuals_ITS_dz", 1000, 1000, "x:y:dcaz");

    // Anomalous TPC sectors
    Draw_1D_histo_and_canvas(anomalous_TPC_sector_events, "anomalous_TPC_sector_events", 720, 720, 720, 720, "");
    for (int i = 0; i < anomalous_sectors.size(); i++)
    {
        histo_name.Form("residual_dca3D_sector%i_histo", anomalous_sectors[i]);
        Draw_1D_histo_and_canvas(anomalous_TPC_dca3D_histos[i], histo_name, 720, 720, 720, 720, "");
        anomalous_TPC_dca3D_histos[i]->SetTitle(Form("Sector %i: Residual DCA in x and y", anomalous_sectors[i]));
        anomalous_TPC_dca3D_histos[i]->GetXaxis()->SetTitle("sqrt(res_x^2 + res_y^2) (cm)");
        anomalous_TPC_dca3D_histos[i]->GetYaxis()->SetTitle("entries [#]");
        histo_name.Form("residual_dcaRad_sector%i_histo", anomalous_sectors[i]);
        Draw_1D_histo_and_canvas(anomalous_TPC_dcaRad_histos[i], histo_name, 720, 720, 720, 720, "");
        anomalous_TPC_dcaRad_histos[i]->SetTitle(Form("Sector %i: Proj. on sector angle-halfing", anomalous_sectors[i]));
        anomalous_TPC_dcaRad_histos[i]->GetXaxis()->SetTitle("d (cm)");
        anomalous_TPC_dcaRad_histos[i]->GetYaxis()->SetTitle("entries [#]");
    }
    // Draw_1D_histo_and_canvas(residual_dca3D_sector4_histo, "residual_dca3D_sector4_histo", 720, 720, 720, 720, "");
    // Draw_1D_histo_and_canvas(residual_dca3D_sector8_histo, "residual_dca3D_sector8_histo", 720, 720, 720, 720, "");
    // Draw_1D_histo_and_canvas(residual_dca3D_sector17_histo, "residual_dca3D_sector17_histo", 720, 720, 720, 720, "");
    // Draw_1D_histo_and_canvas(residual_dca3D_sector26_histo, "residual_dca3D_sector26_histo", 720, 720, 720, 720, "");
    // Draw_1D_histo_and_canvas(residual_dca3D_sector35_histo, "residual_dca3D_sector35_histo", 720, 720, 720, 720, "");

    // residual_dca3D_sector4_histo->SetTitle("Sector 4: Residual DCA in x and y");
    // residual_dca3D_sector4_histo->GetXaxis()->SetTitle("sqrt(res_x^2 + res_y^2) (cm)");
    // residual_dca3D_sector4_histo->GetYaxis()->SetTitle("entries [#]");

    // residual_dca3D_sector8_histo->SetTitle("Sector 8: Residual DCA in x and y");
    // residual_dca3D_sector8_histo->GetXaxis()->SetTitle("sqrt(res_x^2 + res_y^2) (cm)");
    // residual_dca3D_sector8_histo->GetYaxis()->SetTitle("entries [#]");

    // residual_dca3D_sector17_histo->SetTitle("Sector 17: Residual DCA in x and y");
    // residual_dca3D_sector17_histo->GetXaxis()->SetTitle("sqrt(res_x^2 + res_y^2) (cm)");
    // residual_dca3D_sector17_histo->GetYaxis()->SetTitle("entries [#]");

    // residual_dca3D_sector26_histo->SetTitle("Sector 26: Residual DCA in x and y");
    // residual_dca3D_sector26_histo->GetXaxis()->SetTitle("sqrt(res_x^2 + res_y^2) (cm)");
    // residual_dca3D_sector26_histo->GetYaxis()->SetTitle("entries [#]");

    // residual_dca3D_sector35_histo->SetTitle("Sector 35: Residual DCA in x and y");
    // residual_dca3D_sector35_histo->GetXaxis()->SetTitle("sqrt(res_x^2 + res_y^2) (cm)");
    // residual_dca3D_sector35_histo->GetYaxis()->SetTitle("entries [#]");

    h_cls_y_vs_x->GetXaxis()->SetTitle("x (cm)");
    h_cls_y_vs_x->GetYaxis()->SetTitle("y (cm)");
    h_cls_y_vs_x->GetZaxis()->SetTitle("entries");
    TCanvas *can_h_cls_y_vs_x = Draw_2D_histo_and_canvas((TH2D *)h_cls_y_vs_x, "can_h_cls_y_vs_x", 1100, 800, 0.0, 0.0, "colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_h_cls_y_vs_x->SetLogz(0);
    Draw_Circle_Detector_2D(85.225, 250, 2, 18, kBlack, 1, 2, 0.0, 0.0); // (Float_t radius_in = 1, Float_t radius_out = 2,const Int_t n_radii = 1, const Int_t n_delta_phi = 2, Float_t color = 2, Int_t line_style = 1, Int_t line_width = 1, Float_t x_offset = 0.0, Float_t y_offset = 0.0)

    vector<TPolyLine *> vec_TPL_track_y_vs_x;
    Int_t arr_color[10] = {kRed, kGreen + 2, kMagenta, kCyan, kOrange + 2, kAzure - 2, kYellow + 2, kPink + 2, kViolet + 5, kTeal + 7};
    for (Int_t i_track = 0; i_track < (Int_t)vec_track_par_helix.size(); i_track++)
    {
        // printf("i_track: %d \n",i_track);
        vec_TPL_track_y_vs_x.push_back(new TPolyLine());
        // for (Int_t i_par = 0; i_par < 6; i_par++)
        // {
        //     fHelix[i_par] = vec_track_par_helix[i_track][i_par];
        //     // printf("i_par: %d, par: %4.3f \n",i_par,fHelix[i_par]);
        // }
        // for (Float_t track_path = -85.0; track_path < 350.0; track_path += 1.0)
        // {
        //     evaluate_helix(track_path, track_pos);
        //     Double_t radius_track = TMath::Sqrt(track_pos[0] * track_pos[0] + track_pos[1] * track_pos[1]);
        //     if (radius_track > 250.0)
        //         continue;
        //     vec_TPL_track_y_vs_x[i_track]->SetNextPoint(track_pos[0], track_pos[1]);
        // }
        // if (i_track < 10)
        //     vec_TPL_track_y_vs_x[i_track]->SetLineColorAlpha(arr_color[i_track], 0.2);
        // else
        //     vec_TPL_track_y_vs_x[i_track]->SetLineColorAlpha(kOrange + 2, 0.2);
        // vec_TPL_track_y_vs_x[i_track]->SetLineWidth(3);
        // vec_TPL_track_y_vs_x[i_track]->DrawClone();
        TVector3 track_pos;
        for (Float_t track_path = -1000.0; track_path < 1000.0; track_path += .1)
        {
            track_pos = direction_PCA[i_track] * track_path + base_PCA[i_track];
            // cout << "track_pos: [" << track_pos[0] << ", " << track_pos[1] << ", " << track_pos[2] << "]" << endl;
            if (TMath::Sqrt(track_pos[0] * track_pos[0] + track_pos[1] * track_pos[1]) > 250. || TMath::Sqrt(track_pos[0] * track_pos[0] + track_pos[1] * track_pos[1]) < 70.)
                continue;
            vec_TPL_track_y_vs_x[i_track]->SetNextPoint(track_pos[0], track_pos[1]);
        }
        if (i_track < 10)
            vec_TPL_track_y_vs_x[i_track]->SetLineColorAlpha(arr_color[i_track], 0.2);
        else
            vec_TPL_track_y_vs_x[i_track]->SetLineColorAlpha(kOrange + 2, 0.2);

        vec_TPL_track_y_vs_x[i_track]->SetLineWidth(3);
        vec_TPL_track_y_vs_x[i_track]->DrawClone();

        vec_tg_cls_y_vs_x[i_track]->SetMarkerStyle(24);
        vec_tg_cls_y_vs_x[i_track]->SetMarkerSize(0.5);
        if (i_track < 10)
            vec_tg_cls_y_vs_x[i_track]->SetMarkerColor(arr_color[i_track]);
        else
            vec_tg_cls_y_vs_x[i_track]->SetMarkerColor(kGray + 2);
        vec_tg_cls_y_vs_x[i_track]->Draw("same p");
    }

    tg_ITS_hit_x_vs_y->SetMarkerStyle(20);
    tg_ITS_hit_x_vs_y->SetMarkerSize(0.5);
    tg_ITS_hit_x_vs_y->SetMarkerColor(kRed);
    tg_ITS_hit_x_vs_y->Draw("same p");
    //------------------------------------------------------------

#if 0
    //------------------------------------------------------------
    h_cls_y_vs_x_coarse ->GetXaxis()->SetTitle("x (cm)");
    h_cls_y_vs_x_coarse ->GetYaxis()->SetTitle("y (cm)");
    h_cls_y_vs_x_coarse ->GetZaxis()->SetTitle("entries");
    TCanvas* can_h_cls_y_vs_x_coarse = Draw_2D_histo_and_canvas((TH2D*)h_cls_y_vs_x_coarse,"can_h_cls_y_vs_x_coarse",1100,800,0.0,0.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_h_cls_y_vs_x_coarse->SetLogz(0);
    Draw_Circle_Detector_2D(85.225,250,2,18,kBlack,1,2,0.0,0.0); // (Float_t radius_in = 1, Float_t radius_out = 2,const Int_t n_radii = 1, const Int_t n_delta_phi = 2, Float_t color = 2, Int_t line_style = 1, Int_t line_width = 1, Float_t x_offset = 0.0, Float_t y_offset = 0.0)

    for(Int_t i_track = 0; i_track < (Int_t)vec_track_par_helix.size(); i_track++)
    {
        vec_TPL_track_y_vs_x[i_track] ->SetLineColorAlpha(kOrange+2,0.2);
        vec_TPL_track_y_vs_x[i_track] ->SetLineWidth(3);
        vec_TPL_track_y_vs_x[i_track] ->DrawClone();

        vec_tg_cls_y_vs_x[i_track] ->SetMarkerStyle(24);
        vec_tg_cls_y_vs_x[i_track] ->SetMarkerSize(0.5);
        if(i_track < 10) vec_tg_cls_y_vs_x[i_track] ->SetMarkerColor(arr_color[i_track]);
        else vec_tg_cls_y_vs_x[i_track] ->SetMarkerColor(kGray+2);
        vec_tg_cls_y_vs_x[i_track] ->Draw("same p");
    }
    //------------------------------------------------------------
#endif

    //------------------------------------------------------------
    h_cls_x_vs_time->GetXaxis()->SetTitle("time (tb)");
    h_cls_x_vs_time->GetYaxis()->SetTitle("x (cm)");
    h_cls_x_vs_time->GetZaxis()->SetTitle("entries");
    TCanvas *can_h_cls_x_vs_time = Draw_2D_histo_and_canvas((TH2D *)h_cls_x_vs_time, "can_h_cls_x_vs_time", 1100, 800, 0.0, 0.0, "colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_h_cls_x_vs_time->SetLogz(0);

    for (Int_t i_track = 0; i_track < (Int_t)vec_track_par_helix.size(); i_track++)
    {
        vec_tg_cls_x_vs_time[i_track]->SetMarkerStyle(24);
        vec_tg_cls_x_vs_time[i_track]->SetMarkerSize(0.5);
        if (i_track < 10)
            vec_tg_cls_x_vs_time[i_track]->SetMarkerColor(arr_color[i_track]);
        else
            vec_tg_cls_x_vs_time[i_track]->SetMarkerColor(kGray + 2);
        vec_tg_cls_x_vs_time[i_track]->Draw("same p");
    }

    tg_ITS_hit_x_vs_time->SetMarkerStyle(20);
    tg_ITS_hit_x_vs_time->SetMarkerSize(0.5);
    tg_ITS_hit_x_vs_time->SetMarkerColor(kRed);
    tg_ITS_hit_x_vs_time->Draw("same p");
    //------------------------------------------------------------

    //------------------------------------------------------------
    h_cls_y_vs_time->GetXaxis()->SetTitle("time (tb)");
    h_cls_y_vs_time->GetYaxis()->SetTitle("y (cm)");
    h_cls_y_vs_time->GetZaxis()->SetTitle("entries");
    TCanvas *can_h_cls_y_vs_time = Draw_2D_histo_and_canvas((TH2D *)h_cls_y_vs_time, "can_h_cls_y_vs_time", 1100, 800, 0.0, 0.0, "colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_h_cls_y_vs_time->SetLogz(0);

    for (Int_t i_track = 0; i_track < (Int_t)vec_track_par_helix.size(); i_track++)
    {
        vec_tg_cls_y_vs_time[i_track]->SetMarkerStyle(24);
        vec_tg_cls_y_vs_time[i_track]->SetMarkerSize(0.5);
        if (i_track < 10)
            vec_tg_cls_y_vs_time[i_track]->SetMarkerColor(arr_color[i_track]);
        else
            vec_tg_cls_y_vs_time[i_track]->SetMarkerColor(kGray + 2);
        vec_tg_cls_y_vs_time[i_track]->Draw("same p");
    }

    tg_ITS_hit_y_vs_time->SetMarkerStyle(20);
    tg_ITS_hit_y_vs_time->SetMarkerSize(0.5);
    tg_ITS_hit_y_vs_time->SetMarkerColor(kRed);
    tg_ITS_hit_y_vs_time->Draw("same p");
    //------------------------------------------------------------

    //------------------------------------------------------------
    if (flag_ITS_noisy)
    {
        file_ITS_noisy_pixels->cd();
        output_tree.Write();
        file_ITS_noisy_pixels->Close();
    }
    //------------------------------------------------------------

    if (saveToTree)
    {
        cout << "Writing to tree..." << endl;
        data_tree_TPC->Write();
        delete data_tree_TPC;
        data_tree_ITS->Write();
        delete data_tree_ITS;
        // data_tree_TPC_raw_all->Write();
        // delete data_tree_TPC_raw_all;
        data_tree_TPC_raw_wCrit->Write();
        delete data_tree_TPC_raw_wCrit;
        data_tree_TPC_raw_woCrit->Write();
        delete data_tree_TPC_raw_woCrit;
        qafileITSTPC->Write();
        qafileITSTPC->Close();
        delete qafileITSTPC;
    }
    delete PCA;
}

void AnaResiduals()
{
    TFile *inputfile = TFile::Open("QA_tree.root");

    TTree *input_tree = (TTree *)inputfile->Get("residuals_ITS_full");
    Float_t x, y, z, dx, dy, dz, chi2, npoints;

    input_tree->SetBranchAddress("x", &x);
    input_tree->SetBranchAddress("y", &y);
    input_tree->SetBranchAddress("z", &z);
    input_tree->SetBranchAddress("dca_x", &dx);
    input_tree->SetBranchAddress("dca_y", &dy);
    input_tree->SetBranchAddress("dca_z", &dz);
    input_tree->SetBranchAddress("chi2", &chi2);
    input_tree->SetBranchAddress("n_points", &npoints);

    Int_t nentries = (Int_t)input_tree->GetEntries();
    printf("nentries: %d \n", nentries);

    //--------------------------------------
    TH2D *h2D_rad_res_vs_radius = new TH2D("h2D_rad_res_vs_radius", "h2D_rad_res_vs_radius", 100, 0, 45, 100, -1.0, 1.0);
    TProfile2D *TP2D_rad_res_vs_xy = new TProfile2D("TP2D_rad_res_vs_xy", "TP2D_rad_res_vs_xy", 100, -45, 45, 100, -45, 45);

    TH2D *h2D_azim_dcax = new TH2D("h2D_azim_dcax", "h2D_azim_dcax", 50, 0, TMath::Pi(), 50, -0.05, 0.05);
    TH2D *h2D_azim_dcay = new TH2D("h2D_azim_dcay", "h2D_azim_dcay", 50, 0, TMath::Pi(), 50, -0.05, 0.05);
    TH2D *h2D_azim_dcaz = new TH2D("h2D_azim_dcaz", "h2D_azim_dcaz", 50, 0, TMath::Pi(), 50, -0.05, 0.05);
    //--------------------------------------

    //--------------------------------------
    TVector3 y_axis;
    y_axis.SetXYZ(0, 1, 0);

    for (Int_t ientry = 0; ientry < nentries; ientry++)
    {
        input_tree->GetEntry(ientry);

        if (fabs(dx) > 1.0 || fabs(dy) > 1.0 || fabs(dz) > 1.0)
            continue;
        if (chi2 > 5)
            continue;

        Float_t radius = TMath::Sqrt(x * x + y * y);

        TVector3 TV3_point;
        TV3_point.SetXYZ(x, y, 0);
        Double_t mag_TV3_point = TV3_point.Mag();
        TV3_point *= 1.0 / TV3_point.Mag();
        TVector3 TV3_res;
        TV3_res.SetXYZ(dx, dy, 0);
        Double_t rad_res = TV3_res.Dot(TV3_point);

        // if(fabs(rad_res) > 5.0)
        //{
        // printf(" \n");
        // printf("rad_res: %4.3f, dx: %4.3f, dy: %4.3f, mag_TV3_point: %4.3f \n",rad_res,dx,dy,mag_TV3_point);
        // TV3_res.Print();
        // TV3_point.Print();
        //}

        h2D_rad_res_vs_radius->Fill(radius, rad_res);
        TP2D_rad_res_vs_xy->Fill(x, y, rad_res);

        h2D_azim_dcax->Fill(TMath::ACos(y_axis.Dot(TV3_point) / TV3_point.Mag()), dx);
        h2D_azim_dcay->Fill(TMath::ACos(y_axis.Dot(TV3_point) / TV3_point.Mag()), dy);
        h2D_azim_dcaz->Fill(TMath::ACos(y_axis.Dot(TV3_point) / TV3_point.Mag()), dz);
    }
    //--------------------------------------

    //--------------------------------------
    h2D_rad_res_vs_radius->GetXaxis()->SetTitle("radius (cm)");
    h2D_rad_res_vs_radius->GetYaxis()->SetTitle("dr (cm)");
    h2D_rad_res_vs_radius->GetZaxis()->SetTitle("entries");
    TProfile *TP_rad_res_vs_radius = (TProfile *)h2D_rad_res_vs_radius->ProfileX("TP_rad_res_vs_radius", 1, -1);
    TCanvas *can_h2D_rad_res_vs_radius = Draw_2D_histo_and_canvas(h2D_rad_res_vs_radius, "can_h2D_rad_res_vs_radius", 800, 600, 0.0, 0.0, "colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_h2D_rad_res_vs_radius->SetLogz(0);
    TP_rad_res_vs_radius->SetLineColor(kRed);
    TP_rad_res_vs_radius->SetLineWidth(3);
    TP_rad_res_vs_radius->DrawCopy("same hist");
    //--------------------------------------

    //--------------------------------------
    TP2D_rad_res_vs_xy->GetXaxis()->SetTitle("x (cm)");
    TP2D_rad_res_vs_xy->GetYaxis()->SetTitle("y (cm)");
    TP2D_rad_res_vs_xy->GetZaxis()->SetTitle("dr (cm)");
    TCanvas *can_TP2D_rad_res_vs_xy = Draw_2D_histo_and_canvas((TH2D *)TP2D_rad_res_vs_xy, "can_TP2D_rad_res_vs_xy", 800, 600, 0.0, 0.0, "colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_TP2D_rad_res_vs_xy->SetLogz(0);
    //--------------------------------------

    //--------------------------------------
    h2D_azim_dcax->GetXaxis()->SetTitle("#phi (rad)");
    h2D_azim_dcax->GetYaxis()->SetTitle("dca_x (cm)");
    h2D_azim_dcax->GetZaxis()->SetTitle("entries");
    TProfile *prof_dcax = h2D_azim_dcax->ProfileX("prof_dcax", 1, -1);
    TCanvas *can_azim_dcax = Draw_2D_histo_and_canvas(h2D_azim_dcax, "can_azim_dcax", 800, 600, 0.0, 0.0, "colz");
    can_azim_dcax->SetLogz(0);
    prof_dcax->SetLineColor(kRed);
    prof_dcax->SetLineWidth(3);
    prof_dcax->DrawCopy("same hist");

    h2D_azim_dcay->GetXaxis()->SetTitle("#phi (rad)");
    h2D_azim_dcay->GetYaxis()->SetTitle("dca_y (cm)");
    h2D_azim_dcay->GetZaxis()->SetTitle("entries");
    TProfile *prof_dcay = h2D_azim_dcaz->ProfileX("prof_dcay", 1, -1);
    TCanvas *can_azim_dcay = Draw_2D_histo_and_canvas(h2D_azim_dcay, "can_azim_dcay", 800, 600, 0.0, 0.0, "colz");
    can_azim_dcay->SetLogz(0);
    prof_dcay->SetLineColor(kRed);
    prof_dcay->SetLineWidth(3);
    prof_dcay->DrawCopy("same hist");

    h2D_azim_dcaz->GetXaxis()->SetTitle("#phi (rad)");
    h2D_azim_dcaz->GetYaxis()->SetTitle("dca_z (cm)");
    h2D_azim_dcaz->GetZaxis()->SetTitle("entries");
    TProfile *prof_dcaz = h2D_azim_dcaz->ProfileX("prof_dcaz", 1, -1);
    TCanvas *can_azim_dcaz = Draw_2D_histo_and_canvas(h2D_azim_dcaz, "can_azim_dcaz", 800, 600, 0.0, 0.0, "colz");
    can_azim_dcaz->SetLogz(0);
    prof_dcaz->SetLineColor(kRed);
    prof_dcaz->SetLineWidth(3);
    prof_dcaz->DrawCopy("same hist");
    //--------------------------------------
}