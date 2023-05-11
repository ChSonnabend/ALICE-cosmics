
//static const Double_t kAlmost1=1. - Double_t(FLT_EPSILON);
//static const Double_t kAlmost0=Double_t(FLT_MIN);
//static const Double_t kB2C=-0.299792458e-3;
//static Float_t fHelix[9];
//static Float_t B_field = -0.001; // -5.0 kG
//static Float_t track_pos[3];

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
vector<float> sort_indices(const vector<T> &v) {

  vector<float> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

void AnaTrackHitEvent(Long64_t N_events = -1, Int_t event_plot = 0, Int_t flag_ITS_noisy = 0)
{
    // To analyze cosmics events, output from dumpClusters.C
    //.L AnaTrackHitEvent.cc++
    // AnaTrackHitEvent(-2,0,0)

    // AnaTrackHitEvent(-2,1,0) -> central membrane, misalignment?
    // AnaTrackHitEvent(-2,23,0) -> central membrane, misalignment?
    // AnaTrackHitEvent(-2,26,0) -> central membrane, misalignment?


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
    TPrincipal* PCA = new TPrincipal(3, "");
    int fill_index = 0;
    vector<int> tpc_ncls, event_num, num_ITS_hits;
    vector<float> chi2_pca, chi2_pca_xy, chi2_pca_time;
    vector<float> sigma_values = {0.22,0.22,0.7}; /// sigma_x = sigma_y = 1mm, sigma_time=2 
    TChain* input_chain = NULL;
    input_chain = new TChain("TrackHitEvent" , "TrackHitEvent");
    //TString addfile = "./Data/Tree_Cluster_points_cosmics_V3.root";
    //TString addfile = "./Data/Tree_Cluster_points_cosmics_V5c.root";
    //TString addfile = "./Data/Tree_Cluster_points_cosmics_V5e.root";
    TString addfile = "./Data/IchFindAlexToll.root";
    input_chain ->AddFile(addfile.Data(),-1,"TrackHitEvent");
    Long64_t file_entries = input_chain->GetEntries();
    printf("entries: %lld \n",file_entries);

    AlTrackHitEvent* TrackHitEvent = new AlTrackHitEvent();
    input_chain  ->SetBranchAddress("Events", &TrackHitEvent);

    AlTrack      *TPCTrack   = new AlTrack();
    AlTPCCluster *TPCCluster = new AlTPCCluster();
    AlITSHit     *ITSHit     = new AlITSHit();
    Float_t  track_p[5];
    vector< vector<Float_t> > vec_track_par_helix;
    vector<Float_t> vec_helix_par;
    vec_helix_par.resize(6);
    //------------------------------------------------------------




    //------------------------------------------------------------
    TFile* file_ITS_noisy_pixels = NULL;
    Float_t NITS_id, NITS_row, NITS_col;
    TTree output_tree("ITS_noisy","a simple Tree with simple variables");
    if(flag_ITS_noisy)
    {
        file_ITS_noisy_pixels = new TFile("file_ITS_noisy_pixels_V2.root","RECREATE");
        output_tree.Branch("id",&NITS_id);
        output_tree.Branch("row",&NITS_row);
        output_tree.Branch("col",&NITS_col);
    }


    TFile* inputfile = TFile::Open("file_ITS_noisy_pixels.root");
    TTree *input_tree = (TTree*)inputfile->Get("ITS_noisy");
    Float_t NITS_in_id, NITS_in_row, NITS_in_col;
    input_tree->SetBranchAddress("id",&NITS_in_id);
    input_tree->SetBranchAddress("row",&NITS_in_row);
    input_tree->SetBranchAddress("col",&NITS_in_col);
    //------------------------------------------------------------




    //------------------------------------------------------------
    TH2D*   h_cls_y_vs_x  = new TH2D("h_cls_y_vs_x","h_cls_y_vs_x",1000,-250,250,1000,-250,250);
    TH2D*   h_cls_y_vs_x_coarse  = new TH2D("h_cls_y_vs_x_coarse","h_cls_y_vs_x_coarse",300,-250,250,300,-250,250);
    TH1D*   angle_dist = new TH1D("angle_dist", "angle_dist", 10000, 0, TMath::Pi());
    TH1D*   outlier_events = new TH1D("outlier_events", "outlier_events", 50000, 0, 50000);
    TH1D*   chi2_tracks = new TH1D("chi2_tracks", "chi2_tracks", 500, 0, 50);
    TH1D*   chi2_tracks_xy = new TH1D("chi2_tracks_xy", "chi2_tracks_xy", 500, 0, 50);
    TH1D*   chi2_tracks_time = new TH1D("chi2_tracks_time", "chi2_tracks_time", 500, 0, 50);
    TH1D*   residuals_DCA_ITS_histo = new TH1D("residuals_DCA_ITS", "residuals_DCA_ITS", 5000, 0, 10);
    TH1D*   residuals_outliers_ITS_histo = new TH1D("residual_outliers_ITS", "residual_outliers_ITS", 5, 0, 1.2);
    TGraph* tg_cls_y_vs_x = new TGraph();
    TGraph* tg_cls_x_vs_time = new TGraph();
    TGraph* tg_cls_y_vs_time = new TGraph();
    TGraph* tg_ITS_hit_x_vs_y = new TGraph();
    TGraph* tg_ITS_hit_x_vs_time = new TGraph();
    TGraph* tg_ITS_hit_y_vs_time = new TGraph();
    vector<TGraph*> vec_tg_cls_y_vs_x;
    vector<TGraph*> vec_tg_cls_x_vs_time;
    vector<TGraph*> vec_tg_cls_y_vs_time;
    //------------------------------------------------------------




    //------------------------------------------------------------
    vec_track_par_helix.clear();
    Long64_t start_events = 0;
    if(N_events == -1) N_events = file_entries;
    if(N_events == -2)
    {
        start_events = event_plot;
        N_events     = event_plot+1;
    }
    Float_t time_min = +9999999.0;
    Float_t time_max = -9999999.0;
    vector< vector<Float_t> > vec_cls_points;
    vector<Float_t> vec_cls_point;
    vec_cls_point.resize(3);
    Int_t max_ITS_chip_id = 0;
    Int_t max_ITS_row     = 0;
    Int_t max_ITS_col     = 0;

    printf("Define ITS pixels vectors \n");
    // 108 144 180 2688 3360 8232 9408 chips for each ITS layer, total 24120
    // 7 layers
    // 24120 chips, 512 rows, 1024 columns

    vector <vector< vector<Int_t> > > vec_ITS_id_row_column;
    vector <vector< vector<Int_t> > > vec_ITS_id_row_column_noisy;
    vector< vector< vector<Int_t> > > vec_ITS_id_row_counter;
    std::vector<Int_t>::iterator itterator;
    vec_ITS_id_row_column.resize(24120);
    vec_ITS_id_row_column_noisy.resize(24120);
    vec_ITS_id_row_counter.resize(24120);
    for(Int_t i_id = 0; i_id < 24120; i_id++)
    {
        vec_ITS_id_row_column[i_id].resize(512);
        vec_ITS_id_row_column_noisy[i_id].resize(512);
        vec_ITS_id_row_counter[i_id].resize(512);
    }

    Long64_t N_pixels_hit = 0;


    //------------------------------------------------------------------------------
    Int_t N_noisy_ITS_pixels = 0;
    if(flag_ITS_noisy)
    {
        printf("Start time frame loop -> determine noisy ITS pixels \n");
        for(Long64_t counter = 0; counter < N_events; counter++)
        {
            if (counter != 0  &&  counter % 10 == 0)
                cout << "." << flush;
            if (counter != 0  &&  counter % 100 == 0)
            {
                if((file_entries-0) > 0)
                {
                    Double_t event_percent = 100.0*((Double_t)(counter-0))/((Double_t)(file_entries-0));
                    cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data " << flush;
                }
            }

            if (!input_chain->GetEntry( counter )) // take the event -> information is stored in event
                break;  // end of data chunk


            //--------------------
            // Event information
            UShort_t NTPCTracks  = TrackHitEvent ->getNumTPCTrack();
            Int_t    NTPCCluster = TrackHitEvent ->getNumTPCCluster();
            Int_t    NITSHit     = TrackHitEvent ->getNumITSHit();
            //printf("NTPCTracks: %d, NTPCCluster: %d, NITSHit: %d \n",NTPCTracks,NTPCCluster,NITSHit);

            // ITS loop
            for(Int_t i_hit = 0; i_hit < NITSHit; i_hit++)
            {
                ITSHit = TrackHitEvent->getITSHit(i_hit);
                Float_t  x_cls        = ITSHit ->get_cluster_x();
                Float_t  y_cls        = ITSHit ->get_cluster_y();
                Float_t  z_cls        = ITSHit ->get_cluster_z();
                int64_t  BC_cls       = ITSHit ->get_BC();
                UShort_t row          = ITSHit ->get_row();
                UShort_t col          = ITSHit ->get_col();
                UShort_t id           = ITSHit ->get_id();
                Float_t  time_bin     = ((Float_t)BC_cls)/8.0;
                if(id > max_ITS_chip_id) max_ITS_chip_id = id;
                if(row > max_ITS_row) max_ITS_row = row;
                if(col > max_ITS_col) max_ITS_col = col;

                itterator = std::find(vec_ITS_id_row_column[id][row].begin(), vec_ITS_id_row_column[id][row].end(), col);
                if(itterator != vec_ITS_id_row_column[id][row].end() )
                {
                    // pixel already in list
                    Int_t index = itterator - vec_ITS_id_row_column[id][row].begin();
                    vec_ITS_id_row_counter[id][row][index]++;
                    //printf("pixel at id: %d, row: %d, col: %d, has %d hits \n",id,row,col,vec_ITS_id_row_counter[id][row][index]);
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


        for(Int_t i_id = 0; i_id < 24120; i_id++) // total number of ITS-2 chips
        {
            for(Int_t i_row = 0; i_row < 512; i_row++)
            {
                Int_t N_columns_hit = (Int_t)vec_ITS_id_row_column[i_id][i_row].size();
                for(Int_t i_hit = 0; i_hit < N_columns_hit; i_hit++)
                {
                    Int_t i_col  = vec_ITS_id_row_column[i_id][i_row][i_hit];
                    Int_t N_hits = vec_ITS_id_row_counter[i_id][i_row][i_hit];
                    //printf("id: %d, row: %d, hit: %d, col: %d, N_hits: %d \n",i_id,i_row,i_hit,i_col,N_hits);
                    if(N_hits > 20)
                    {
                        N_noisy_ITS_pixels++;
                        vec_ITS_id_row_column_noisy[i_id][i_row].push_back(i_col);

                        NITS_id  = i_id;
                        NITS_row = i_row;
                        NITS_col = i_col;
                        output_tree.Fill();
                    }
                }
            }
        }
        printf("Noisy ITS pixels: %d \n",N_noisy_ITS_pixels);
    }

    Int_t nentries_noisy_ITS = (Int_t)input_tree->GetEntries();
    printf("Noisy ITS pixels: %d \n",nentries_noisy_ITS);
    for(Int_t ientry = 0; ientry < nentries_noisy_ITS; ientry++)
    {
        input_tree->GetEntry(ientry);
        vec_ITS_id_row_column_noisy[(Int_t)NITS_in_id][(Int_t)NITS_in_row].push_back((Int_t)NITS_in_col);
    }
    //------------------------------------------------------------------------------




    //------------------------------------------------------------------------------
    printf("Start time frame loop \n");
    for(Long64_t counter = start_events; counter < N_events; counter++)
    {
        if (counter != 0  &&  counter % 10 == 0)
            cout << "." << flush;
        if (counter != 0  &&  counter % 100 == 0)
        {
            if((file_entries-0) > 0)
            {
                Double_t event_percent = 100.0*((Double_t)(counter-0))/((Double_t)(file_entries-0));
                cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data " << flush;
            }
        }

        if (!input_chain->GetEntry( counter )) // take the event -> information is stored in event
            break;  // end of data chunk


        //--------------------
        // Event information
        UShort_t NTPCTracks  = TrackHitEvent ->getNumTPCTrack();
        Int_t    NTPCCluster = TrackHitEvent ->getNumTPCCluster();
        Int_t    NITSHit     = TrackHitEvent ->getNumITSHit();
        //printf("NTPCTracks: %d, NTPCCluster: %d, NITSHit: %d \n",NTPCTracks,NTPCCluster,NITSHit);

        // TPC track loop
        vector<TVector3> vec_TV3_dir;
        vector<TVector3> vec_TV3_base;
        vector<Int_t>    vec_N_TPC_clusters;
        for(Int_t i_track = 0; i_track < NTPCTracks; i_track++)
        {
            int current_elem = counter*NTPCTracks + i_track;
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

            Int_t    NTPCCluster_track =  TPCTrack->getNumTPCCluster();
            tpc_ncls.push_back(NTPCCluster_track);
            event_num.push_back(counter);

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

            //printf("t_cls first: %4.3f, t_cls last: %4.3f \n",t_cls_A,t_cls_B);

            // vec_TV3_dir.push_back(TV3_dir_track);
            // vec_TV3_base.push_back(TV3_base_track);

            // PCA fit for straight line
            for(int i=0; i<NTPCCluster_track; i++){
                TPCCluster = TPCTrack->getTPCCluster(i);
                // cout << "Cluster " << i << ":" << TPCCluster->get_cluster_x() << " " << TPCCluster->get_cluster_y() << " " << TPCCluster->get_cluster_time() << endl;
                Double_t cluster[3] = {static_cast<Double_t>(TPCCluster->get_cluster_x()), static_cast<Double_t>(TPCCluster ->get_cluster_y()), static_cast<Double_t>(TPCCluster ->get_cluster_time())};
                PCA->AddRow(static_cast<Double_t*>(cluster));
            }
            PCA->MakePrincipals();
            TMatrixD eigen = *(PCA->GetEigenVectors());
            TVectorD mean = *(PCA->GetMeanValues());
            TVector3 base; TVector3 direction;
            direction.SetXYZ(eigen[0][0],eigen[1][0],eigen[2][0]);
            base.SetXYZ(mean[0],mean[1],mean[2]);

            chi2_pca.push_back(0);
            chi2_pca_xy.push_back(0);
            chi2_pca_time.push_back(0);
            TVector3 current_point;
            for(int i=0; i<NTPCCluster_track; i++){
                TPCCluster = TPCTrack->getTPCCluster(i);
                current_point.SetXYZ(TPCCluster->get_cluster_x(), TPCCluster ->get_cluster_y(), TPCCluster ->get_cluster_time());
                current_point = calculateDCA_vec_StraightToPoint(TV3_base_track,TV3_dir_track,current_point);
                current_point[0] /= sigma_values[0];
                current_point[1] /= sigma_values[1];
                current_point[2] /= sigma_values[2];
                chi2_pca[fill_index] += (TMath::Power(current_point[0], 2) + TMath::Power(current_point[1], 2) + TMath::Power(current_point[2], 2))/(NTPCCluster_track-3);
                chi2_pca_xy[fill_index] += (TMath::Power(current_point[0], 2) + TMath::Power(current_point[1], 2))/(NTPCCluster_track-2);
                chi2_pca_time[fill_index] += (TMath::Power(current_point[2], 2))/(NTPCCluster_track-1);
            }
            chi2_tracks->Fill(chi2_pca[fill_index]);
            chi2_tracks_xy->Fill(chi2_pca_xy[fill_index]);
            chi2_tracks_time->Fill(chi2_pca_time[fill_index]);
            fill_index++;


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
            TVector3 TV3_dca_center = calculateDCA_vec_StraightToPoint(TV3_base_track,TV3_dir_track,TV3_beam_center);
            Double_t radius_xy = TMath::Sqrt(TV3_dca_center.X()*TV3_dca_center.X() + TV3_dca_center.Y()*TV3_dca_center.Y());

            if(radius_xy < 20.0 && fabs(phi_track_at_inner_wall_sector - 10.0) < 7.0 && NTPCCluster_track > 50)
            {
                //printf("TF: %lld, track: %d, phi_track_at_inner_wall: %4.3f, phi_track_at_inner_wall_sector: %4.3f, TV3_dca_center: {%4.3f, %4.3f, %4.3f}, radius_xy: %4.3f, NTPCCluster_track: %d \n",counter,i_track,phi_track_at_inner_wall,phi_track_at_inner_wall_sector,TV3_dca_center.X(),TV3_dca_center.Y(),TV3_dca_center.Z(),radius_xy,NTPCCluster_track);
            }
            if(NTPCCluster_track > 100)
            {
                //printf("TF: %lld, track: %d, phi_track_at_inner_wall: %4.3f, phi_track_at_inner_wall_sector: %4.3f, TV3_dca_center: {%4.3f, %4.3f, %4.3f}, radius_xy: %4.3f, NTPCCluster_track: %d \n",counter,i_track,phi_track_at_inner_wall,phi_track_at_inner_wall_sector,TV3_dca_center.X(),TV3_dca_center.Y(),TV3_dca_center.Z(),radius_xy,NTPCCluster_track);
            }

            //printf("i_track: %d, par: {%4.3f, %4.3f, %4.3f, %4.3f, %4.3f} \n",i_track,vec_helix_par[0],vec_helix_par[1],vec_helix_par[2],vec_helix_par[3],vec_helix_par[4]);
            if(counter == event_plot)
            {
                vec_track_par_helix.push_back(vec_helix_par);

                // clusters attached to track
                Int_t i_point = 0;
                for(Int_t i_cls = 0; i_cls < NTPCCluster_track; i_cls++)
                {
                    TPCCluster = TPCTrack->getTPCCluster(i_cls);
                    Float_t x_cls        = TPCCluster ->get_cluster_x();
                    Float_t y_cls        = TPCCluster ->get_cluster_y();
                    Float_t t_cls        = TPCCluster ->get_cluster_time();
                    Int_t   row_cls      = TPCCluster ->get_cluster_row();
                    Int_t   sec_cls      = TPCCluster ->get_cluster_sector();
                    Int_t   track_id_cls = TPCCluster ->get_cluster_track_id();
                    tg_cls_y_vs_x    ->SetPoint(i_point,x_cls,y_cls);
                    tg_cls_x_vs_time ->SetPoint(i_point,t_cls,x_cls);
                    tg_cls_y_vs_time ->SetPoint(i_point,t_cls,y_cls);
                    i_point++;

                    if(t_cls < time_min) time_min = t_cls;
                    if(t_cls > time_max) time_max = t_cls;
                }
                vec_tg_cls_y_vs_x.push_back((TGraph*)tg_cls_y_vs_x->Clone());
                vec_tg_cls_x_vs_time.push_back((TGraph*)tg_cls_x_vs_time->Clone());
                vec_tg_cls_y_vs_time.push_back((TGraph*)tg_cls_y_vs_time->Clone());
            }

            // PCA->MakeHistograms();
            PCA->Clear();

        }

        // TPC cluster loop -> not attached to tracks
        for(Int_t i_cls = 0; i_cls < NTPCCluster; i_cls++)
        {
            TPCCluster = TrackHitEvent->getTPCCluster(i_cls);
            Float_t x_cls        = TPCCluster ->get_cluster_x();
            Float_t y_cls        = TPCCluster ->get_cluster_y();
            Float_t t_cls        = TPCCluster ->get_cluster_time();
            Int_t   row_cls      = TPCCluster ->get_cluster_row();
            Int_t   sec_cls      = TPCCluster ->get_cluster_sector();
            Int_t   track_id_cls = TPCCluster ->get_cluster_track_id();
            vec_cls_point[0] = x_cls;
            vec_cls_point[1] = y_cls;
            vec_cls_point[2] = t_cls;
            if(counter == event_plot)
            {
                h_cls_y_vs_x ->Fill(x_cls,y_cls);
                h_cls_y_vs_x_coarse ->Fill(x_cls,y_cls);
                vec_cls_points.push_back(vec_cls_point);
                //printf("TPC i_cls: %d, time_bin: %4.3f \n",i_cls,t_cls);
            }
        }

        // ITS loop
        Int_t N_good_ITS_hits   = 0;
        Int_t i_point_ITS_match = 0;
        vector<TVector3> vec_TV3_ITS_hits, vec_TV3_ITS_hits_raw;
        TVector3 TV3_ITS_hit, ITS_raw_hit;
        for(Int_t i_hit = 0; i_hit < NITSHit; i_hit++)
        {
            ITSHit = TrackHitEvent->getITSHit(i_hit);
            Float_t  x_cls        = ITSHit ->get_cluster_x();
            Float_t  y_cls        = ITSHit ->get_cluster_y();
            Float_t  z_cls        = ITSHit ->get_cluster_z();
            int64_t  BC_cls       = ITSHit ->get_BC();
            UShort_t row          = ITSHit ->get_row();
            UShort_t col          = ITSHit ->get_col();
            UShort_t id           = ITSHit ->get_id();
            Float_t  time_bin     = 400.0 + ((Float_t)BC_cls)/8.0;  // 400.0 ???
            if(id > max_ITS_chip_id) max_ITS_chip_id = id;
            if(row > max_ITS_row) max_ITS_row = row;
            if(col > max_ITS_col) max_ITS_col = col;

            itterator = std::find(vec_ITS_id_row_column_noisy[id][row].begin(), vec_ITS_id_row_column_noisy[id][row].end(), col);
            if(itterator != vec_ITS_id_row_column_noisy[id][row].end() )
            {
                // pixel flagged as noisy
                continue;
            }
            else
            {
                N_good_ITS_hits++;
            }
            num_ITS_hits.push_back(N_good_ITS_hits);

            vec_cls_point[0] = x_cls;
            vec_cls_point[1] = y_cls;
            vec_cls_point[2] = time_bin;

            TVector3 TV3_ITS_hit;
            TV3_ITS_hit.SetXYZ(x_cls,y_cls,time_bin);
            ITS_raw_hit.SetXYZ(x_cls,y_cls,z_cls);
            vec_TV3_ITS_hits.push_back(TV3_ITS_hit);
            vec_TV3_ITS_hits_raw.push_back(ITS_raw_hit);
            for(Int_t i_track = 0; i_track < (Int_t)vec_TV3_dir.size(); i_track++)
            {
                TVector3 TV3_dca_ITS_hit = calculateDCA_vec_StraightToPoint(vec_TV3_base[i_track],vec_TV3_dir[i_track],TV3_ITS_hit);
                if(fabs(TV3_dca_ITS_hit.Z()) < 100.0 && TV3_dca_ITS_hit.Perp() < 3.0)
                {
                    if(vec_N_TPC_clusters[i_track] > 30)
                    {
                        printf("Match at TF: %lld, track: %d, ITS hit: {%4.3f, %4.3f, %4.3f}, id,row,col: {%d, %d, %d} \n",counter,i_track,x_cls,y_cls,time_bin,id,row,col);
                    }
                    if(counter == event_plot)
                    {
                        tg_ITS_hit_x_vs_y ->SetPoint(i_point_ITS_match,x_cls,y_cls);
                        //if(i_track == 5) vec_TV3_ITS_hits.push_back(TV3_ITS_hit);
                        i_point_ITS_match++;
                    }
                }
            }


            if(counter == event_plot)
            {
                h_cls_y_vs_x ->Fill(x_cls,y_cls);
                h_cls_y_vs_x_coarse ->Fill(x_cls,y_cls);
                vec_cls_points.push_back(vec_cls_point);
                tg_ITS_hit_x_vs_time ->SetPoint(N_good_ITS_hits-1,time_bin,x_cls);
                tg_ITS_hit_y_vs_time ->SetPoint(N_good_ITS_hits-1,time_bin,y_cls);
                //printf("ITS i_hit: %d, time_bin: %4.3f, row/col/chipid: {%d, %d, %d} \n",i_hit,time_bin,row,col,id);
            }
        }

        // cout << "Processing ITS hits" << endl; 
        // vector<float> radii_ITS_hits, sorted_radii_indices;
        TVector3 current_its_cluster;
        vector<TVector3> temp_its_clusters;
        vector<float> residuals_DCA_ITS, residual_DCA_XYZ_ITS, residuals_outliers_ITS;
        float max_distance=0, its_residual_dca=0; int id_hit1=0, id_hit2=0, counter_hit1=0, counter_hit2=0;
        bool print_event = false;
        Int_t found_outlier = 0;
        TVector3 ITS_outer_base, ITS_outer_dir;
        if((Int_t)vec_TV3_ITS_hits_raw.size() > 6){
            // for(int i_hit = 0; i_hit < (Int_t)vec_TV3_ITS_hits_raw.size(); i_hit++){
            //     radii_ITS_hits.push_back(TMath::Sqrt(TMath::Power(vec_TV3_ITS_hits_raw[i_hit][0],2) + TMath::Power(vec_TV3_ITS_hits_raw[i_hit][1],2)));
            // }
            // sorted_radii_indices = sort_indices(radii_ITS_hits);
            // cout << "Vector and sorting indices:" << endl;
            // for(auto elem : sorted_radii_indices){
            //     cout << elem << " ";
            // }
            // cout << endl;
            // for(auto elem : radii_ITS_hits){
            //     cout << elem << " ";
            // }
            // cout << endl;

            // Calculate mutual distance and find two most distant points
            for(auto hit1 : vec_TV3_ITS_hits_raw){
                for(auto hit2 : vec_TV3_ITS_hits_raw){
                    if(max_distance<(hit1-hit2).Mag()){
                        max_distance = (hit1-hit2).Mag();
                        id_hit1 = counter_hit1;
                        id_hit2 = counter_hit2;
                    }
                    counter_hit2++;
                }
                counter_hit2=0;
                counter_hit1++;
            }
            ITS_outer_dir.SetXYZ(vec_TV3_ITS_hits_raw[id_hit1][0]-vec_TV3_ITS_hits_raw[id_hit2][0],
                                 vec_TV3_ITS_hits_raw[id_hit1][1]-vec_TV3_ITS_hits_raw[id_hit2][1],
                                 vec_TV3_ITS_hits_raw[id_hit1][2]-vec_TV3_ITS_hits_raw[id_hit2][2]);
            ITS_outer_base.SetXYZ(vec_TV3_ITS_hits_raw[id_hit1][0],
                                  vec_TV3_ITS_hits_raw[id_hit1][1],
                                  vec_TV3_ITS_hits_raw[id_hit1][2]);

            for(Int_t hits_in_between = 0; hits_in_between<(Int_t)vec_TV3_ITS_hits_raw.size(); hits_in_between++){
                current_its_cluster = vec_TV3_ITS_hits_raw[hits_in_between];
                temp_its_clusters.push_back(current_its_cluster);
                // current_its_cluster.SetZ(0);
                its_residual_dca = calculateMinimumDistanceStraightToPoint(ITS_outer_base,ITS_outer_dir,current_its_cluster);
                residual_DCA_XYZ_ITS.push_back(TMath::Abs(calculateDCA_vec_StraightToPoint(ITS_outer_base,ITS_outer_dir,current_its_cluster)[0]));
                residual_DCA_XYZ_ITS.push_back(TMath::Abs(calculateDCA_vec_StraightToPoint(ITS_outer_base,ITS_outer_dir,current_its_cluster)[1]));
                residual_DCA_XYZ_ITS.push_back(TMath::Abs(calculateDCA_vec_StraightToPoint(ITS_outer_base,ITS_outer_dir,current_its_cluster)[2]));
                residuals_DCA_ITS.push_back(its_residual_dca);
                residuals_DCA_ITS_histo->Fill(its_residual_dca);
                if(its_residual_dca>0.04 && its_residual_dca<0.06){
                    found_outlier=1; // cut can be adjusted
                    print_event=true;
                }
            }
            if(print_event){
                Int_t count_loop = 0;
                for(Int_t hits_in_between = 0; hits_in_between<(Int_t)vec_TV3_ITS_hits_raw.size(); hits_in_between++){
                    cout << "DCA_total = " << residuals_DCA_ITS[count_loop] << "; DCA_X = " << residual_DCA_XYZ_ITS[3*count_loop + 0] << "; DCA_Y = " << residual_DCA_XYZ_ITS[3*count_loop +1] << "; DCA_Z = " << residual_DCA_XYZ_ITS[3*count_loop +2] << "; X = " << temp_its_clusters[count_loop][0] << "; Y = " << temp_its_clusters[count_loop][1] << "; Z = " << temp_its_clusters[count_loop][2] << endl;
                    count_loop++;
                }
                cout << "Total clusters ITS: " << count_loop/3 << endl;
            }
        }
        else{
            residuals_DCA_ITS.push_back(-1);
            // residuals_angles_ITS.push_back(10);
            residuals_DCA_ITS_histo->Fill(-1);
            // sresiduals_angles_ITS_histo->Fill(10);
        }
        residuals_outliers_ITS.push_back(found_outlier);
        residuals_outliers_ITS_histo->Fill(found_outlier);
        // if(found_outlier){
        //     cout << its_residual_dca << endl;
        //     for(auto elem : vec_TV3_ITS_hits_raw[id_hit1]){
        //         cout << elem << " ";
        //     }
        //     cout << endl;
        //     for(auto elem : vec_TV3_ITS_hits_raw[id_hit2]){
        //         cout << elem << " ";
        //     }
        //     cout << endl;
        // }
        // radii_ITS_hits.clear(); sorted_radii_indices.clear();
        //printf("N_good_ITS_hits: %d \n",N_good_ITS_hits);

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
        for(Int_t i_hit = 0; i_hit < (Int_t)vec_TV3_ITS_hits.size(); i_hit++)
        {
            vec_ITS_hits_used[i_hit] = 0;
        }
        for(Int_t i_hit = 0; i_hit < (Int_t)vec_TV3_ITS_hits.size()-1; i_hit++)
        {
            if(vec_ITS_hits_used[i_hit]) continue;
            vector<TVector3> vec_TV3_ITS_hits_same_time;
            Double_t time_hit = vec_TV3_ITS_hits[i_hit].Z();
            vec_TV3_ITS_hits_same_time.push_back(vec_TV3_ITS_hits[i_hit]);
            vec_ITS_hits_used[i_hit] = 1;
            for(Int_t i_hitB = (i_hit+1); i_hitB < (Int_t)vec_TV3_ITS_hits.size(); i_hitB++)
            {
                if(vec_ITS_hits_used[i_hitB]) continue;
                Double_t time_hitB = vec_TV3_ITS_hits[i_hitB].Z();
                if(time_hit != time_hitB) continue;
                vec_TV3_ITS_hits_same_time.push_back(vec_TV3_ITS_hits[i_hitB]);
                vec_ITS_hits_used[i_hitB] = 1;
            }

            Int_t N_hits_same_time = (Int_t)vec_TV3_ITS_hits_same_time.size();
            Int_t N_hits_same_time_close = 0;
            if(N_hits_same_time > 2)
            {
                TVector3 TV3_base = vec_TV3_ITS_hits_same_time[0];
                TVector3 TV3_dir  = vec_TV3_ITS_hits_same_time[0] - vec_TV3_ITS_hits_same_time[1];
                for(Int_t i_hitC = 2; i_hitC < N_hits_same_time; i_hitC++)
                {
                    TVector3 TV3_dca_ITS_hit = calculateDCA_vec_StraightToPoint(TV3_base,TV3_dir,vec_TV3_ITS_hits_same_time[i_hitC]);
                    if(TV3_dca_ITS_hit.Perp() < 0.3)
                    {
                        N_hits_same_time_close++;
                        //printf("TF: %lld \n",counter);
                        //TV3_dca_ITS_hit.Print();

                    }
                }
            }
            if(N_hits_same_time_close > 3) printf("TF: %lld, N_hits_same_time_close: %d \n",counter,N_hits_same_time_close);
        }
        //--------------------



        //--------------------
        // Check for TPC to TPC track matches
        Int_t flag_good_TPC_match = 0;
        for(Int_t i_trackA = 0; i_trackA < (Int_t)vec_TV3_dir.size(); i_trackA++)
        {
            for(Int_t i_trackB = (i_trackA+1); i_trackB < (Int_t)vec_TV3_dir.size(); i_trackB++)
            {
                if(i_trackA == i_trackB) continue;
                TVector3 TV3_dca_TPC_hit = calculateDCA_vec_StraightToPoint(vec_TV3_base[i_trackA],vec_TV3_dir[i_trackA],vec_TV3_base[i_trackB]);
                Float_t dca_3D_lines = calculateDCA_vec_StraightToStraight(vec_TV3_base[i_trackA],vec_TV3_dir[i_trackA],vec_TV3_base[i_trackB],vec_TV3_dir[i_trackB]);
                printf("tracks: {%d, %d}, perp dist: %4.3f, 3D dca: %4.3f \n",i_trackA,i_trackB,TV3_dca_TPC_hit.Perp(),dca_3D_lines);
                //if(fabs(TV3_dca_TPC_hit.Z()) < 6.0 && TV3_dca_TPC_hit.Perp() < 3.0)
                // if(dca_3D_lines < 5.0 && tpc_ncls[i_trackA]>40 && tpc_ncls[i_trackB]>40 && (chi2_pca[i_trackA] < 5 || chi2_pca[i_trackB] < 5) && (num_ITS_hits[i_trackA] >= 5 || num_ITS_hits[i_trackB] >= 5))
                if((residuals_outliers_ITS[i_trackA] == 1 || residuals_outliers_ITS[i_trackB] == 1))
                {
                    flag_good_TPC_match = 1;
                    printf("Good TPC-to-TPC match found! tracks: {%d, %d} \n",i_trackA,i_trackB);
                    Double_t angle = TMath::ACos((vec_TV3_dir[i_trackA].Dot(vec_TV3_dir[i_trackB]))/(vec_TV3_dir[i_trackA].Mag()*vec_TV3_dir[i_trackB].Mag()));
                    cout<< "Angle: " << angle << endl;
                    // if(angle < 0.1 && angle > 0.){
                    //     outlier_events->Fill(counter);
                    // }
                    outlier_events->Fill(counter);
                    angle_dist->Fill(angle);
                }
            }
        }
        //--------------------



        //printf("Time frame: %lld, N_good_ITS_hits: %d \n",counter,N_good_ITS_hits);
    } // end if time frame loop
    //------------------------------------------------------------------------------




    printf("time: {%4.3f, %4.3f}, max_ITS_chip_id: %d, max_ITS_row: %d, max_ITS_col: %d, N_pixels_hit: %lld, N_noisy_ITS_pixels: %d \n",time_min,time_max,max_ITS_chip_id,max_ITS_row,max_ITS_col,N_pixels_hit,N_noisy_ITS_pixels);
    TH2D*   h_cls_x_vs_time  = new TH2D("h_cls_x_vs_time","h_cls_x_vs_time",1000,time_min,time_max,1000,-250,250);
    TH2D*   h_cls_y_vs_time  = new TH2D("h_cls_y_vs_time","h_cls_y_vs_time",1000,time_min,time_max,1000,-250,250);

    for(Int_t i_cls = 0; i_cls < (Int_t)vec_cls_points.size(); i_cls++)
    {
        h_cls_x_vs_time ->Fill(vec_cls_points[i_cls][2],vec_cls_points[i_cls][0]);
        h_cls_y_vs_time ->Fill(vec_cls_points[i_cls][2],vec_cls_points[i_cls][1]);
    }
    //------------------------------------------------------------




    //------------------------------------------------------------
    Draw_1D_histo_and_canvas(angle_dist, "angular_distribution", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(outlier_events, "outlier_events", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(chi2_tracks, "chi2_tracks", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(chi2_tracks_xy, "chi2_tracks_xy", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(chi2_tracks_time, "chi2_tracks_time", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(residuals_DCA_ITS_histo, "residuals_DCA_ITS_histo", 720, 720, 720, 720, "");
    Draw_1D_histo_and_canvas(residuals_outliers_ITS_histo, "residuals_outliers_ITS_histo", 720, 720, 720, 720, "");

    h_cls_y_vs_x ->GetXaxis()->SetTitle("x (cm)");
    h_cls_y_vs_x ->GetYaxis()->SetTitle("y (cm)");
    h_cls_y_vs_x ->GetZaxis()->SetTitle("entries");
    TCanvas* can_h_cls_y_vs_x = Draw_2D_histo_and_canvas((TH2D*)h_cls_y_vs_x,"can_h_cls_y_vs_x",1100,800,0.0,0.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_h_cls_y_vs_x->SetLogz(0);
    Draw_Circle_Detector_2D(85.225,250,2,18,kBlack,1,2,0.0,0.0); // (Float_t radius_in = 1, Float_t radius_out = 2,const Int_t n_radii = 1, const Int_t n_delta_phi = 2, Float_t color = 2, Int_t line_style = 1, Int_t line_width = 1, Float_t x_offset = 0.0, Float_t y_offset = 0.0)

    vector<TPolyLine*> vec_TPL_track_y_vs_x;
    Int_t arr_color[10] = {kRed,kGreen+2,kMagenta,kCyan,kOrange+2,kAzure-2,kYellow+2,kPink+2,kViolet+5,kTeal+7};
    for(Int_t i_track = 0; i_track < (Int_t)vec_track_par_helix.size(); i_track++)
    {
        //printf("i_track: %d \n",i_track);
        vec_TPL_track_y_vs_x.push_back(new TPolyLine());
        for(Int_t i_par = 0; i_par < 6; i_par++)
        {
            fHelix[i_par] = vec_track_par_helix[i_track][i_par];
            //printf("i_par: %d, par: %4.3f \n",i_par,fHelix[i_par]);
        }
        for(Float_t track_path = -85.0; track_path < 350.0; track_path += 1.0)
        {
            evaluate_helix(track_path,track_pos);
            Double_t radius_track = TMath::Sqrt(track_pos[0]*track_pos[0] + track_pos[1]*track_pos[1]);
            if(radius_track > 250.0) continue;
            vec_TPL_track_y_vs_x[i_track] ->SetNextPoint(track_pos[0],track_pos[1]);
        }
        if(i_track < 10) vec_TPL_track_y_vs_x[i_track] ->SetLineColorAlpha(arr_color[i_track],0.2);
        else  vec_TPL_track_y_vs_x[i_track] ->SetLineColorAlpha(kOrange+2,0.2);
        vec_TPL_track_y_vs_x[i_track] ->SetLineWidth(3);
        vec_TPL_track_y_vs_x[i_track] ->DrawClone();

        vec_tg_cls_y_vs_x[i_track] ->SetMarkerStyle(24);
        vec_tg_cls_y_vs_x[i_track] ->SetMarkerSize(0.5);
        if(i_track < 10) vec_tg_cls_y_vs_x[i_track] ->SetMarkerColor(arr_color[i_track]);
        else vec_tg_cls_y_vs_x[i_track] ->SetMarkerColor(kGray+2);
        vec_tg_cls_y_vs_x[i_track] ->Draw("same p");
    }

    tg_ITS_hit_x_vs_y ->SetMarkerStyle(20);
    tg_ITS_hit_x_vs_y ->SetMarkerSize(0.5);
    tg_ITS_hit_x_vs_y ->SetMarkerColor(kRed);
    tg_ITS_hit_x_vs_y ->Draw("same p");
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
    h_cls_x_vs_time ->GetXaxis()->SetTitle("time (tb)");
    h_cls_x_vs_time ->GetYaxis()->SetTitle("x (cm)");
    h_cls_x_vs_time ->GetZaxis()->SetTitle("entries");
    TCanvas* can_h_cls_x_vs_time = Draw_2D_histo_and_canvas((TH2D*)h_cls_x_vs_time,"can_h_cls_x_vs_time",1100,800,0.0,0.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_h_cls_x_vs_time->SetLogz(0);

    for(Int_t i_track = 0; i_track < (Int_t)vec_track_par_helix.size(); i_track++)
    {
        vec_tg_cls_x_vs_time[i_track] ->SetMarkerStyle(24);
        vec_tg_cls_x_vs_time[i_track] ->SetMarkerSize(0.5);
        if(i_track < 10) vec_tg_cls_x_vs_time[i_track] ->SetMarkerColor(arr_color[i_track]);
        else vec_tg_cls_x_vs_time[i_track] ->SetMarkerColor(kGray+2);
        vec_tg_cls_x_vs_time[i_track] ->Draw("same p");
    }


   tg_ITS_hit_x_vs_time  ->SetMarkerStyle(20);
   tg_ITS_hit_x_vs_time  ->SetMarkerSize(0.5);
   tg_ITS_hit_x_vs_time  ->SetMarkerColor(kRed);
   tg_ITS_hit_x_vs_time  ->Draw("same p");
    //------------------------------------------------------------


    //------------------------------------------------------------
    h_cls_y_vs_time ->GetXaxis()->SetTitle("time (tb)");
    h_cls_y_vs_time ->GetYaxis()->SetTitle("y (cm)");
    h_cls_y_vs_time ->GetZaxis()->SetTitle("entries");
    TCanvas* can_h_cls_y_vs_time = Draw_2D_histo_and_canvas((TH2D*)h_cls_y_vs_time,"can_h_cls_y_vs_time",1100,800,0.0,0.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_h_cls_y_vs_time->SetLogz(0);

    for(Int_t i_track = 0; i_track < (Int_t)vec_track_par_helix.size(); i_track++)
    {
        vec_tg_cls_y_vs_time[i_track] ->SetMarkerStyle(24);
        vec_tg_cls_y_vs_time[i_track] ->SetMarkerSize(0.5);
        if(i_track < 10) vec_tg_cls_y_vs_time[i_track] ->SetMarkerColor(arr_color[i_track]);
        else vec_tg_cls_y_vs_time[i_track] ->SetMarkerColor(kGray+2);
        vec_tg_cls_y_vs_time[i_track] ->Draw("same p");
    }


    tg_ITS_hit_y_vs_time  ->SetMarkerStyle(20);
    tg_ITS_hit_y_vs_time  ->SetMarkerSize(0.5);
    tg_ITS_hit_y_vs_time  ->SetMarkerColor(kRed);
    tg_ITS_hit_y_vs_time  ->Draw("same p");
   //------------------------------------------------------------



    //------------------------------------------------------------
    if(flag_ITS_noisy)
    {
        file_ITS_noisy_pixels ->cd();
        output_tree.Write();
        file_ITS_noisy_pixels ->Close();
    }
    //------------------------------------------------------------

    delete PCA;

}