
#include "functions_2.h"

void Ana_cosmics()
{
    printf("Ana_cosmics started \n");




    //--------------------------------------
    //TFile* inputfile = TFile::Open("./QA_tree.root");
    TFile* inputfile = TFile::Open("./QA_tree.root");

    //TTree *input_tree = (TTree*)inputfile->Get("tpc_raw_dca_wCrit");
    TTree *input_tree = (TTree*)inputfile->Get("tpc_raw_dca_wCrit");
    Float_t x,y,time,sector,row,dca_x,dca_y,dca_time;

    input_tree->SetBranchAddress("x",&x);
    input_tree->SetBranchAddress("y",&y);
    input_tree->SetBranchAddress("time",&time);
    input_tree->SetBranchAddress("dca_x",&dca_x);
    input_tree->SetBranchAddress("dca_y",&dca_y);
    input_tree->SetBranchAddress("sector",&sector);
    input_tree->SetBranchAddress("row",&row);
    input_tree->SetBranchAddress("dca_time",&dca_time);

    Int_t nentries = (Int_t)input_tree->GetEntries();
    printf("nentries: %d \n",nentries);
    //--------------------------------------




    //--------------------------------------
    const Int_t N_sectors = 36;
    vector<TH2D*> vec_h2D_DX_vs_row;
    vector<TProfile*> vec_TP_DX_vs_row;
    vector<TCanvas*> vec_can_vec_h2D_DX_vs_row;
    vec_h2D_DX_vs_row.resize(N_sectors);
    vec_TP_DX_vs_row.resize(N_sectors);
    vec_can_vec_h2D_DX_vs_row.resize(N_sectors);
    for(Int_t i_sec = 0; i_sec < N_sectors; i_sec++)
    {
        HistName = "vec_h2D_DX_vs_row_";
        HistName += i_sec;
        vec_h2D_DX_vs_row[i_sec] = new TH2D(HistName.Data(),HistName.Data(),161,0,160,400,-2,2);

        HistName = "vec_hTP_DX_vs_row_";
        HistName += i_sec;
        vec_TP_DX_vs_row[i_sec] = new TProfile(HistName.Data(),HistName.Data(),161,0,160);
    }
    //--------------------------------------




    //--------------------------------------
    TVector3 TV3_sec_center, TV3_sec_center_perp, TV3_dca, TV3_z_axis;
    TV3_z_axis.SetXYZ(0.0,0.0,1.0);
    for(Int_t ientry = 0; ientry < nentries; ientry++)
    {
        if(ientry % 1000000 == 0) printf("ientry: %d, out of %d \n",ientry,nentries);
        input_tree->GetEntry(ientry);
        if(dca_x != dca_x || dca_y != dca_y) continue;
        Float_t phi_sector = (10.0 + (((Int_t)sector)%18)*20.0)*TMath::DegToRad();
        TV3_sec_center.SetPtEtaPhi(1.0,1.0,phi_sector);
        TV3_sec_center_perp = TV3_sec_center.Cross(TV3_z_axis);
        TV3_dca.SetXYZ(dca_x,dca_y,0.0);
        Float_t DX = TV3_dca.Dot(TV3_sec_center);
        Float_t DY = TV3_dca.Dot(TV3_sec_center_perp);
        if(fabs(DX) > 5.0) continue;
        //printf("ientry: %d, row: %4.1f, DX: %4.3f dca: {%4.3f, %4.3f} \n",ientry,row,DX,dca_x,dca_y);
        vec_h2D_DX_vs_row[sector] ->Fill(row,DX);
        vec_TP_DX_vs_row[sector]  ->Fill(row,DX);
    }
    //--------------------------------------




    //--------------------------------------
    Int_t sector_plot = 9;

    for(Int_t i_sec = 0; i_sec < N_sectors; i_sec++)
    {
        sector_plot = i_sec;
        vec_h2D_DX_vs_row[sector_plot] ->GetXaxis()->SetTitle("pad row");
        vec_h2D_DX_vs_row[sector_plot] ->GetYaxis()->SetTitle("#DeltaX (cm)");
        vec_h2D_DX_vs_row[sector_plot] ->GetZaxis()->SetTitle("entries");

        HistName = "vec_can_vec_h2D_DX_vs_row_";
        HistName += i_sec;
        vec_can_vec_h2D_DX_vs_row[i_sec] = Draw_2D_histo_and_canvas((TH2D*)vec_h2D_DX_vs_row[sector_plot],HistName.Data(),1100,800,0.0,0.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
        vec_TP_DX_vs_row[sector_plot] ->SetMarkerStyle(20);
        vec_TP_DX_vs_row[sector_plot] ->SetMarkerSize(0.65);
        vec_TP_DX_vs_row[sector_plot] ->SetMarkerColor(kRed);
        vec_TP_DX_vs_row[sector_plot] ->DrawCopy("same P");
        vec_can_vec_h2D_DX_vs_row[i_sec]->SetLogz(0);
        HistName = "sector: ";
        HistName += i_sec;
        plotTopLegend((char*)HistName.Data(),0.25,0.85,0.05,1,0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    }
    //--------------------------------------




}