
#include "functions.h"

void AnaResiduals()
{
    TFile* inputfile = TFile::Open("./Data/QA_tree.root");

    TTree *input_tree = (TTree*)inputfile->Get("residuals_ITS_full");
    Float_t x,y,z,dx,dy,dz,chi2,npoints;

    input_tree->SetBranchAddress("x",&x);
    input_tree->SetBranchAddress("y",&y);
    input_tree->SetBranchAddress("z",&z);
    input_tree->SetBranchAddress("dca_x",&dx);
    input_tree->SetBranchAddress("dca_y",&dy);
    input_tree->SetBranchAddress("dca_z",&dz);
    input_tree->SetBranchAddress("chi2",&chi2);
    input_tree->SetBranchAddress("n_points",&npoints);

    Int_t nentries = (Int_t)input_tree->GetEntries();
    printf("nentries: %d \n",nentries);


    //--------------------------------------
    TH2D* h2D_rad_res_vs_radius = new TH2D("h2D_rad_res_vs_radius","h2D_rad_res_vs_radius",100,0,45,100,-1.0,1.0);
    TProfile2D* TP2D_rad_res_vs_xy = new TProfile2D("TP2D_rad_res_vs_xy","TP2D_rad_res_vs_xy",100,-45,45,100,-45,45);
    //--------------------------------------



    //--------------------------------------
    for(Int_t ientry = 0; ientry < nentries; ientry++)
    {
        input_tree->GetEntry(ientry);

        if(fabs(dx) > 1.0 || fabs(dy) > 1.0 || fabs(dz) > 1.0) continue;
        if(chi2 > 5) continue;

        Float_t radius = TMath::Sqrt(x*x + y*y);

        TVector3 TV3_point;
        TV3_point.SetXYZ(x,y,0);
        Double_t mag_TV3_point = TV3_point.Mag();
        TV3_point *= 1.0/TV3_point.Mag();
        TVector3 TV3_res;
        TV3_res.SetXYZ(dx,dy,0);
        Double_t rad_res = TV3_res.Dot(TV3_point);

        //if(fabs(rad_res) > 5.0)
        //{
            //printf(" \n");
            //printf("rad_res: %4.3f, dx: %4.3f, dy: %4.3f, mag_TV3_point: %4.3f \n",rad_res,dx,dy,mag_TV3_point);
            //TV3_res.Print();
            //TV3_point.Print();
        //}

        h2D_rad_res_vs_radius ->Fill(radius,rad_res);
        TP2D_rad_res_vs_xy    ->Fill(x,y,rad_res);
    }
    //--------------------------------------


    //--------------------------------------
    h2D_rad_res_vs_radius ->GetXaxis()->SetTitle("radius (cm)");
    h2D_rad_res_vs_radius ->GetYaxis()->SetTitle("dr (cm)");
    h2D_rad_res_vs_radius ->GetZaxis()->SetTitle("entries");
    TProfile* TP_rad_res_vs_radius = (TProfile*)h2D_rad_res_vs_radius->ProfileX("TP_rad_res_vs_radius",1,-1);
    TCanvas* can_h2D_rad_res_vs_radius = Draw_2D_histo_and_canvas(h2D_rad_res_vs_radius,"can_h2D_rad_res_vs_radius",800,600,0.0,0.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_h2D_rad_res_vs_radius->SetLogz(0);
    TP_rad_res_vs_radius ->SetLineColor(kRed);
    TP_rad_res_vs_radius ->SetLineWidth(3);
    TP_rad_res_vs_radius ->DrawCopy("same hist");
    //--------------------------------------


    //--------------------------------------
    TP2D_rad_res_vs_xy ->GetXaxis()->SetTitle("x (cm)");
    TP2D_rad_res_vs_xy ->GetYaxis()->SetTitle("y (cm)");
    TP2D_rad_res_vs_xy ->GetZaxis()->SetTitle("dr (cm)");
    TCanvas* can_TP2D_rad_res_vs_xy = Draw_2D_histo_and_canvas((TH2D*)TP2D_rad_res_vs_xy,"can_TP2D_rad_res_vs_xy",800,600,0.0,0.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_TP2D_rad_res_vs_xy->SetLogz(0);
    //--------------------------------------




}