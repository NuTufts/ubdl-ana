#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLine.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"

#include "larflow/Reco/NuSelectionVariables.h"
#include "larflow/Reco/NuVertexCandidate.h"
#include "larflow/Reco/cluster_functions.h"

// define helper functions here
//start by loading all kpsreco files
TFile* in_f;
// load into the branches
TTree* in;
int _run = -1;
int _subrun= -1;
int _event =-1;
int _ispi0=-1;
int _best_passing_vtx_index=-1;
std::vector<float>* _leading_shower_ADCU = NULL;
std::vector<float>* _leading_shower_ADCV= NULL;
std::vector<float>* _leading_shower_ADCY= NULL;
std::vector<float>* _subleading_shower_ADCU= NULL;
std::vector<float>* _subleading_shower_ADCV= NULL;
std::vector<float>* _subleading_shower_ADCY= NULL;
std::vector<float>* _leading_shower_Xdir= NULL;
std::vector<float>* _subleading_shower_Xdir= NULL;
std::vector<float>* _leading_shower_Ydir= NULL;
std::vector<float>* _subleading_shower_Ydir= NULL;
std::vector<float>* _leading_shower_Zdir= NULL;
std::vector<float>* _subleading_shower_Zdir= NULL;
std::vector<float>* _overlap_frac= NULL;

int main( int nargs, char** argv ) {
  //script to load in gen2 reco files, get shower calibration variables, and produce plots
  std::cout << "Analyze Energy Calibration Plots" << std::endl;
  std::cout << "Args: anafile outdir"<<std::endl;
  in_f = new TFile(argv[1],"read");
  // load into the branches
  in = (TTree*)in_f->Get("pi0reco_anatree");

  in->SetBranchAddress("Run",&_run);
  in->SetBranchAddress("Subrun",&_subrun);
  in->SetBranchAddress("Event",&_event);
  in->SetBranchAddress("ispi0",&_ispi0);
  in->SetBranchAddress("best_passing_vtx_index",&_best_passing_vtx_index);
  in->SetBranchAddress("overlap_frac",&_overlap_frac);
  in->SetBranchAddress("leading_shower_ADCU",&_leading_shower_ADCU);
  in->SetBranchAddress("leading_shower_ADCV",&_leading_shower_ADCV);
  in->SetBranchAddress("leading_shower_ADCY",&_leading_shower_ADCY);
  in->SetBranchAddress("subleading_shower_ADCU",&_subleading_shower_ADCU);
  in->SetBranchAddress("subleading_shower_ADCV",&_subleading_shower_ADCV);
  in->SetBranchAddress("subleading_shower_ADCY",&_subleading_shower_ADCY);
  in->SetBranchAddress("leading_shower_Xdir",&_leading_shower_Xdir);
  in->SetBranchAddress("subleading_shower_Xdir",&_subleading_shower_Xdir);
  in->SetBranchAddress("leading_shower_Ydir",&_leading_shower_Ydir);
  in->SetBranchAddress("subleading_shower_Ydir",&_subleading_shower_Ydir);
  in->SetBranchAddress("leading_shower_Zdir",&_leading_shower_Zdir);
  in->SetBranchAddress("subleading_shower_Zdir",&_subleading_shower_Zdir);

  // define outputdir
  std::string outdir = argv[2];

  // initialize histograms
  TH1F* overlap_frac_h = new TH1F("overlap_frac","overlap_frac",25,0.0,0.5);
  TH1F* leading_shower_ADCY_h = new TH1F("leading_shower_ADCY","leading_shower_ADCY",50,0.0,200000);
  TH1F* subleading_shower_ADCY_h = new TH1F("subleading_shower_ADCY","subleading_shower_ADCY",50,0.0,200000);
  TH1F* leading_shower_EY_h = new TH1F("leading_shower_EY","leading_shower_EY",20,0.0,1200);
  TH1F* subleading_shower_EY_h = new TH1F("subleading_shower_EY","subleading_shower_EY",20,0.0,1200);


  // loop through entries
  int nentries = in->GetEntries();
  for (int ientry=0; ientry<nentries; ientry++) {
  // for (int ientry=0; ientry<1; ientry++) {
    std::cout<<"Entry: "<<ientry<<std::endl;
    // get the entry for each io manager
    in->GetEntry(ientry);
    // loop over vertex
    if (_overlap_frac->size() >0){
      for (int idx = 0; idx <_overlap_frac->size();idx++){

        if (_overlap_frac->at(idx) >= 0.0 && _best_passing_vtx_index == idx && _ispi0==1){
          overlap_frac_h->Fill(_overlap_frac->at(idx));
          leading_shower_ADCY_h->Fill(_leading_shower_ADCY->at(idx));
          subleading_shower_ADCY_h->Fill(_subleading_shower_ADCY->at(idx));
          leading_shower_EY_h->Fill(_leading_shower_ADCY->at(idx)*0.012607+22.1);
          std::cout<<_leading_shower_ADCY->at(idx)*0.012607+22.1<<std::endl;
          subleading_shower_EY_h->Fill(_subleading_shower_ADCY->at(idx)*0.012607+22.1);
        }
      }//end of vertex loop
    }
  }//end of entry loop
  // save to a png and format
  // gStyle->SetOptStat(0);

  TCanvas can("can", "histograms ", 1500, 1000);
  can.cd();
  overlap_frac_h->SetTitle("Overlap Fraction between showers");
  overlap_frac_h->SetXTitle("Overlap Fraction");
  overlap_frac_h->SetYTitle("Number of Entries");
  overlap_frac_h->Draw("");
  can.SaveAs(Form("%s/pi0_overlapfrac.png",outdir.c_str()));

  TCanvas can1("can", "histograms ", 1500, 1000);
  can1.cd();
  leading_shower_ADCY_h->SetTitle("Shower ADC");
  leading_shower_ADCY_h->SetXTitle("ADC Sum");
  leading_shower_ADCY_h->SetYTitle("Number of Entries");
  leading_shower_ADCY_h->SetLineColor(kRed);
  leading_shower_ADCY_h->Draw("");
  subleading_shower_ADCY_h->SetLineColor(kBlue);
  subleading_shower_ADCY_h->Draw("SAME");
  can1.SaveAs(Form("%s/pi0_showerADCY.png",outdir.c_str()));

  TCanvas can2("can", "histograms ", 1500, 1000);
  can2.cd();
  subleading_shower_EY_h->SetTitle("Shower Energy");
  subleading_shower_EY_h->SetXTitle("Reco Energy (MeV)");
  subleading_shower_EY_h->SetYTitle("Number of Entries");
  leading_shower_EY_h->SetLineColor(kRed);

  subleading_shower_EY_h->SetLineColor(kBlue);
  subleading_shower_EY_h->Draw("");
  leading_shower_EY_h->Draw("SAME");
  can2.SaveAs(Form("%s/pi0_showerEY.png",outdir.c_str()));

  // delete pointers
  delete overlap_frac_h;
  return 0;
}
