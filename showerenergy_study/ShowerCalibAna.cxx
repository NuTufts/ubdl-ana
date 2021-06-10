#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
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

int main( int nargs, char** argv ) {
  //script to load in gen2 reco files, get shower calibration variables, and produce plots
  std::cout << "Analyze Energy Calibration Plots" << std::endl;
  std::cout << "Args: plotfile outdir"<<std::endl;

  //start by loading all kpsreco files
  TFile* plots_f = new TFile(argv[1],"read");
  // load into the histograms we need
  TH2D* ADCHist_plane02 = (TH2D*)plots_f->Get("hnshowerADC_vs_evislep_1e1p_2");

  // define outputdir
  std::string outdir = argv[2];

  // first step, for each vertical bin get a vector of z values
  // loop through x bins
  int xbins = 50;
  int ybins = 50;
  // initialize root points
  TH1D* verticalbins[xbins] = {0};
  TF1* f1[xbins] = {0};
  Double_t xarr[xbins], yarr[xbins];
  Double_t xarr_error[xbins], yarr_error[xbins];
  int nbinsnew=0;
  for(int x =0;x<xbins;x++){
    std::stringstream ss;
    ss << "vericalbin_" << x;
    verticalbins[x] = new TH1D(ss.str().c_str(),ss.str().c_str(),ybins,0,2500);
    // loop through ybins - first find peak
    // int maxidx = 0;
    // int maxval = 0;
    // for (int y=1;y<=ybins;y++){
    //   int zval = ADCHist_plane02->GetBinContent(x,y);
    //   if (zval > maxval){
    //     maxidx = y;
    //     maxval = zval;
    //   }
    // }
    int totalinbin=0;
    for (int y=0;y<ybins;y++){
      int zval = ADCHist_plane02->GetBinContent(x,y);
      verticalbins[x]->SetBinContent(y,zval);
      totalinbin+=zval;
    }
    std::cout<<"bin number: "<< x << " total in bin: "<< totalinbin <<" stat error: "<<sqrt(totalinbin)/(float)totalinbin<<std::endl;
    if (totalinbin>400){
      f1[x] = new TF1("f1","gaus",0,2500);
      // f1[x]->SetParLimits(1,200,2500);
      verticalbins[x]->Fit("f1");
      std::cout<<"bin: "<<x<<" mean "<<f1[x]->GetParameters()[1]<<std::endl;
      // fill arrays with fit points
      xarr[x] = x*float(200000)/float(xbins)+float(200000)/float(2*xbins);
      xarr_error[x] = float(200000)/float(2.0*xbins);
      yarr[x] = f1[x]->GetParameters()[1];
      yarr_error[x] = (f1[x]->GetParameters()[2])/2.0;
      nbinsnew+=1;
    }
    else{
      xarr[x] = -1;
      xarr_error[x] = -1;
      yarr[x] = -1;
      yarr_error[x] = -1;
    }


  }

  // go through to remove indexes without a fit
  Double_t xarr_reduced[nbinsnew], yarr_reduced[nbinsnew];
  Double_t xarr_error_reduced[nbinsnew], yarr_error_reduced[nbinsnew];
  int i =1;
  for (int x =1;x<=xbins;x++){
    if (xarr_error[x] >= 0 && yarr[x] > 0){
      xarr_reduced[i] = xarr[x];
      yarr_reduced[i] = yarr[x];
      xarr_error_reduced[i] = xarr_error[x];
      yarr_error_reduced[i] = yarr_error[x];
      i+=1;
    }
  }

  // make scatter plot of fit points
  // TGraphErrors* fitpts = new TGraphErrors(xbins,xarr,yarr,xarr_error,yarr_error);
  TGraphErrors* fitpts = new TGraphErrors(nbinsnew,xarr_reduced,yarr_reduced,xarr_error_reduced,yarr_error_reduced);
  // TGraph* fitline = new TGraph(xbins,xarr,yarr);
  TGraph* fitline = new TGraph(nbinsnew,xarr_reduced,yarr_reduced);
  TF1* p1 = new TF1("p1","pol1",0,200000);
  fitline->Fit("p1");

  // now I need a gaussian fit to each vertical bin

  // make gen 1 line
  TLine* gen1line = new TLine(0.0, 0.0, 2500.0/0.01256,2500.0);
  // make gen 1 line
  TLine* gen2line = new TLine(0.0, p1->GetParameters()[0], 2500.0/p1->GetParameters()[1],2500.0);
  // get all points to plot
  // TGraphErrors* fitpts_all = new TGraphErrors(xbins,xarr,yarr,xarr_error,yarr_error);
  // save to a png and format
  gStyle->SetOptStat(0);

  TCanvas can("can", "histograms ", 1500, 1000);
  can.cd();
  ADCHist_plane02->SetTitle("Shower Calibration Plot");
  ADCHist_plane02->SetXTitle("ADC Sum of Reco Shower");
  ADCHist_plane02->SetYTitle("True Electron Energy (MeV)");
  ADCHist_plane02->SetOption("COLZ");
  ADCHist_plane02->Draw("");
  gen1line->SetLineColor(kWhite);
  gen1line->SetLineStyle(2);
  gen1line->SetLineWidth(4);
  gen2line->SetLineColor(kRed);
  gen2line->SetLineStyle(1);
  gen2line->SetLineWidth(6);
  // fitpts_all->SetMarkerStyle(20);
  // fitpts_all->SetMarkerColor(kMagenta);
  // fitpts_all->Draw("SAME P");
  fitpts->SetMarkerStyle(20);
  fitpts->SetMarkerColor(kBlack);
  fitpts->Draw("SAME P");
  gen1line->Draw("SAME");
  gen2line->Draw("SAME");
  // can.SetLogz(1);

  can.SaveAs(Form("%s/showercalib_plane2.png",outdir.c_str()));

  // for(int x =0;x<xbins;x++){
  //   TCanvas can("can", "histograms ", 1500, 1000);
  //   can.cd();
  //   float minE = (x)*(200000.0/float(xbins));
  //   float maxE = (x+1)*(200000.0/float(xbins));
  //   std::stringstream minEs;
  //   minEs << minE;
  //   std::stringstream maxEs;
  //   maxEs << maxE;
  //
  //   verticalbins[x]->SetTitle(Form("ADCBin_%s_to_%s",minEs.str().c_str(),maxEs.str().c_str()));
  //   verticalbins[x]->SetXTitle("True Shower Energy(MeV)");
  //   verticalbins[x]->SetYTitle("Number of Events");
  //   verticalbins[x]->SetMarkerStyle(kDot);
  //   verticalbins[x]->SetMarkerSize(5);
  //   verticalbins[x]->Draw("SCATTER");
  //   f1[x]->SetLineColor(kRed);
  //   f1[x]->Draw("SAME L");
  //   can.SaveAs(Form("%s/ADCBin_%s_to_%s.png",outdir.c_str(),minEs.str().c_str(),maxEs.str().c_str()));
  // }//end of looping through bins


  // delete pointers
  delete ADCHist_plane02;
  return 0;
}
