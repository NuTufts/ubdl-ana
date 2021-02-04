#include <iostream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

int main( int nargs, char** argv ) {

  std::cout << "Plot VA-reco Ana" << std::endl;

  TChain* in = new TChain("vtxactivityana");
  in->Add("out_vareco_ana.root");

  // true MC variables
  int is1l1p0pi;
  int is1l0p0pi;
  float evis_had;
  float evis_lep;
  float Enu_true;
  in->SetBranchAddress( "Enu_true", &Enu_true );
  in->SetBranchAddress( "evis_had", &evis_had );
  in->SetBranchAddress( "is1l1p0pi", &is1l1p0pi );
  in->SetBranchAddress( "is1l0p0pi", &is1l0p0pi );  
  
  // VA reco variables
  std::vector< std::vector<float> >* pca_dir_vv = 0;
  std::vector<int>* nbackwards_shower_pts = 0;
  std::vector<int>* nbackwards_track_pts = 0;
  std::vector<int>* nforwards_shower_pts = 0;
  std::vector<int>* nforwards_track_pts = 0;
  std::vector<int>* npix_on_cosmic_v = 0;
  std::vector<int>* attcluster_nall_v = 0;
  std::vector<int>* attcluster_nshower_v = 0;
  std::vector<int>* attcluster_ntrack_v = 0;
  std::vector<float>* dist_closest_forwardshower = 0;
  std::vector<float>* shower_likelihood = 0;
  std::vector<float>* dist2truescevtx = 0;
  in->SetBranchAddress( "npix_on_cosmic_v", &npix_on_cosmic_v );
  in->SetBranchAddress( "attcluster_nshower_v",  &attcluster_nshower_v );
  in->SetBranchAddress( "attcluster_ntrack_v",  &attcluster_ntrack_v );
  in->SetBranchAddress( "attcluster_nall_v",  &attcluster_nall_v );
  in->SetBranchAddress( "dist2truescevtx_v", &dist2truescevtx );
  in->SetBranchAddress( "shower_likelihood_v", &shower_likelihood );
  in->SetBranchAddress( "nbackwards_shower_pts_v", &nbackwards_shower_pts );
  in->SetBranchAddress( "nbackwards_track_pts_v",  &nbackwards_track_pts );
  in->SetBranchAddress( "nforwards_track_pts_v",  &nforwards_track_pts );    
                       
  int nentries = in->GetEntries();
  
  TFile* out = new TFile("plots_vareco_ana.root","recreate");
  
  
  const int nsamples = 3;
  enum { k1eVA=0, k1e1p, kAll };
  std::vector<std::string> sample_names = {"is1eVA","1e1p","all"};
  std::vector<std::string> sample_cuts = {"is1l0p0pi==1 && evis_had>30.0",
                                          "is1l1p0pi==1",
                                          "1==1"};

  // provides way to decide on remaining VA candidates per event
  // sort on shower-likelihood
  struct EventVA_t {
    int index;
    float shower_likelihood;
    bool operator<( EventVA_t& rhs ) {
      if ( shower_likelihood<rhs.shower_likelihood )
        return true;
      return false;
    }
  };
  
  // Efficiency versus Enu
  TH1D* henu[nsamples] = {0};
  TH1D* henu_pass[nsamples] = {0};
  TH1D* henu_eff[nsamples] = {0};      
  for (int i=0; i<nsamples; i++) {
    std::stringstream ss;
    ss << "hEnu_" << sample_names[i];
    henu[i] = new TH1D(ss.str().c_str(),";true E_{#nu} (MeV)", 50,0,500);
    ss << "_pass";
    henu_pass[i] = new TH1D(ss.str().c_str(),";true E_{#nu} (MeV)",50,0,500);

    std::stringstream sseff;
    sseff << "hEff_" << sample_names[i];
    henu_eff[i] = new TH1D(sseff.str().c_str(),";true E_{#nu} (MeV)",50,0,500);
  }

  // Distance to true vertex
  TH1D* hdist[nsamples] = {0};
  for (int i=0; i<nsamples; i++) {
    std::stringstream ss;
    ss << "hdist_" << sample_names[i];
    hdist[i] = new TH1D(ss.str().c_str(),";cm;",500,0,500);
  }
    
  for (int ientry=0; ientry<nentries; ientry++) {

    std::cout << "[ ENTRY " << ientry << "]" << std::endl;
    
    in->GetEntry(ientry);
    
    // define cut
    // && attcluster_nall_v>10 && (attcluster_nshower_v/attcluster_nall_v<0.95 || attcluster_nshower_v>400) && attcluster_ntrack_v<400
    int nvacand = (int)shower_likelihood->size();
    std::vector<EventVA_t> passing_va;
    for (int iva=0; iva<nvacand; iva++) {
      bool pass_cosmic = ((*npix_on_cosmic_v)[iva]<10);
      bool pass_nall   = ((*attcluster_nall_v)[iva]>10);
      float shower_frac = ((float)(*attcluster_nshower_v)[iva]/(float)(*attcluster_nall_v)[iva]);
      bool pass_nshower = shower_frac<0.95 || (*attcluster_nshower_v)[iva]>500;
      bool pass_ntrack  = (*attcluster_ntrack_v)[iva]<400;
      
      // std::cout << " pass cosmic: " << pass_cosmic << std::endl;
      // std::cout << " pass nall: " << pass_nall << std::endl;
      // std::cout << " pass shower: " << pass_nshower << std::endl;
      // std::cout << " pass track: " << pass_ntrack << std::endl;
      
      bool pass_cand = pass_cosmic & pass_nall & pass_nshower & pass_ntrack;
      if ( pass_cand ) {
        EventVA_t ev;
        ev.index = iva;
        ev.shower_likelihood = (*shower_likelihood)[iva];
        passing_va.push_back(ev);
      }
    }

    std::sort( passing_va.begin(), passing_va.end() );
    bool pass = (passing_va.size()>0) ? true : false;

    float minll_dist2truevtx = -1.0;
    bool true_positive = false;
    if ( pass ) {
      minll_dist2truevtx = (*dist2truescevtx)[ passing_va.front().index ];
      if ( minll_dist2truevtx<5.0 )
        true_positive = true;
    }
    
    // fill sample hists
    henu[kAll]->Fill( Enu_true );
    if ( true_positive ) {
      henu_pass[kAll]->Fill( Enu_true );
      henu_eff[kAll]->Fill( Enu_true );      
    }
    if ( pass ) {
      hdist[kAll]->Fill( minll_dist2truevtx );
    }

    // 1e0p
    if ( is1l0p0pi==0 && evis_had>30.0 ) {
      henu[k1eVA]->Fill( Enu_true );
      if ( pass ) {
        hdist[k1eVA]->Fill( minll_dist2truevtx );
      }
      if ( true_positive ) {
        henu_pass[k1eVA]->Fill( Enu_true );
        henu_eff[k1eVA]->Fill( Enu_true );              
      }
    }


    if ( is1l1p0pi==1 ) {
      henu[k1e1p]->Fill( Enu_true );
      if ( pass ) {
        hdist[k1e1p]->Fill( minll_dist2truevtx );
      }
      if ( true_positive ) {
        henu_pass[k1e1p]->Fill( Enu_true );
        henu_eff[k1e1p]->Fill( Enu_true );        
      }
    }
  }

  // efficiency
  for (int i=0; i<nsamples; i++)
    henu_eff[i]->Divide( henu[i] );

  // Stacks

  // bad rate
  int cutoff = hdist[kAll]->GetXaxis()->FindBin(5.0);
  int nbins = hdist[kAll]->GetXaxis()->GetNbins();
  float badrate = hdist[kAll]->Integral(cutoff+1,nbins+1)/henu[kAll]->Integral();
  std::cout << "Bad rate: " << badrate << " >5 cm false positive per event" << std::endl;
  
  cutoff = hdist[kAll]->GetXaxis()->FindBin(20.0);  
  badrate = hdist[kAll]->Integral(cutoff+1,nbins+1)/henu[kAll]->Integral();
  std::cout << "Bad rate: " << badrate << " >20 cm false positive per event" << std::endl;
  
  out->Write();
  
  return 0;
}
