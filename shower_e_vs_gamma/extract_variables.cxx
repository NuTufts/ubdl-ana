#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include "larflow/Reco/NuSelectionVariables.h"
#include "larflow/Reco/NuVertexCandidate.h"
#include "ublarcvapp/ubdllee/dwall.h"

#include "ShowerLikelihoodBuilder.h"

int main( int nargs, char** argv ) {

  std::cout << "Apply Gen-2 selection and save isolated events" << std::endl;
  std::cout << "KPSRECO: "<< argv[1] << std::endl;
  std::cout << "LARCV: " << argv[2] << std::endl;
  std::cout << "OUTFILE: " << argv[3] << std::endl;  
  
  std::string input_kpsreco   = argv[1];
  std::string input_larcv     = argv[2];
  std::string output_filename = argv[3];

  bool DEBUG = true;
  bool HAS_MC = true;

  // source files
  TFile* kpsreco_file = new TFile( input_kpsreco.c_str(), "open" );

  // KPRECO TREE
  //TTree* kprecotree = (TTree*)kpsreco_file->Get("KPSRecoManagerTree");
  TTree* kprecotree = (TTree*)kpsreco_file->Get("perfect_reco");  
  std::vector< larflow::reco::NuVertexCandidate >* pnu_perfect_v = nullptr;
  kprecotree->SetBranchAddress( "nu_perfect_v", &pnu_perfect_v );

  // IO managers
  larcv::IOManager iolcv( larcv::IOManager::kREAD );
  iolcv.add_in_file( input_larcv );
  iolcv.initialize();
    
  int nentries = kprecotree->GetEntries();

  // output file
  TFile* out = new TFile(output_filename.c_str(),"new");
  larflow::reco::ShowerLikelihoodBuilder builder;

  // Analysis classes

  // shower likelihood

  for (int ientry=0; ientry<nentries; ientry++) {

    std::cout << "//////// ENTRY " << ientry << "  //////////////" << std::endl;
    
    // load entry
    kprecotree->GetEntry(ientry);
    iolcv.read_entry(ientry);

    larcv::EventImage2D* ev_img
      = (larcv::EventImage2D*)iolcv.get_data(larcv::kProductImage2D,"wire");
    auto const& adc_v = ev_img->as_vector();
    
    // get perfect shower reco
    if ( pnu_perfect_v->size()==1 ) {
      auto const& nuperfect = pnu_perfect_v->at(0);
      builder.fill( nuperfect.shower_v, nuperfect.shower_trunk_v, adc_v );
    }

    // add data into shower likelihood class

    
  }


  // save shower likelihood class
  out->Write();
  
  return 0;
}
