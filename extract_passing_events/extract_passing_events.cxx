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

int main( int nargs, char** argv ) {

  std::cout << "Apply Gen-2 selection and save isolated events" << std::endl;
  std::cout << "KPSRECO: "<< argv[1] << std::endl;
  std::cout << "LARCV: " << argv[2] << std::endl;
  std::cout << "KPS LARLITE: " << argv[3] << std::endl;
  std::cout << "OUTFILE: " << argv[4] << std::endl;  
  
  std::string input_kpsreco   = argv[1];
  std::string input_larcv     = argv[2];
  std::string input_larlite   = argv[3];
  std::string output_filename = argv[4];

  bool DEBUG = true;
  

  TFile* larcv_file   = new TFile( input_larcv.c_str(), "open" );
  TFile* kpsreco_file = new TFile( input_kpsreco.c_str(), "open" );
  TFile* larlite_file = new TFile( input_larlite.c_str(), "open" );

  // source trees
  TTree* kprecotree = (TTree*)kpsreco_file->Get("KPSRecoManagerTree");
  TTree* img2dtree  = (TTree*)larcv_file->Get("image2d_wire_tree");
  TTree* thrumutree = (TTree*)larcv_file->Get("image2d_thrumu_tree");
  TTree* larliteid  = (TTree*)larlite_file->Get("larlite_id_tree");
  std::vector< TTree* > cosmic_tree_v;
  std::vector< TTree* > cosmic_outtree_v;
  std::vector< std::string > cosmic_tree_name_v
    = {"boundarycosmicnoshift",
       "containedcosmic" };
  for ( auto& name : cosmic_tree_name_v ) {
    std::stringstream treename;
    treename << "track_" << name << "_tree";
    TTree* cosmic_track = (TTree*)larlite_file->Get( treename.str().c_str() );
    TTree* copy_track   = (TTree*)cosmic_track->CloneTree(0);
    cosmic_tree_v.push_back( cosmic_track );
  }
  std::vector< TTree* > larlite_mcinfo_tree_v;
  std::vector< std::string > mcinfo_tree_name_v
    = {"mctruth_generator_tree",
       "mcshower_mcreco_tree",
       "mctrack_mcreco_tree" };
  for ( auto& name : mcinfo_tree_name_v ) {
    TTree* lltree = (TTree*)larlite_file->Get( name.c_str() );
    larlite_mcinfo_tree_v.push_back(lltree);
  }

  // output file and cloned trees
  TFile* output_file  = new TFile( output_filename.c_str(), "new" );  
  TTree* out_kpreco_tree = (TTree*)kprecotree->CloneTree(0);
  TTree* out_img2d_tree  = (TTree*)img2dtree->CloneTree(0);
  TTree* out_thrumu_tree = (TTree*)thrumutree->CloneTree(0);
  TTree* out_llid_tree   = (TTree*)larliteid->CloneTree(0);
  std::vector< TTree* > out_mcinfo_v;

  for ( auto& ptree : cosmic_tree_v ) {
    TTree* copy_track = (TTree*)ptree->CloneTree(0);
    cosmic_outtree_v.push_back( copy_track );
  }
  for ( auto& ptree : larlite_mcinfo_tree_v )  {
    TTree* copy_tree = (TTree*)ptree->CloneTree(0);
    out_mcinfo_v.push_back( copy_tree );
  }


  // cut stages
  enum { kFV=0,           // [0] true vertex in FV (10 cm from TPC boundary), sets baseline for efficiency study
         kVertexCand3cm,  // [1] reco candidate formed within 3 cm of vertex
         kMinShowerSize,  // [2] min shower size cut
         kNShowerProngs,  // [3] number of shower prongs
         kNTrackProngs,   // [4] number of track prongs         
         kHadronic,       // [5] see hadronic particles (proton or vertex activity)
         kShowerGap,      // [6] shower gap
         kVertexAct,      // [7] vertex activity cut
         kShowerLLCut,    // [8] shower likelihood cut
         kRecoFV,         // [9] reco fv cut
         kAllCuts,        // [10] All cuts applied except FV -- represents reco pass rate
         kNumCuts };      // [11] Number in enum

  std::vector<std::string> selcut_names
    = { "fv",             // [0]
        "vertexcand",     // [1]
        "minshower",      // [2]
        "nshowerprongs",  // [3]
        "ntrackprongs",   // [4]
        "hadronic",       // [5]
        "showergap",      // [6]
        "vertexact",      // [7]
        "showerll",       // [8]
        "recofv",         // [9]
        "allreco",        // [10]
        "numcuts"};       // [11]
  
  
  const float vtx_cutoff  = 5.0;
  const float vtx_cutoff2 = 20.0;

  const int cut_cosmic_npixels = 3;
  
  // true MC variables
  int is1l1p0pi;
  int is1l0p0pi;
  float evis_lep;
  float evis_had;
  float Enu_true;
  float vtx_dwall;
  int genie_mode;
  int interactionType;
  int ccnc;
  int nu_pdg;
  kprecotree->SetBranchAddress( "Enu_true",  &Enu_true );
  kprecotree->SetBranchAddress( "evis_lep",  &evis_lep );  
  kprecotree->SetBranchAddress( "evis_had",  &evis_had );
  kprecotree->SetBranchAddress( "is1l1p0pi", &is1l1p0pi );
  kprecotree->SetBranchAddress( "is1l0p0pi", &is1l0p0pi );
  kprecotree->SetBranchAddress( "vtx_dwall", &vtx_dwall );
  kprecotree->SetBranchAddress( "currentType", &ccnc );
  kprecotree->SetBranchAddress( "genieMode", &genie_mode );
  kprecotree->SetBranchAddress( "interactionType", &interactionType );
  kprecotree->SetBranchAddress( "nu_pdg", &nu_pdg );
  
  // Selection variables
  /**
     max_shower_nhits>500 
     && nshowers>0 && nshowers<=2 && ntracks<=2 
     && ( (max_proton_pid<0) || (vertex_hip_fraction>0.5)) 
     && min_shower_gap<5.0 && max_shower_gap<5.0 
     && (max_track_length>3.0 || vertex_charge_per_pixel>50.0)
  */
  std::vector< larflow::reco::NuSelectionVariables >* pnu_sel_v = nullptr;
  std::vector< larflow::reco::NuVertexCandidate >* pnu_fitted_v = nullptr;  
  kprecotree->SetBranchAddress( "nu_sel_v",   &pnu_sel_v );    ///< selection variables per vertex
  kprecotree->SetBranchAddress( "nufitted_v", &pnu_fitted_v ); ///< reconstruction objects per vertex

  
  int nentries = kprecotree->GetEntries();

  int TOTAL_PASSING = 0;
  
  for (int ientry=0; ientry<nentries; ientry++) {

    if ( ientry%100==0 )
      std::cout << "[ ENTRY " << ientry << "]" << std::endl;
    
    kprecotree->GetEntry(ientry);
    img2dtree->GetEntry(ientry);
    thrumutree->GetEntry(ientry);
    larliteid->GetEntry(ientry);
    for ( auto& ptree : cosmic_tree_v )
      ptree->GetEntry(ientry);
    for ( auto& ptree : larlite_mcinfo_tree_v )
      ptree->GetEntry(ientry);

    std::vector<larflow::reco::NuVertexCandidate>    pass_vtx_v;
    std::vector<larflow::reco::NuSelectionVariables> pass_var_v;
    int num_pass = 0;
    
    // truth cuts
    bool cut_fv = vtx_dwall>10.0;

    std::vector<bool> event_passes_cut( kNumCuts, false ); //< indicates if event has a vertex that passes this cut stage
    int nvtx = (int)(*pnu_sel_v).size();
    for (int ivtx=0; ivtx<nvtx; ivtx++) {

      auto& nusel = (*pnu_sel_v)[ivtx];
      auto& nuvtx = (*pnu_fitted_v)[ivtx];

      // dwall-reco
      int reco_boundary = 0;
      float reco_dwall = ublarcvapp::dwall( nuvtx.pos, reco_boundary );

      // selection cuts
      std::vector<bool> vtx_pass( kNumCuts, false );
      vtx_pass[kFV] = cut_fv;
      vtx_pass[kVertexCand3cm] = nusel.dist2truevtx<3.0;
      vtx_pass[kMinShowerSize] = nusel.max_shower_nhits>500;
      vtx_pass[kNShowerProngs] = ( nusel.nshowers>0 && nusel.nshowers<=2 );
      vtx_pass[kNTrackProngs]  = ( nusel.ntracks<=2 );
      vtx_pass[kHadronic]      = (nusel.max_proton_pid<0 || nusel.vertex_hip_fraction>0.5);
      vtx_pass[kShowerGap]     = (nusel.min_shower_gap<2.0 && nusel.max_shower_gap<2.0);
      vtx_pass[kVertexAct]     = (nusel.max_track_length>3.0 || nusel.vertex_charge_per_pixel>50.0);
      //vtx_pass[kShowerLLCut]   = (nusel.largest_shower_ll < 0.0 || nusel.closest_shower_ll < 0.0 );
      //vtx_pass[kRecoFV]        = (reco_dwall>5.0);
      vtx_pass[kShowerLLCut]   = true;
      vtx_pass[kRecoFV]        = true;
      vtx_pass[kAllCuts]       = true;

      // reco variable cuts only
      for ( int i=kMinShowerSize; i<kAllCuts; i++)
        vtx_pass[kAllCuts] = vtx_pass[kAllCuts] && vtx_pass[i];

      // for debug
      if (DEBUG) {
        std::cout << "[entry " << ientry << ", vtx" << ivtx << "] " << std::endl;
        bool vtx_seq = true;
        for (int i=0; i<=kAllCuts; i++) {
          vtx_seq = vtx_seq && vtx_pass[i]; // follows sequence        
          std::cout << "  " << selcut_names[i] << ": this=" << vtx_pass[i] << " chain=" << vtx_seq << std::endl;
        }
      }

      if ( vtx_pass[kAllCuts] ) {

        pass_vtx_v.emplace_back( std::move( nuvtx ) );
        pass_var_v.emplace_back( std::move( nusel ) );
        num_pass++;
        TOTAL_PASSING++;
      }
      
    }//end of VERTEX LOOP


    std::cout << " Number vertices passing: " << num_pass << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    
    if ( num_pass>0 ) {

      // clear the vector
      pnu_sel_v->clear();
      pnu_fitted_v->clear();

      for ( size_t ivtx=0; ivtx<pass_vtx_v.size(); ivtx++ ) {
        auto& pass_nuvtx = pass_vtx_v[ivtx];
        auto& pass_nusel = pass_var_v[ivtx];

        pnu_fitted_v->emplace_back( std::move(pass_nuvtx) );        
        pnu_sel_v->emplace_back( std::move(pass_nusel) );

      }

      out_kpreco_tree->Fill();
      out_img2d_tree->Fill();
      out_thrumu_tree->Fill();
      out_llid_tree->Fill();
      for ( auto& ptree : cosmic_outtree_v )
        ptree->Fill();
      for ( auto& ptree : out_mcinfo_v )
        ptree->Fill();
    }

  }//end of entry loop

  std::cout << "TOTAL PASSING: " << TOTAL_PASSING << std::endl;
  
  out_kpreco_tree->AutoSave();
  out_img2d_tree->AutoSave();
  out_thrumu_tree->AutoSave();
  out_llid_tree->AutoSave();
  for ( auto& ptree : cosmic_outtree_v ) { 
    ptree->AutoSave();
    delete ptree;
    ptree = nullptr;
  }
  for ( auto& ptree : out_mcinfo_v ) {
    ptree->AutoSave();
    delete ptree;
    ptree = nullptr;
  }
  
  delete out_kpreco_tree;
  delete out_img2d_tree;
  delete out_thrumu_tree;
  delete out_llid_tree;
  
  return 0;
}
