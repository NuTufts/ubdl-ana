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

int main( int nargs, char** argv ) {

  std::cout << "Gen-2 1e1p Selection Plots" << std::endl;

  std::ifstream input_list(argv[1]);
  char zinput[1028];  
  input_list >> zinput;

  int is_bnbnu = std::atoi(argv[2]);
  
  std::vector<std::string> input_v;
  do {
    input_v.push_back( std::string(zinput) );
    input_list >> zinput;
  } while ( input_list.good() && !input_list.eof() );
  
  TChain* in = new TChain("KPSRecoManagerTree");
  int ninputs = 0;
  for ( auto& finput : input_v ) {
    in->Add( finput.c_str() );
    std::cout << finput << std::endl;
    ninputs++;
  }
  std::cout << "number of files loaded: " << ninputs << std::endl;

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
  in->SetBranchAddress( "Enu_true",  &Enu_true );
  in->SetBranchAddress( "evis_lep",  &evis_lep );  
  in->SetBranchAddress( "evis_had",  &evis_had );
  in->SetBranchAddress( "is1l1p0pi", &is1l1p0pi );
  in->SetBranchAddress( "is1l0p0pi", &is1l0p0pi );
  in->SetBranchAddress( "vtx_dwall", &vtx_dwall );
  in->SetBranchAddress( "currentType", &ccnc );
  in->SetBranchAddress( "genieMode", &genie_mode );
  in->SetBranchAddress( "interactionType", &interactionType );
  in->SetBranchAddress( "nu_pdg", &nu_pdg );
  
  // Selection variables
  /**
     max_shower_nhits>500 
     && nshowers>0 && nshowers<=2 && ntracks<=2 
     && ( (max_proton_pid<0) || (vertex_hip_fraction>0.5)) 
     && min_shower_gap<5.0 && max_shower_gap<5.0 
     && (max_track_length>3.0 || vertex_charge_per_pixel>50.0)
  */
  // int max_shower_nhits;
  // int nshowers;
  // int ntracks;
  // float max_proton_pid;
  // float vertex_hip_fraction;
  // float min_shower_gap;
  // float max_shower_gap;
  // float max_track_length;
  // float vertex_charge_per_pixel;
  std::vector< larflow::reco::NuSelectionVariables >* pnu_sel_v = nullptr;
  
  // good ol' c++ boilerplate
  // in->SetBranchAddress( "max_shower_nhits", &max_shower_nhits );
  // in->SetBranchAddress( "nshowers", &nshowers );
  // in->SetBranchAddress( "ntracks", &ntracks );
  // in->SetBranchAddress( "max_proton_pid", &max_proton_pid );
  // in->SetBranchAddress( "vertex_hip_fraction", &vertex_hip_fraction );
  // in->SetBranchAddress( "min_shower_gap", &min_shower_gap );
  // in->SetBranchAddress( "max_shower_gap", &max_shower_gap );
  // in->SetBranchAddress( "max_track_length", &max_track_length );
  // in->SetBranchAddress( "vertex_charge_per_pixel", &vertex_charge_per_pixel );
  in->SetBranchAddress( "nu_sel_v", &pnu_sel_v ); ///< selection variables per vertex

  
  int nentries = in->GetEntries();
  
  TFile* out = new TFile("plots_1e1p_sel.root","recreate");

  // truth categories
  const int nsamples = 3;
  enum { k1eVA=0, k1e1p, kAll };
  std::vector<std::string> sample_names = {"is1eVA","1e1p","all"};
  std::vector<std::string> sample_cuts  = {"is1l0p0pi==1 && evis_had>30.0",
                                           "is1l1p0pi==1",
                                           "1==1"};

  // truth modes
  const int nmodes = 7;
  enum { kCCQE=0, kCCRes, kCCOther, kNCQE, kNCRes, kNCOther, kAllModes, kNumModes };
  std::vector<std::string> modename
    = { "ccqe",
        "ccres",
        "ccother",
        "ncqe",
        "ncres",
        "ncother",
        "all" };

  // true reco state
  enum { kOnVertex=0, //on nu + within 3 cm of true vertex
         kOnNu, // on nu + > 3 cm of true vertex
         kOffNu, // not on nu pixels
         kNumRecoStatus
  };
  std::vector<std::string> reco_status_names
    = {"onvtx",
       "offvtxnu",
       "offvtxcosmic"};
  
  // cut stages
  enum { kFV=0,           // true vertex in FV (10 cm from TPC boundary), sets baseline for efficiency study
         kVertexCand3cm,  // reco candidate formed within 3 cm of vertex
         kMinShowerSize,  // min shower size cut
         kNShowerProngs,  // number of shower prongs
         kNTrackProngs,   // number of track prongs         
         kHadronic,       // see hadronic particles (proton or vertex activity)
         kShowerGap,      // shower gap
         kVertexAct,      // vertex activity cut
         kAllCuts,        // All cuts applied except FV -- represents reco pass rate
         kNumCuts };      // Number in enum
  std::vector<std::string> selcut_names
    = { "fv",
        "vertexcand",
        "minshower",
        "nshowerprongs",
        "ntrackprongs",        
        "hadronic",
        "showergap",
        "vertexact",
        "allreco",
        "numcuts"};

  // Cut variables for studying optimal cuts
  enum { kdwall=0,
         kdist2true,
         kmaxshowerhits,
         knshowerprongs,
         kntrackprongs,
         kllpid,
         khipfraction,
         kminshowergap,
         kmaxshowergap,
         kmaxtracklen,
         kvertexact,
         kNumCutVariables };
         
  std::vector<std::string> cutvar_names
    = { "dwall",
        "dist2true",
        "nmaxshowerhits",
        "nshowerprongs",
        "ntrackprongs",
        "llpid",
        "hipfraction",
        "minshowergap",
        "maxshowergap",
        "maxtracklen",
        "vertexcharge" };
  float cutvar_range[11][2] = { {-10,200},  // dwall
                                {0, 50 },   // distance to true vertex
                                {0, 2000},  // hits in largest shower
                                {0, 10},    // num shower prongs
                                {0, 10},    // num track prongs
                                {-100,100}, // proton likelihood
                                {0,1.01},   // hip fraction
                                {0,50.0},   // minshowergap
                                {0,50.0},   // maxshowergap
                                {0,500},    // maxtracklen
                                {0,150.0}   // vertex activity: charge per pixel around reco vertex
  };
  int cutvar_nbins[11] = { 210, // dwall
                           150, // dist 2 true
                           100, // hits in largest shower
                           10,  // num shower prongs
                           10,  // num track prongs
                           200, // proton likelihood
                           50,  // hip fraction
                           50,  // minshowergap
                           50,  // maxshowergap
                           500, // maxtracklen
                           50 };// vertex activity

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

  struct EventDist2True_t {
    int index;
    float dist;
    EventDist2True_t()
      : index(-1),
        dist(1000.0)
    {};
    EventDist2True_t( float d, int idx )
      : index(idx),
        dist(d)
    {};    
    bool operator<( EventDist2True_t& rhs ) {
      if ( dist<rhs.dist )
        return true;
      return false;
    };
  };
  
  // Efficiency versus Enu
  TH1D* henu[nsamples][kNumCuts][kNumModes]      = {0};
  TH1D* henu_eff[nsamples][kNumCuts][kNumModes]  = {0};      
  for (int isample=0; isample<nsamples; isample++) {
    for (int icut=0; icut<kNumCuts; icut++) {
      for (int imode=0; imode<kNumModes; imode++) {
        std::stringstream ss;
        ss << "hEnu_" << sample_names[isample] << "_" << selcut_names[icut] << "cut" << "_" << modename[imode];
        henu[isample][icut][imode] = new TH1D(ss.str().c_str(),";true E_{#nu} (MeV)", 30,0,3000);
        
        std::stringstream sseff;
        sseff << "hEff_" << sample_names[isample] << "_" << selcut_names[icut] << "cut" << "_" << modename[imode];
        henu_eff[isample][icut][imode] = new TH1D(sseff.str().c_str(),";true E_{#nu} (MeV)",30,0,3000);
      }
    }
  }

  // Distance to true vertex
  TH1D* hdist[nsamples] = {0};
  for (int i=0; i<nsamples; i++) {
    std::stringstream ss;
    ss << "hdist_" << sample_names[i];
    hdist[i] = new TH1D(ss.str().c_str(),";cm;",500,0,500);
  }

  // number of shower pixels forward: proxy for "reco energy" until we build one
  TH1D* hnshower[nsamples] = {0};
  TH2D* hshower_vs_enu[nsamples] = {0};
  TH2D* hshower_vs_evislep[nsamples] = {0};
  for (int isample=0; isample<nsamples; isample++) {
    std::stringstream ss;
    ss << "hnshower_" << sample_names[isample];
    hnshower[isample] = new TH1D(ss.str().c_str(),";num shower hits;",100,0,10000);

    std::stringstream ss2d;
    ss2d << "hnshower_vs_enu_" << sample_names[isample];
    hshower_vs_enu[isample] = new TH2D(ss2d.str().c_str(), "",30,0,3000,100,0,20000);

    std::stringstream ss_vs_elep;
    ss_vs_elep << "hnshower_vs_evislep_" << sample_names[isample];
    hshower_vs_evislep[isample] = new TH2D(ss_vs_elep.str().c_str(), "",150,0,3000,200,0,10000);
  }
                                     
  // cut variables, with truth tag. for study selection power.
  std::vector<TH1D*> hvariable_bad_v( cutvar_names.size(), 0 );
  std::vector<TH1D*> hvariable_good_v( cutvar_names.size(), 0 );
  for (int icut=0; icut<kNumCutVariables; icut++) {
    std::stringstream ss_bad;
    ss_bad << "hCutVar_" << cutvar_names[icut] << "_bad";
    std::stringstream ss_good;
    ss_good << "hCutVar_" << cutvar_names[icut] << "_good";
    hvariable_bad_v[icut]  = new TH1D(ss_bad.str().c_str(),"",
                                      cutvar_nbins[icut], cutvar_range[icut][0], cutvar_range[icut][1] );
    hvariable_good_v[icut] = new TH1D(ss_good.str().c_str(),";",
                                      cutvar_nbins[icut], cutvar_range[icut][0], cutvar_range[icut][1] );

  }
  
  // plot cut variables that currently pass all cut, to see if on-nu or off-nu
  // this is to inform directions for reco improvement
  TH1D* hvar_onnu[kNumRecoStatus][kNumCutVariables] = {0};
  for (int icut=0; icut<kNumCutVariables; icut++) {  
    for (int istat=0; istat<(int)kNumRecoStatus; istat++) {
      std::stringstream ss;
      ss << "hCutVar_" << cutvar_names[icut] << "_" << reco_status_names[istat];
      hvar_onnu[istat][icut] = new TH1D(ss.str().c_str(),"", cutvar_nbins[icut],
                                        cutvar_range[icut][0], cutvar_range[icut][1] );
    }
  }
  
  for (int ientry=0; ientry<nentries; ientry++) {
  //for (int ientry=0; ientry<100; ientry++) {  

    if ( ientry%100==0 )
      std::cout << "[ ENTRY " << ientry << "]" << std::endl;
    
    in->GetEntry(ientry);

    // truth cuts
    bool cut_fv = vtx_dwall>10.0;

    // which mode are we
    int event_mode = 0;
    if ( ccnc==0 ) {
      // charged current
      if ( interactionType==1001 )
        event_mode = kCCQE;
      else if (  interactionType>=1003 && interactionType<=1090 )
        event_mode = kCCRes;
      else
        event_mode = kCCOther;
    }
    else {
      if ( interactionType==1002 )
        event_mode = kNCQE;
      else if (  interactionType>=1003 && interactionType<=1090 )
        event_mode = kNCRes;
      else
        event_mode = kNCOther;
    }

    if ( is_bnbnu==1 && abs(nu_pdg)==12 && ccnc==0 )
      continue; // omit nue-CC events from bnb nu histograms
      
    
    // find best reco vertex at each cut stage, measured by closeness to true vertex
    std::vector<EventDist2True_t> index_by_dist_v;
    std::vector<bool> event_passes_cut( kNumCuts, false ); //< indicates if event has a vertex that passes this cut stage
    int nvtx = (int)(*pnu_sel_v).size();
    for (int ivtx=0; ivtx<nvtx; ivtx++) {

      auto const& nusel = (*pnu_sel_v)[ivtx];

      EventDist2True_t idx( nusel.dist2truevtx, ivtx );

      // selection cuts
      std::vector<bool> vtx_pass( kNumCuts, false );
      vtx_pass[kFV] = cut_fv;
      vtx_pass[kVertexCand3cm] = nusel.dist2truevtx<3.0;
      vtx_pass[kMinShowerSize] = nusel.max_shower_nhits>500;
      vtx_pass[kNShowerProngs] = ( nusel.nshowers>0 && nusel.nshowers<=2 );
      vtx_pass[kNTrackProngs]  = ( nusel.ntracks<=2 );
      vtx_pass[kHadronic]      = (nusel.max_proton_pid<0 || nusel.vertex_hip_fraction>0.5);
      vtx_pass[kShowerGap]     = (nusel.min_shower_gap<5.0 && nusel.max_shower_gap<5.0);
      vtx_pass[kVertexAct]     = (nusel.max_track_length>3.0 || nusel.vertex_charge_per_pixel>50.0);
      vtx_pass[kAllCuts]       = true;

      // reco variable cuts only
      for ( int i=kMinShowerSize; i<=kVertexAct; i++) 
        vtx_pass[kAllCuts] = vtx_pass[kAllCuts] && vtx_pass[i];

      bool vtx_seq = true;
      for (int icut=0; icut<kNumCuts; icut++) {
        vtx_seq = vtx_seq && vtx_pass[icut];
        event_passes_cut[ icut ] = event_passes_cut[icut] || vtx_seq;
      }

      if ( vtx_pass[kAllCuts] )
        index_by_dist_v.push_back( idx );

      // for debug
      // std::cout << "[entry " << ientry << ", vtx" << ivtx << "] " << std::endl;
      // for (int i=0; i<=kAllCuts; i++) {
      //   std::cout << "  " << selcut_names[i] << ": " << vtx_pass[i] << std::endl;
      // }

      // Cut variables: study between "good" or "bad" vertex
      if ( nusel.dist2truevtx<3.0 ) {
        
        hvariable_good_v[kdwall]->Fill( vtx_dwall );
        hvariable_good_v[kdist2true]->Fill( nusel.dist2truevtx );
        hvariable_good_v[kmaxshowerhits]->Fill( nusel.max_shower_nhits );
        hvariable_good_v[knshowerprongs]->Fill( nusel.nshowers );
        hvariable_good_v[kntrackprongs]->Fill( nusel.ntracks );
        hvariable_good_v[kllpid]->Fill( nusel.max_proton_pid );
        hvariable_good_v[khipfraction]->Fill( nusel.vertex_hip_fraction );
        hvariable_good_v[kminshowergap]->Fill( nusel.min_shower_gap );
        hvariable_good_v[kmaxshowergap]->Fill( nusel.max_shower_gap );
        hvariable_good_v[kmaxtracklen]->Fill( nusel.max_track_length );
        hvariable_good_v[kvertexact]->Fill( nusel.vertex_charge_per_pixel );
        
      }
      else {

        hvariable_bad_v[kdwall]->Fill( vtx_dwall );
        hvariable_bad_v[kdist2true]->Fill( nusel.dist2truevtx );
        hvariable_bad_v[kmaxshowerhits]->Fill( nusel.max_shower_nhits );
        hvariable_bad_v[knshowerprongs]->Fill( nusel.nshowers );
        hvariable_bad_v[kntrackprongs]->Fill( nusel.ntracks );
        hvariable_bad_v[kllpid]->Fill( nusel.max_proton_pid );
        hvariable_bad_v[khipfraction]->Fill( nusel.vertex_hip_fraction );
        hvariable_bad_v[kminshowergap]->Fill( nusel.min_shower_gap );
        hvariable_bad_v[kmaxshowergap]->Fill( nusel.max_shower_gap );
        hvariable_bad_v[kmaxtracklen]->Fill( nusel.max_track_length );
        hvariable_bad_v[kvertexact]->Fill( nusel.vertex_charge_per_pixel );
        
      }

      // Cut variables: study correlation between reco state:
      // on vertex, off-vertex nu, off-vertex cosmic
      // determine state
      int vtx_reco_state = 0;
      if ( nusel.dist2truevtx<3.0 ) {
        vtx_reco_state  = (int)kOnVertex;
      }
      else {
        if ( nusel.truth_vtxFracNu>0.65 )
          vtx_reco_state = (int)kOnNu;
        else
          vtx_reco_state = (int)kOffNu;
      }

      if ( vtx_pass[kAllCuts] ) {
        hvar_onnu[vtx_reco_state][kdwall]->Fill( vtx_dwall );
        hvar_onnu[vtx_reco_state][kdist2true]->Fill( nusel.dist2truevtx );
        hvar_onnu[vtx_reco_state][kmaxshowerhits]->Fill( nusel.max_shower_nhits );
        hvar_onnu[vtx_reco_state][knshowerprongs]->Fill( nusel.nshowers );
        hvar_onnu[vtx_reco_state][kntrackprongs]->Fill( nusel.ntracks );
        hvar_onnu[vtx_reco_state][kllpid]->Fill( nusel.max_proton_pid );
        hvar_onnu[vtx_reco_state][khipfraction]->Fill( nusel.vertex_hip_fraction );
        hvar_onnu[vtx_reco_state][kminshowergap]->Fill( nusel.min_shower_gap );
        hvar_onnu[vtx_reco_state][kmaxshowergap]->Fill( nusel.max_shower_gap );
        hvar_onnu[vtx_reco_state][kmaxtracklen]->Fill( nusel.max_track_length );
        hvar_onnu[vtx_reco_state][kvertexact]->Fill( nusel.vertex_charge_per_pixel );
      }
      
    }
    std::sort( index_by_dist_v.begin(), index_by_dist_v.end() );

    // for debug
    // std::cout << "[entry results] ------------" <<  std::endl;
    // for (int i=0; i<=kAllCuts; i++) {
    //   std::cout << "  " << selcut_names[i] << ": " << event_passes_cut[i] << std::endl;
    // }
    // std::cout << "----------------------------" << std::endl;
    
    // Event-based plots

    // Enu filled based on if vertex passes cut stage
    bool still_passing = true;
    for (int icut=0; icut<=(int)kVertexAct; icut++) {
      still_passing = still_passing & event_passes_cut[icut];

      if ( !still_passing )
        break;

      // 1eVA
      if ( is1l0p0pi==1 && evis_had>30.0 ) {
        henu[k1eVA][icut][kAllModes]->Fill( Enu_true );
        henu_eff[k1eVA][icut][kAllModes]->Fill( Enu_true );
        henu[k1eVA][icut][event_mode]->Fill( Enu_true );
        henu_eff[k1eVA][icut][event_mode]->Fill( Enu_true );        
      }

      // 1e1p
      if ( is1l1p0pi==1 ) {
        henu[k1e1p][icut][kAllModes]->Fill( Enu_true );
        henu_eff[k1e1p][icut][kAllModes]->Fill( Enu_true );
        henu[k1e1p][icut][event_mode]->Fill( Enu_true );
        henu_eff[k1e1p][icut][event_mode]->Fill( Enu_true );                
      }

      // All
      henu[kAll][icut][kAllModes]->Fill( Enu_true );
      henu_eff[kAll][icut][kAllModes]->Fill( Enu_true );
      henu[kAll][icut][event_mode]->Fill( Enu_true );
      henu_eff[kAll][icut][event_mode]->Fill( Enu_true );                      
    }

    // nshowerhits: temporary proxy for neutrino energy
    if ( index_by_dist_v.size()>0 ) {
      // this means event had passing vertex
      
      int best_passing_vtx_index = index_by_dist_v.front().index;
      float num_shr_hits = (*pnu_sel_v)[best_passing_vtx_index].max_shower_nhits;
      hnshower[kAll]->Fill( num_shr_hits );
      hshower_vs_enu[kAll]->Fill( Enu_true, num_shr_hits );
      hshower_vs_evislep[kAll]->Fill( evis_lep, num_shr_hits );
      henu[kAll][kAllCuts][kAllModes]->Fill( Enu_true );
      henu[kAll][kAllCuts][event_mode]->Fill( Enu_true );      
				 
      if ( is1l0p0pi==1 && evis_had>30.0 ) {
        hnshower[k1eVA]->Fill( num_shr_hits );
	hshower_vs_enu[k1eVA]->Fill( Enu_true, num_shr_hits );
	hshower_vs_evislep[k1eVA]->Fill( evis_lep, num_shr_hits );
        henu[k1eVA][kAllCuts][kAllModes]->Fill( Enu_true );        
        henu[k1eVA][kAllCuts][event_mode]->Fill( Enu_true );
      }
      
      if ( is1l1p0pi==1 ) {
        hnshower[k1e1p]->Fill( (*pnu_sel_v)[best_passing_vtx_index].max_shower_nhits );
	hshower_vs_enu[k1e1p]->Fill( Enu_true, num_shr_hits );
	hshower_vs_evislep[k1e1p]->Fill( evis_lep, num_shr_hits );
        henu[k1e1p][kAllCuts][kAllModes]->Fill( Enu_true );        
        henu[k1e1p][kAllCuts][event_mode]->Fill( Enu_true );        
      }

      
    }
    
    // std::cout << "[entry] to continue." << std::endl;
    // std::cin.get();
    
  }//end of entry loop

  // efficiency
  for (int i=0; i<nsamples; i++) {
    for (int j=0; j<kAllCuts; j++) {
      henu_eff[i][j][kAllModes]->Divide( henu[i][kFV][kAllModes] );
    }
  }
  
  out->Write();
  delete in;
  
  return 0;
}
