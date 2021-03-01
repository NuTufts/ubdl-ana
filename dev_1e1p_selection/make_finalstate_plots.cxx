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

  std::cout << "Gen-2 Final State Selection Plots" << std::endl;

  std::ifstream input_list(argv[1]);
  char zinput[1028];  
  input_list >> zinput;

  int is_bnbnu = std::atoi(argv[2]);

  bool DEBUG = false;
  
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
  int nproton_60mev;
  int nmeson_35mev;
  int npi0;
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
  in->SetBranchAddress( "nproton_60mev", &nproton_60mev );
  in->SetBranchAddress( "nmeson_35mev",  &nmeson_35mev );
  in->SetBranchAddress( "npi0",  &npi0 );
  
  // Selection variables
  std::vector< larflow::reco::NuSelectionVariables >* pnu_sel_v = nullptr;
  std::vector< larflow::reco::NuVertexCandidate >* pnu_fitted_v = nullptr;  
  
  in->SetBranchAddress( "nu_sel_v",   &pnu_sel_v );    ///< selection variables per vertex
  in->SetBranchAddress( "nufitted_v", &pnu_fitted_v ); ///< reconstruction objects per vertex

  
  int nentries = in->GetEntries();
  
  TFile* out = new TFile("plots_finalstates.root","recreate");

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

  // final states
  enum { k1eVA=0,
         k1e_1p_0pi_nopi0,  k1e_0p_1pi_nopi0, k1e_0p_0pi_wpi0, 
         k1e_Np_0pi_nopi0,  k1e_1p_Npi_nopi0, k1e_1p_0pi_wpi0,
         k1e_Np_Npi_wpi0,
         k1mVA, //[8]
         k1m_1p_0pi_nopi0,  k1m_0p_1pi_nopi0, k1m_0p_0pi_wpi0,
         k1m_Np_0pi_nopi0,  k1m_1p_Npi_nopi0, k1m_1p_0pi_wpi0,
         k1m_Np_Npi_wpi0,
         // NC // [16]
         knc_1p_0pi_nopi0,  knc_0p_1pi_nopi0, knc_0p_0pi_wpi0,
         knc_Np_0pi_nopi0,  knc_1p_Npi_nopi0, knc_1p_0pi_wpi0,
         knc_Np_Npi_wpi0,
         kOffNu,
         kNumFinalStates };

  std::vector<std::string> final_state_names
    = { "1eVA",
        "1e_1p_0pi_nopi0",  "1e_0p_1pi_nopi0", "1e_0p_0pi_wpi0",
        "1e_Np_0pi_nopi0",  "1e_1p_Npi_nopi0", "1e_1p_0pi_wpi0",
        "1e_Np_Npi_wpi0",
        "1mVA",
        "1m_1p_0pi_nopi0",  "1m_0p_1pi_nopi0", "1m_0p_0pi_wpi0",
        "1m_Np_0pi_nopi0",  "1m_1p_Npi_nopi0", "1m_1p_0pi_wpi0",
        "1m_Np_Npi_wpi0",
        // NC
        "nc_1p_0pi_nopi0",  "nc_0p_1pi_nopi0", "nc_0p_0pi_wpi0",
        "nc_Np_0pi_nopi0",  "nc_1p_Npi_nopi0", "nc_1p_0pi_wpi0",
        "nc_Np_Npi_wpi0",
        // cosmic
        "offnu" };
        
    
  // cut stages
  enum { kFV=0,           // [0] true vertex in FV (10 cm from TPC boundary), sets baseline for efficiency study
         kVertexCand3cm,  // [1] reco candidate formed within 3 cm of vertex
         kMinShowerSize,  // [2] min shower size cut (might want to loosen)
         kNShowerProngs,  // [3] number of shower prongs
         kNTrackProngs,   // [4] number of track prongs         
         kShowerGap,      // [5] shower gap
         kTrackGap,       // [6] track gap
         kMaxTrackLen,    // [7] max track len
         kSecondShower,   // [8] second shower size
         kVertexAct,      // [9] vertex activity cut
         kRecoFV,         // [10] reco fv cut
         kShowerLLCut,    // [11] shower likelihood cut         
         kWCPixel,        // [12] Wire-Cell pixel cut
         kHadronic,       // [13] see hadronic particles (proton or vertex activity)         
         kAllCuts,        // [14] All cuts applied except FV -- represents reco pass rate
         kNumCuts };      // [15] Number in enum
  std::vector<std::string> selcut_names
    = { "fv",             // [0]
        "vertexcand",     // [1]
        "minshower",      // [2]
        "nshowerprongs",  // [3]
        "ntrackprongs",   // [4]
        "showergap",      // [5]
        "trackgap",       // [6]
        "maxtracklen",    // [7]
        "secondshower",   // [8]
        "vertexact",      // [9]
        "showerll",       // [10]        
        "recofv",         // [11]
        "wcpixel",        // [12]        
        "hadronic",       // [13]        
        "allreco",        // [14]
        "numcuts"};       // [15]

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
  TH1D* henu[kNumFinalStates][kNumCuts][kNumModes]      = {0};
  for (int istate=0; istate<(int)final_state_names.size(); istate++) {
    for (int icut=0; icut<kNumCuts; icut++) {
      for (int imode=0; imode<kNumModes; imode++) {
        std::stringstream ss;
        ss << "hEnu_" << final_state_names[istate] << "_" << selcut_names[icut] << "cut" << "_" << modename[imode];
        henu[istate][icut][imode] = new TH1D(ss.str().c_str(),";true E_{#nu} (MeV)", 30,0,3000);        
      }
    }
  }
  
  for (int ientry=0; ientry<nentries; ientry++) {

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
    std::vector<bool>  cutpass_onnu( kNumCuts, false );
    std::vector<float> cutpass_mindist( kNumCuts, 1e6 );
      
    int nvtx = (int)(*pnu_sel_v).size();
    for (int ivtx=0; ivtx<nvtx; ivtx++) {

      auto const& nusel = (*pnu_sel_v)[ivtx];
      auto const& nuvtx = (*pnu_fitted_v)[ivtx];

      EventDist2True_t idx( nusel.dist2truevtx, ivtx );

      // dwall-reco
      int reco_boundary = 0;
      float reco_dwall = ublarcvapp::dwall( nuvtx.pos, reco_boundary );

      // second shower size
      int nhits_second_shower = 0;
      if ( nuvtx.shower_v.size()>1 ) {
        std::vector<int> nhit_shower_v(nuvtx.shower_v.size(),0);
        for (size_t ishr=0; ishr<nuvtx.shower_v.size(); ishr++)
          nhit_shower_v[ishr] = (int)nuvtx.shower_v[ishr].size();
        std::sort( nhit_shower_v.begin(), nhit_shower_v.end() );
        nhits_second_shower = nhit_shower_v[1];
      }

      // selection cuts
      std::vector<bool> vtx_pass( kNumCuts, false );
      vtx_pass[kFV] = cut_fv; // [0]
      vtx_pass[kVertexCand3cm] = nusel.dist2truevtx<3.0; // [1]
      vtx_pass[kMinShowerSize] = nusel.max_shower_nhits>500; // [2]
      vtx_pass[kNShowerProngs] = ( nusel.nshowers>0 && nusel.nshowers<=2 ); // [3]
      vtx_pass[kNTrackProngs]  = ( nusel.ntracks<=2 ); // [4]
      vtx_pass[kShowerGap]     = nusel.nplanes_connected>=2; // [5]
      vtx_pass[kTrackGap]      = (nusel.ntracks==0 || nusel.min_track_gap<3.0); // [6]
      vtx_pass[kMaxTrackLen]   = (nusel.ntracks==0 || nusel.max_track_length<300.0); // [7]
      vtx_pass[kSecondShower]  = (nhits_second_shower<100); // [8]
      vtx_pass[kVertexAct]     = (nusel.max_track_length>3.0 || nusel.vertex_charge_per_pixel>50.0); // [9]      
      vtx_pass[kRecoFV]        = (reco_dwall>5.0); // [10]

      vtx_pass[kShowerLLCut]   = (nusel.largest_shower_avedqdx > 20.0 && nusel.largest_shower_avedqdx>20); // [11]

      vtx_pass[kWCPixel]       = (nusel.frac_allhits_on_cosmic<0.5); // [12]
      
      vtx_pass[kHadronic]      = (nusel.max_proton_pid<40 && (nusel.max_track_length>3.0 || nusel.vertex_hip_fraction>0.2)); // [13]

      vtx_pass[kAllCuts]       = true; // [14]
      
      // reco variable cuts only
      for ( int i=kMinShowerSize; i<kAllCuts; i++)
        vtx_pass[kAllCuts] = vtx_pass[kAllCuts] && vtx_pass[i];

      bool vtx_seq = true;
      for (int icut=0; icut<kNumCuts; icut++) {
        vtx_seq = vtx_seq && vtx_pass[icut]; // follows sequence
        // if still true mark as passing
        // or if previously passed, another vertex had passed this stage
        event_passes_cut[ icut ] = event_passes_cut[icut] || vtx_seq;
        if ( event_passes_cut[ icut ] && (nusel.dist2truevtx < cutpass_mindist[icut]) ) {
          cutpass_mindist[icut] = nusel.dist2truevtx;
          if ( nusel.truth_vtxFracNu>0.65 ) 
            cutpass_onnu[icut] = true;
        }
      }

      // for debug
      if (DEBUG) {
        std::cout << "[entry " << ientry << ", vtx" << ivtx << "] " << std::endl;
        vtx_seq = true;
        for (int i=0; i<=kAllCuts; i++) {
          vtx_seq = vtx_seq && vtx_pass[i]; // follows sequence        
          std::cout << "  " << selcut_names[i] << ": this=" << vtx_pass[i] << " chain=" << vtx_seq << std::endl;
        }
      }

      int min_connected_pass = 3;
      for (int p=0; p<(int)nusel.plane_connected_on_pass.size(); p++) {
        if ( nusel.plane_connected_on_pass[p]>0 && nusel.plane_connected_on_pass[p]<min_connected_pass )
          min_connected_pass = nusel.plane_connected_on_pass[p];
      }
      
    }//end of VERTEX LOOP
    std::sort( index_by_dist_v.begin(), index_by_dist_v.end() );

    // for debug
    if ( DEBUG) {
      std::cout << "[entry results] ------------" <<  std::endl;
      for (int i=0; i<=kAllCuts; i++) {
        std::cout << "  " << selcut_names[i] << ": " << event_passes_cut[i] << std::endl;
      }
      std::cout << "----------------------------" << std::endl;
      std::cout << "[ENTER] to continue" << std::endl;
      std::cin.get();
    }

    // FINAL STATE
    int event_final_state = -1;
    if ( is1l0p0pi==1 && evis_had>30.0 ) {
      // [0] 0 proton, 0 charged pion, 0 pi0
      if ( ccnc==0 && abs(nu_pdg)==12 )
        event_final_state = k1eVA;
      else if ( ccnc==0 && abs(nu_pdg)==14 )
        event_final_state = k1mVA;
    }
    else if ( is1l1p0pi==1 ) {
      // [1] 1 proton, 0 charged pion, 0 pi0
      if ( ccnc==0 && abs(nu_pdg)==12 )
        event_final_state = k1e_1p_0pi_nopi0;
      else if ( ccnc==0 && abs(nu_pdg)==14 )
        event_final_state = k1m_1p_0pi_nopi0;
      else if ( ccnc==1 )
        event_final_state = knc_1p_0pi_nopi0;
    }
    else if ( nproton_60mev==0 && nmeson_35mev==0 && npi0==0 ) {
      // [2] 0 proton, 1 charged pion, 0 pi0
      if ( ccnc==0 && abs(nu_pdg)==12 )
        event_final_state = k1e_0p_1pi_nopi0;
      else if ( ccnc==0 && abs(nu_pdg)==14 )
        event_final_state = k1m_0p_1pi_nopi0;
      else if ( ccnc==1 )
        event_final_state = knc_0p_1pi_nopi0;
    }
    else if ( nproton_60mev==0 && nmeson_35mev==0 && npi0>=1 ) {
      // [3] 0 proton, 0 charged pion, N pi0
      if ( ccnc==0 && abs(nu_pdg)==12 )
        event_final_state = k1e_0p_0pi_wpi0;
      else if ( ccnc==0 && abs(nu_pdg)==14 )
        event_final_state = k1m_0p_0pi_wpi0;
      else if ( ccnc==1 )
        event_final_state = knc_0p_0pi_wpi0;
    }
    else if ( nproton_60mev>=1 && nmeson_35mev==0 && npi0==0 ) {
      // [4] N proton, 0 charged pion, 0 neutral pion
      if ( ccnc==0 && abs(nu_pdg)==12 )
        event_final_state = k1e_Np_0pi_nopi0;
      else if ( ccnc==0 && abs(nu_pdg)==14 )
        event_final_state = k1m_Np_0pi_nopi0;
      else if ( ccnc==1 )
        event_final_state = knc_Np_0pi_nopi0;
    }
    else if ( nproton_60mev==1 && nmeson_35mev>=1 && npi0==0 ) {
      // [5] 1 proton, N charged pion, 0 neutral pion
      if ( ccnc==0 && abs(nu_pdg)==12 )
        event_final_state = k1e_1p_Npi_nopi0;
      else if ( ccnc==0 && abs(nu_pdg)==14 )
        event_final_state = k1m_1p_Npi_nopi0;
      else if ( ccnc==1 )
        event_final_state = knc_1p_Npi_nopi0;
    }
    else if ( nproton_60mev==1 && nmeson_35mev==0 && npi0>=1 ) {
      // [6] 1 proton, 0 charged pion, N neutral pion
      if ( ccnc==0 && abs(nu_pdg)==12 )
        event_final_state = k1e_1p_0pi_wpi0;
      else if ( ccnc==0 && abs(nu_pdg)==14 )
        event_final_state = k1m_1p_0pi_wpi0;
      else if ( ccnc==1 )
        event_final_state = knc_1p_0pi_wpi0;
    }
    else if ( nproton_60mev>=1 && nmeson_35mev>=1 && npi0>=1 ) {
      // [7] N proton, N charged pion, N neutral pion
      if ( ccnc==0 && abs(nu_pdg)==12 )
        event_final_state = k1e_Np_Npi_wpi0;
      else if ( ccnc==0 && abs(nu_pdg)==14 )
        event_final_state = k1m_Np_Npi_wpi0;
      else if ( ccnc==1 )
        event_final_state = knc_Np_Npi_wpi0;
    }

    // Enu filled based on if vertex passes cut stage
    bool still_passing = true;
    for (int icut=0; icut<=(int)kAllCuts; icut++) {
      still_passing = still_passing & event_passes_cut[icut];

      if ( !still_passing || event_final_state==-1 )
        break;

      if ( cutpass_onnu[icut] ) {
        henu[event_final_state][icut][kAllModes]->Fill( Enu_true );
        henu[event_final_state][icut][event_mode]->Fill( Enu_true );
      }
      else {
        henu[kOffNu][icut][kAllModes]->Fill( Enu_true );
        henu[kOffNu][icut][event_mode]->Fill( Enu_true );
      }
      
    }

    // std::cout << "[entry] to continue." << std::endl;
    // std::cin.get();
    
  }//end of entry loop

  out->Write();
  delete in;
  
  return 0;
}
