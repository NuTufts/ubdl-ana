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
#include "larflow/Reco/NuSelCosmicTagger.h"
#include "larflow/Reco/TrackForwardBackwardLL.h"
#include "ublarcvapp/ubdllee/dwall.h"

int main( int nargs, char** argv ) {

  std::cout << "Gen-2 1e1p Selection Plots" << std::endl;

  std::ifstream input_list(argv[1]);
  char zinput[1028];  
  input_list >> zinput;

  int is_bnbnu = std::atoi(argv[2]);
  int is_mc = 0;
  if ( nargs>=4 )
    is_mc = std::atoi(argv[3]);

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
  std::vector< larflow::reco::NuVertexCandidate >* pnu_fitted_v = nullptr;  
  
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

  in->SetBranchAddress( "nu_sel_v",   &pnu_sel_v );    ///< selection variables per vertex
  in->SetBranchAddress( "nufitted_v", &pnu_fitted_v ); ///< reconstruction objects per vertex

  
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
         kUnrecoQ,        // [14] how much of the in-time pixels have been used
         kCosmicTag,      // [15] cosmic tagging cuts
         kBackMuon,       // [16] backward muon track check
         kAllCuts,        // [17] All cuts applied except FV -- represents reco pass rate
         kNumCuts };      // [18] Number in enum
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
        "unrecoq",        // [14]
        "cosmictag",      // [15]
        "backmuon",       // [16]
        "allreco",        // [17]
        "numcuts"};       // [18]

  // Cut variables for studying optimal cuts
  enum { kdwall=0, // [0]
         kdist2true,  // [1]
         kmaxshowerhits,  // [2]
         knshowerprongs,  // [3]
         kntrackprongs,   // [4]
         kllpid,          // [5]
         khipfraction,    // [6]
         kminshowergap,   // [7]
         kmaxshowergap,   // [8]
         kmintrackgap,    // [9]         
         kmaxtracklen,    // [10]         
         kvertexact,      // [12]
         klargestshowerll, // [13]
         kclosestshowerll, // [14]
         klargestshoweravedqdx, // [15]
         kclosestshoweravedqdx, // [16]
         knplanesconnected,     // [17]
         kminconnectpass,       // [18]
         ksecondshowersize,     // [19]
         kunrecoqmedfrac,       // [20]
         kbackwardmuonllr,      // [21]
         kNumCutVariables };    // [22]         
         
  std::vector<std::string> cutvar_names
    = { "dwall",                // [0]
        "dist2true",            // [1]
        "nmaxshowerhits",       // [2]
        "nshowerprongs",        // [3]
        "ntrackprongs",         // [4]
        "llpid",                // [5]
        "hipfraction",          // [6]
        "minshowergap",         // [7]
        "maxshowergap",         // [8]
        "mintrackgap",          // [9]
        "maxtracklen",          // [10]
        "vertexcharge",         // [11]
        "largestshowerll",      // [12]
        "closestshowerll",      // [13]
        "largestshoweravedqdx", // [14]
        "closestshoweravedqdx", // [15]
        "nplanesconnected",     // [16]
        "minconnectpass",       // [17]
        "secondshowersize",     // [18]
        "unrecoqmedfrac",       // [19]
        "backwardmuonllr"       // [20]
  };
  float cutvar_range[21][2] = { {-10,200},  // [0] dwall
                                {0, 50 },   // [1] distance to true vertex
                                {0, 10000}, // [2] hits in largest shower
                                {0, 10},    // [3] num shower prongs
                                {0, 10},    // [4] num track prongs
                                {-100,100}, // [5] proton likelihood
                                {0,1.01},   // [6] hip fraction
                                {0,50.0},   // [7] minshowergap
                                {0,50.0},   // [8] maxshowergap
                                {0,50.0},   // [9] mintrackgap
                                {0,500},    // [10] maxtracklen
                                {0,150.0},  // [11] vertex activity: charge per pixel around reco vertex
                                {-50,110},  // [12] largest shower likelihood
                                {-50,110},  // [13] largest shower ave dqdx
                                {0,200},    // [14] closest shower likelihood
                                {0,200},    // [15] closest shower ave dqdx
                                {0,4},      // [16] num connected planes
                                {0,4},      // [17] num connected planes
                                {0,10000},  // [18] second shower size
                                {0,1.0},    // [19] unreco charge, median fraction
                                {-100,100}  // [20] backward muon LLR
  };
  int cutvar_nbins[21] = { 210, // [0] dwall
                           150, // [1] dist 2 true
                           100, // [2] hits in largest shower
                           10,  // [3] num shower prongs
                           10,  // [4] num track prongs
                           200, // [5] proton likelihood
                           50,  // [6] hip fraction
                           50,  // [7] minshowergap
                           50,  // [8] maxshowergap
                           50,  // [9] mintrackgap
                           500, // [10] maxtracklen
                           50,  // [11] vertex activity
                           161, // [12] largest shower likelihood
                           161, // [13] largest shower ave dqdx
                           100, // [14] closest shower likelihood
                           100, // [15] closest shower ave dqdx
                           4,   // [16] nplanes connected
                           4,   // [17] min connected pass among planes
                           100, // [18] second shower size
                           20,  // [19] unreco charge, median fraction
                           51   // [20] backward muon LLR
  };


  //=========================================
  // Algorithms
  // -----------
  larflow::reco::NuSelCosmicTagger cosmictagger;
  cosmictagger.set_verbosity(larcv::msg::kDEBUG);

  larflow::reco::TrackForwardBackwardLL muvsproton;
  muvsproton.set_verbosity(larcv::msg::kINFO);
  
  // dq/dx plots: we will fill for vtx that passes vertex activity cuts
  // TH2D* hdqdx_shower_good = new TH2D("hdqdx_shower_good","",
  //                                    200, 0, 10.0,
  //                                    300, 0, 300.0 );
  // TH2D* hdqdx_shower_bad  = new TH2D("hdqdx_shower_bad","",
  //                                    200, 0, 10.0,
  //                                    300, 0, 300.0 );

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
  TH1D* henu[nsamples][kNumCuts][kNumModes]      = {0}; // for plotting Enu
  TH1D* henu_eff[nsamples][kNumCuts][kNumModes]  = {0}; // for plotting Eff
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
  
  // plot cut variables for
  //  (1) currently pass all cut,
  // w/ break down to see if on-nu or off-nu
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

    if ( ientry%100==0 )
      std::cout << "[ ENTRY " << ientry << "]" << std::endl;
    
    in->GetEntry(ientry);

    // truth cuts
    bool cut_fv = true;
    if ( is_mc==1 )
      cut_fv = vtx_dwall>10.0;

    // which mode are we
    int event_mode = 0;
    if ( is_mc==1 ) {
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

    }
    else {
      event_mode = kAllModes;
    }
      
    
    // find best reco vertex at each cut stage, measured by closeness to true vertex
    std::vector<EventDist2True_t> index_by_dist_v;
    std::vector<bool> event_passes_cut( kNumCuts, false ); //< indicates if event has a vertex that passes this cut stage
    int nvtx = (int)(*pnu_sel_v).size();
    for (int ivtx=0; ivtx<nvtx; ivtx++) {

      auto & nusel = (*pnu_sel_v)[ivtx];
      auto & nuvtx = (*pnu_fitted_v)[ivtx];

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

      // cosmic tagger
      cosmictagger.analyze( nuvtx, nusel );

      // backward muon versus forward proton
      muvsproton.analyze( nuvtx, nusel );

      // selection cuts
      std::vector<bool> vtx_pass( kNumCuts, false );
      vtx_pass[kFV] = cut_fv; // [0]

      if ( is_mc==1 ) {
        vtx_pass[kVertexCand3cm] = nusel.dist2truevtx<3.0; // [1] mc event
      }
      else {
        vtx_pass[kVertexCand3cm] = true; // [1] data event
      }
      
      vtx_pass[kMinShowerSize] = nusel.max_shower_nhits>500; // [2]
      vtx_pass[kNShowerProngs] = ( nusel.nshowers>0 && nusel.nshowers<=2 ); // [3]
      vtx_pass[kNTrackProngs]  = ( nusel.ntracks<=2 ); // [4]
      vtx_pass[kShowerGap]     = nusel.nplanes_connected>=2; // [5]
      vtx_pass[kTrackGap]      = (nusel.ntracks==0 || nusel.min_track_gap<3.0); // [6]
      vtx_pass[kMaxTrackLen]   = (nusel.ntracks>=1 && nusel.max_track_length<300.0); // [7]
      vtx_pass[kSecondShower]  = (nhits_second_shower<100); // [8]
      vtx_pass[kVertexAct]     = (nusel.max_track_length>3.0 || nusel.vertex_charge_per_pixel>50.0); // [9]      
      vtx_pass[kRecoFV]        = (reco_dwall>5.0); // [10]

      //vtx_pass[kShowerLLCut]   = (nusel.largest_shower_ll < 0.0 || nusel.closest_shower_ll < 0.0 ); // [11]
      vtx_pass[kShowerLLCut]   = (nusel.largest_shower_avedqdx > 20.0 && nusel.largest_shower_avedqdx>20 ); // [11]
      //vtx_pass[kShowerLLCut]   = true; // [11] pass for study      

      vtx_pass[kWCPixel]       = (nusel.frac_allhits_on_cosmic<0.5); // [12]
      
      //vtx_pass[kHadronic]      = (nusel.max_proton_pid<0 || nusel.vertex_hip_fraction>0.5); // [13]
      vtx_pass[kHadronic]      = (nusel.max_proton_pid<100 && nusel.vertex_hip_fraction>0.05) || (nusel.max_track_length<=3.0 && nusel.vertex_charge_per_pixel>50.0 ); // [13]
      //vtx_pass[kHadronic]      = true; // [13] pass for study

      // [14] Unreconstructed pixel fraction cut
      std::vector<float> unrecoq_v = nusel.unreco_fraction_v;
      float unrecoq_medfrac = 0.;
      std::sort( unrecoq_v.begin(), unrecoq_v.end() );
      if ( unrecoq_v.size()<3 )
        vtx_pass[kUnrecoQ]       = true;
      else {
        vtx_pass[kUnrecoQ]       = (unrecoq_v[1]<0.5);
        unrecoq_medfrac = unrecoq_v[1];
      }

      // [15] cosmic tagger
      vtx_pass[kCosmicTag] = true;      
      if ( nuvtx.track_v.size()>=2 ) {
        if ( cosmictagger._showercosmictag_maxbacktoback_dwall < 10.0
             && cosmictagger._showercosmictag_maxbacktoback_costrack < -0.8 ) {
          vtx_pass[kCosmicTag] = false;
        }
      }
      if ( nuvtx.track_v.size()>=1 ) {
        if ( cosmictagger._showercosmictag_maxboundarytrack_length>20.0
             && cosmictagger._showercosmictag_maxboundarytrack_verticalcos<-0.1
             && cosmictagger._showercosmictag_maxboundarytrack_showercos>0.0 ) {
          vtx_pass[kCosmicTag] = false;
        }
      }

      // [16] backward muon
      vtx_pass[kBackMuon] = true;
      float min_backwardmu_llr = 1000;
      for ( auto const& llr : muvsproton.best_llr_v ) {
        if ( llr<0 ) 
          vtx_pass[kBackMuon] = false;
        if ( llr < min_backwardmu_llr )
          min_backwardmu_llr = llr;
      }

      // set the all cuts flag
      // reco variable cuts only
      vtx_pass[kAllCuts]       = true;      
      for ( int i=kMinShowerSize; i<kAllCuts; i++)
        vtx_pass[kAllCuts] = vtx_pass[kAllCuts] && vtx_pass[i];

      bool vtx_seq = true;
      for (int icut=0; icut<kNumCuts; icut++) {
        vtx_seq = vtx_seq && vtx_pass[icut]; // follows sequence
        // if still true mark as passing
        // or if previously passed, another vertex had passed this stage
        event_passes_cut[ icut ] = event_passes_cut[icut] || vtx_seq;
      }

      if ( vtx_pass[kAllCuts] )
        index_by_dist_v.push_back( idx );

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

      // shower reco charge
      

      // Cut variables: study between "good" or "bad" vertex
      if ( nusel.dist2truevtx<2.0 || is_mc==0 ) {
        
        hvariable_good_v[kdwall]->Fill( vtx_dwall );
        if (is_mc==1)
          hvariable_good_v[kdist2true]->Fill( nusel.dist2truevtx );
        hvariable_good_v[kmaxshowerhits]->Fill( nusel.max_shower_nhits );
        hvariable_good_v[knshowerprongs]->Fill( nusel.nshowers );
        hvariable_good_v[kntrackprongs]->Fill( nusel.ntracks );
        hvariable_good_v[kllpid]->Fill( nusel.max_proton_pid );
        hvariable_good_v[khipfraction]->Fill( nusel.vertex_hip_fraction );
        if ( nusel.nshowers>0 ) hvariable_good_v[kminshowergap]->Fill( nusel.min_shower_gap );
        if ( nusel.nshowers>0 ) hvariable_good_v[kmaxshowergap]->Fill( nusel.max_shower_gap );
        if ( nusel.nshowers>0 ) hvariable_good_v[ksecondshowersize]->Fill( (float)nhits_second_shower );
        if ( nusel.ntracks>0 )  hvariable_good_v[kmintrackgap]->Fill( nusel.min_track_gap );        
        hvariable_good_v[kmaxtracklen]->Fill( nusel.max_track_length );
        hvariable_good_v[kvertexact]->Fill( nusel.vertex_charge_per_pixel );        
        hvariable_good_v[klargestshowerll]->Fill( nusel.largest_shower_ll );
        hvariable_good_v[klargestshoweravedqdx]->Fill( nusel.largest_shower_avedqdx );
        hvariable_good_v[kclosestshowerll]->Fill( nusel.closest_shower_ll );
        hvariable_good_v[kclosestshoweravedqdx]->Fill( nusel.closest_shower_avedqdx );
        hvariable_good_v[knplanesconnected]->Fill( nusel.nplanes_connected );
        hvariable_good_v[kminconnectpass]->Fill( min_connected_pass );
        hvariable_good_v[kunrecoqmedfrac]->Fill( unrecoq_medfrac );
        hvariable_good_v[kbackwardmuonllr]->Fill( min_backwardmu_llr );
      }
      else {
        
        hvariable_bad_v[kdwall]->Fill( vtx_dwall );
        hvariable_bad_v[kdist2true]->Fill( nusel.dist2truevtx );
        hvariable_bad_v[kmaxshowerhits]->Fill( nusel.max_shower_nhits );
        hvariable_bad_v[knshowerprongs]->Fill( nusel.nshowers );
        hvariable_bad_v[kntrackprongs]->Fill( nusel.ntracks );
        hvariable_bad_v[kllpid]->Fill( nusel.max_proton_pid );
        hvariable_bad_v[khipfraction]->Fill( nusel.vertex_hip_fraction );
        if ( nusel.nshowers>0 ) hvariable_bad_v[kminshowergap]->Fill( nusel.min_shower_gap );
        if ( nusel.nshowers>0 ) hvariable_bad_v[kmaxshowergap]->Fill( nusel.max_shower_gap );
        if ( nusel.nshowers>0 ) hvariable_bad_v[ksecondshowersize]->Fill( (float)nhits_second_shower );        
        if ( nusel.ntracks>0 )  hvariable_bad_v[kmintrackgap]->Fill( nusel.min_track_gap );          
        hvariable_bad_v[kmaxtracklen]->Fill( nusel.max_track_length );
        hvariable_bad_v[kvertexact]->Fill( nusel.vertex_charge_per_pixel );
        hvariable_bad_v[klargestshowerll]->Fill( nusel.largest_shower_ll );
        hvariable_bad_v[klargestshoweravedqdx]->Fill( nusel.largest_shower_avedqdx );
        hvariable_bad_v[kclosestshowerll]->Fill( nusel.closest_shower_ll );
        hvariable_bad_v[kclosestshoweravedqdx]->Fill( nusel.closest_shower_avedqdx );
        hvariable_bad_v[knplanesconnected]->Fill( nusel.nplanes_connected );
        hvariable_bad_v[kminconnectpass]->Fill( min_connected_pass );
        hvariable_bad_v[kunrecoqmedfrac]->Fill( unrecoq_medfrac );
        hvariable_bad_v[kbackwardmuonllr]->Fill( min_backwardmu_llr );        
      }

      // Cut variables: study correlation between reco state:
      // on vertex, off-vertex nu, off-vertex cosmic
      // determine state
      int vtx_reco_state = 0;
      if ( is_mc==1 ) {
        if ( nusel.dist2truevtx<3.0 ) {
          vtx_reco_state  = (int)kOnVertex;
        }
        else {
          if ( nusel.truth_vtxFracNu>0.65 )
            vtx_reco_state = (int)kOnNu;
          else
            vtx_reco_state = (int)kOffNu;
        }
      }

      if ( vtx_pass[kAllCuts] ) {
        hvar_onnu[vtx_reco_state][kdwall]->Fill( vtx_dwall );
        if ( is_mc==1 )
          hvar_onnu[vtx_reco_state][kdist2true]->Fill( nusel.dist2truevtx );
        hvar_onnu[vtx_reco_state][kmaxshowerhits]->Fill( nusel.max_shower_nhits );
        hvar_onnu[vtx_reco_state][knshowerprongs]->Fill( nusel.nshowers );
        hvar_onnu[vtx_reco_state][kntrackprongs]->Fill( nusel.ntracks );
        hvar_onnu[vtx_reco_state][kllpid]->Fill( nusel.max_proton_pid );
        hvar_onnu[vtx_reco_state][khipfraction]->Fill( nusel.vertex_hip_fraction );
        if ( nusel.nshowers>0 ) hvar_onnu[vtx_reco_state][kminshowergap]->Fill( nusel.min_shower_gap );
        if ( nusel.nshowers>0 ) hvar_onnu[vtx_reco_state][kmaxshowergap]->Fill( nusel.max_shower_gap );
        if ( nusel.nshowers>0 ) hvar_onnu[vtx_reco_state][ksecondshowersize]->Fill( (float)nhits_second_shower );
        if ( nusel.ntracks>0 )  hvar_onnu[vtx_reco_state][kmintrackgap]->Fill( nusel.min_track_gap );
        hvar_onnu[vtx_reco_state][kmaxtracklen]->Fill( nusel.max_track_length );
        hvar_onnu[vtx_reco_state][kvertexact]->Fill( nusel.vertex_charge_per_pixel );
        hvar_onnu[vtx_reco_state][klargestshowerll]->Fill( nusel.largest_shower_ll );
        hvar_onnu[vtx_reco_state][klargestshoweravedqdx]->Fill( nusel.largest_shower_avedqdx );
        hvar_onnu[vtx_reco_state][kclosestshowerll]->Fill( nusel.closest_shower_ll );
        hvar_onnu[vtx_reco_state][kclosestshoweravedqdx]->Fill( nusel.closest_shower_avedqdx );
        hvar_onnu[vtx_reco_state][knplanesconnected]->Fill( nusel.nplanes_connected );
        hvar_onnu[vtx_reco_state][kminconnectpass]->Fill( min_connected_pass );
        hvar_onnu[vtx_reco_state][kunrecoqmedfrac]->Fill( unrecoq_medfrac );
        hvar_onnu[vtx_reco_state][kbackwardmuonllr]->Fill( min_backwardmu_llr );        
      }

      // dQ/dx plots: need to save dqdx data, which didnt ..
      // if ( vtx_dwall>10.0 && vtx_pass[kVertexAct] ) {

      //   // we have to recover the dqdx track (cannot, no image2d ...)
      //   NuSelShowerTrunkAna
      // if ( nusel.dist2truevtx<2.0 && 
      //   int largest_shower_idx = 0;
      //   int max_nhits = 0;
      //   for (int ishower=0; ishower<(int)nuvtx.shower_v.size(); ishower++) {
      //     if ( max_nhits<(int)nuvtx.shower_v.size() ) {
      //       max_nhits = (int)nuvtx.shower_v.size();
      //       largest_shower_idx = ishower;
      //     }
      //   }

        
        
      // }
      // else if ( nusel.dist2truevtx>2.0 && vtx_dwall>10.0 && vtx_pass[kVertexAct] ) {
      // }
      
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
    
    // Event-based plots

    // Eff vs. Enu true filled based on if vertex passes cut stage
    bool still_passing = true;
    for (int icut=0; icut<=(int)kAllCuts; icut++) {
      still_passing = still_passing & event_passes_cut[icut];

      if ( !still_passing )
        break;

      if ( is_mc==1 ) {
        // 1eVA
        if ( is1l0p0pi==1 && evis_had>30.0 ) {
          //henu[k1eVA][icut][kAllModes]->Fill( Enu_true );
          henu_eff[k1eVA][icut][kAllModes]->Fill( Enu_true );
          //henu[k1eVA][icut][event_mode]->Fill( Enu_true );
          henu_eff[k1eVA][icut][event_mode]->Fill( Enu_true );        
        }
        
        // 1e1p
        if ( is1l1p0pi==1 ) {
          //henu[k1e1p][icut][kAllModes]->Fill( Enu_true );
          henu_eff[k1e1p][icut][kAllModes]->Fill( Enu_true );
          //henu[k1e1p][icut][event_mode]->Fill( Enu_true );
          henu_eff[k1e1p][icut][event_mode]->Fill( Enu_true );                
        }
      }

      // All
      //henu[kAll][icut][kAllModes]->Fill( Enu_true );
      henu_eff[kAll][icut][kAllModes]->Fill( Enu_true );
      //henu[kAll][icut][event_mode]->Fill( Enu_true );
      henu_eff[kAll][icut][event_mode]->Fill( Enu_true );                      
    }//end of loop over cut sequence

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

      if ( is_mc==1 ) {

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

      
    }//end of if passing vertex
    
    // std::cout << "[entry] to continue." << std::endl;
    // std::cin.get();
    
  }//end of entry loop

  // efficiency: divide at time of plotting
  // for (int i=0; i<nsamples; i++) {
  //   for (int j=0; j<kAllCuts; j++) {
  //     henu_eff[i][j][kAllModes]->Divide( henu[i][kFV][kAllModes] );
  //   }
  // }
  
  out->Write();
  delete in;
  
  return 0;
}
