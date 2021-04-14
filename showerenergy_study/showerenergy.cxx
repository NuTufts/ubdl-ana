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

#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"

#include "larflow/Reco/NuSelectionVariables.h"
#include "larflow/Reco/NuVertexCandidate.h"
#include "larflow/Reco/cluster_functions.h"

#include "larlite/core/DataFormat/storage_manager.h"
#include "larlite/core/DataFormat/mcshower.h"

// define helper functions
std::vector<float> GetADCSum(larlite::larflowcluster shower,
                  std::vector<larcv::Image2D> wire_img,int threshold);
double GetShowerDir(larlite::pcaxis showerpc, larlite::event_mcshower* ev_mcshower);
double GetShowerStartDist(larlite::track showertrunk,larlite::event_mcshower* ev_mcshower);
int PickShower(std::vector<larlite::larflowcluster> shower_v);
void TestEventDisp(larlite::larflowcluster shower, TH2D** hshower_disp,
                  std::vector<larcv::Image2D> wire_img);
void ADCEventDisp(TH2D** h_disp,std::vector<larcv::Image2D> wire_img);


int main( int nargs, char** argv ) {
  //script to load in gen2 reco files, get shower calibration variables, and produce plots
  std::cout << "Gen-2 Shower Energy Calibration Plots" << std::endl;

  //start by loading all kpsreco files
  TFile* kps_f = new TFile(argv[1],"read");
  TTree* in = (TTree*)kps_f->Get("KPSRecoManagerTree");
  std::cout<<kps_f<<std::endl;

// load in reco file for wire tree
  std::string reco_f = argv[2];
  larcv::IOManager* io_larcv  = new larcv::IOManager(larcv::IOManager::kREAD,"IOManager_Tagger", larcv::IOManager::kTickBackward);
  io_larcv->reverse_all_products();
  io_larcv->add_in_file(reco_f);
  io_larcv->initialize();

  // load in larlite file for mc shower tree
  larlite::storage_manager* ioll = new larlite::storage_manager( larlite::storage_manager::kREAD );
  ioll->add_in_filename(reco_f);
  ioll->open();

  // load in true MC variables
  // need to add in true shower energies
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
  std::vector< larflow::reco::NuSelectionVariables >* pnu_sel_v = nullptr;
  in->SetBranchAddress( "nu_sel_v", &pnu_sel_v ); ///< selection variables per vertex

  // load in shower cluster object
  std::vector< larflow::reco::NuVertexCandidate>* nufitted_v = nullptr;
  in->SetBranchAddress( "nufitted_v",&nufitted_v);


  int nentries = in->GetEntries();

  // define outputfile
  std::string outdir = argv[3];
  std::string outfile = argv[4];
  // std::cout<<outdir+"/"+outfile+".root"<<std::endl;
  TFile* out = new TFile((outdir+"/"+outfile+".root").c_str(),"recreate");

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
         kAllCuts,        // All cuts applied except FV -- represents reco pass rate
         kNumCuts };       // Number in enum
  std::vector<std::string> selcut_names
    = { "fv",
        "vertexcand",
        "allreco",
        "numcuts"};

  // Cut variables for studying optimal cuts
  enum { kdwall=0,
         kdist2true,
         kNumCutVariables };

  std::vector<std::string> cutvar_names
    = { "dwall",
        "dist2true"};
  float cutvar_range[11][2] = { {-10,200},  // dwall
                                {0, 50 },   // distance to true vertex
  };
  int cutvar_nbins[11] = { 210, // dwall
                           150}; // dist 2 true

  // provides way to decide on remaining VA candidates per event
  // sort on shower-likelihood

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


  // Setup root histograms
  TH1D* henu[nsamples][kNumCuts]      = {0};
  TH1D* helep[nsamples][kNumCuts]     = {0};
  for (int isample=0; isample<nsamples; isample++) {
    for (int icut=0; icut<kNumCuts; icut++) {
        std::stringstream ss;
        ss << "hEnu_" << sample_names[isample] << "_" << selcut_names[icut] << "cut";
        henu[isample][icut] = new TH1D(ss.str().c_str(),";true E_{#nu} (MeV)", 30,0,3000);
        std::stringstream ss_lep;
        ss_lep << "hElep_" << sample_names[isample] << "_" << selcut_names[icut] << "cut";
        helep[isample][icut] = new TH1D(ss.str().c_str(),";true E_{#lep} (MeV)", 30,0,3000);
    }
  }

  // Distance to true vertex
  TH1D* hdist[nsamples] = {0};
  for (int i=0; i<nsamples; i++) {
    std::stringstream ss;
    ss << "hdist_" << sample_names[i];
    hdist[i] = new TH1D(ss.str().c_str(),";cm;",500,0,500);
  }

  // Distance to true vertex
  TH1D* cosangle[nsamples] = {0};
  for (int i=0; i<nsamples; i++) {
    std::stringstream ss;
    ss << "cosangle_" << sample_names[i];
    cosangle[i] = new TH1D(ss.str().c_str(),";cm;",90,0,180.0);
  }

  // number of shower pixels forward: proxy for "reco energy" until we build one
  TH1D* hnshower[nsamples] = {0};
  TH2D* hshower_vs_evislep[nsamples] = {0};
  TH2D* hshower_vs_enu[nsamples] = {0};
  for (int isample=0; isample<nsamples; isample++) {
    std::stringstream ss;
    ss << "hnshower_" << sample_names[isample];
    hnshower[isample] = new TH1D(ss.str().c_str(),";num shower hits;",100,0,10000);

    std::stringstream ss2d;
    ss2d << "hnshower_vs_enu_" << sample_names[isample];
    hshower_vs_enu[isample] = new TH2D(ss2d.str().c_str(), "",100,0,20000,30,0,3000);

    std::stringstream ss_vs_elep;
    ss_vs_elep << "hnshower_vs_evislep_" << sample_names[isample];
    hshower_vs_evislep[isample] = new TH2D(ss_vs_elep.str().c_str(), "",200,0,10000,150,0,3000);
  }

  // shower energy calibration plots

  TH2D* hshowerADC_vs_evislep[nsamples][3] = {0};
  TH1D* distance_between_reco_true_start[nsamples] = {0};
  for (int isample=0; isample<nsamples; isample++) {
    std::stringstream ss;
    ss << "distance_between_reco_true_start_" << sample_names[isample];
    distance_between_reco_true_start[isample] = new TH1D(ss.str().c_str(), "",1000,0,100);

    for (int p=0;p<3;p++){
      std::stringstream ss_vs_elep;
      ss_vs_elep << "hnshowerADC_vs_evislep_" << sample_names[isample]<<"_"<<p;
      hshowerADC_vs_evislep[isample][p] = new TH2D(ss_vs_elep.str().c_str(), "",50,0,400000,50,0,2000);
      hshowerADC_vs_evislep[isample][p]->SetOption("COLZ");
      hshowerADC_vs_evislep[isample][p]->SetXTitle("ADC Sum");
      hshowerADC_vs_evislep[isample][p]->SetYTitle("Evis_lepton (MeV)");

    }
  }
  // Event disp for testing
  TH2D* hshower_disp[3] = {0};
  TH2D* h_disp[3] = {0};
  for (int i=0; i<3; i++) {
    std::stringstream ss;
    ss << "hshower_disp_" << i;
    std::stringstream ssadc;
    ssadc << "h_disp_" << i;
    if (i==0){
      hshower_disp[i] = new TH2D(ss.str().c_str(),"",3456,0,3456.,1008,0,1008.);
      h_disp[i] = new TH2D(ssadc.str().c_str(),"",3456,0,3456.,1008,0,1008.);
    }
    if (i==1){
      hshower_disp[i] = new TH2D(ss.str().c_str(),"",3456,0,3456.,1008,0,1008.);
      h_disp[i] = new TH2D(ssadc.str().c_str(),"",3456,0,3456.,1008,0,1008.);
    }
    if (i==2){
      hshower_disp[i] = new TH2D(ss.str().c_str(),"",3456,0,3456.,1008,0,1008.);
      h_disp[i] = new TH2D(ssadc.str().c_str(),"",3456,0,3456.,1008,0,1008.);
    }

  }


  // loop through entries
  // std::cout<<"Loop through entries: "<<nentries<<std::endl;

  for (int ientry=0; ientry<nentries; ientry++) {
  // for (int ientry=0; ientry<1; ientry++) {

    if ( ientry%1==0 )
      std::cout << "[ ENTRY " << ientry << "]" << std::endl;

    // get the entry for each io manager
    in->GetEntry(ientry);
    io_larcv->read_entry(ientry);
    ioll->go_to(ientry);

    // load larcv/larlite branches
    const auto ev_img = (larcv::EventImage2D*)io_larcv->get_data( larcv::kProductImage2D, "wire" );
    auto const& wire_img = ev_img->Image2DArray();
    larlite::event_mcshower* ev_mcshower = (larlite::event_mcshower*)ioll->get_data(larlite::data::kMCShower,"mcreco");

    // truth cuts
    bool cut_fv = vtx_dwall>10.0;

    // find best reco vertex at each cut stage, measured by closeness to true vertex
    std::vector<EventDist2True_t> index_by_dist_v;
    std::vector<bool> event_passes_cut( kNumCuts, false ); //< indicates if event has a vertex that passes this cut stage
    int nvtx = (int)(*pnu_sel_v).size();
    // std::cout<<"nvtx: "<<nvtx<<std::endl;

    std::vector<std::vector<float>> ADCSum_vv;
    std::vector<double> dir_v;
    std::vector<double> dist_v;

    for (int ivtx=0; ivtx<nvtx; ivtx++) {

      auto const& nusel = (*pnu_sel_v)[ivtx];

      EventDist2True_t idx( nusel.dist2truevtx, ivtx );


      // get shower cluster!
      std::vector<larlite::larflowcluster> nufit_shower_v = (*nufitted_v)[ivtx].shower_v;
      std::vector<larlite::pcaxis> nufit_shower_pca_v = (*nufitted_v)[ivtx].shower_pcaxis_v;
      std::vector<larlite::track> nufit_shower_trunk_v = (*nufitted_v)[ivtx].shower_trunk_v;

      // for now, I only want nshowers == 1 for vertices in this study
      // if (nufit_shower_v.size()!=1) continue;
      // std::cout<<"vtx: "<<ivtx<<" number of showers: "<<nufit_shower_v.size()<<std::endl;
      // pick largest shower for now.
      int bestshower=PickShower(nufit_shower_v);

      // get sum of ADC in cluster
      std::vector<float> ADCSum_v;
      double dir;
      double dist;
      if (nufit_shower_v.size()>0){
        ADCSum_v = GetADCSum(nufit_shower_v[bestshower],wire_img,10);
        dir = GetShowerDir(nufit_shower_pca_v[bestshower],ev_mcshower);
        dist = GetShowerStartDist(nufit_shower_trunk_v[bestshower],ev_mcshower);
        // std::cout<<"ADCSums: "<<ADCSum_v[0]<<" "<<ADCSum_v[1]<<" "<<ADCSum_v[2]<<std::endl;
      }
      else {
        ADCSum_v={0.0,0.0,0.0};
      }
      ADCSum_vv.push_back(ADCSum_v);
      dir_v.push_back(dir);
      dist_v.push_back(dist);
      //make test event display for debugging
      // if (ivtx==0 && ientry==0) {
      //   TestEventDisp(nufit_shower_v[bestshower],hshower_disp,wire_img);
      //   ADCEventDisp(h_disp,wire_img);
      // }


      // selection cuts
      std::vector<bool> vtx_pass( kNumCuts, false );
      vtx_pass[kFV] = cut_fv;
      vtx_pass[kVertexCand3cm] = nusel.dist2truevtx<3.0;
      vtx_pass[kAllCuts]       = vtx_pass[kFV] && vtx_pass[kVertexCand3cm];

      bool vtx_seq = true;
      for (int icut=0; icut<kNumCuts; icut++) {
        vtx_seq = vtx_seq && vtx_pass[icut];
        event_passes_cut[ icut ] = event_passes_cut[icut] || vtx_seq;
      }

      if ( vtx_pass[kAllCuts] )
        index_by_dist_v.push_back( idx );

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

    }//end of vtx loop

    std::sort( index_by_dist_v.begin(), index_by_dist_v.end() );

    // Event-based plots

    // Enu filled based on if vertex passes cut stage
    bool still_passing = true;
    for (int icut=0; icut<=1; icut++) {
      still_passing = still_passing & event_passes_cut[icut];

      if ( !still_passing )
        break;

      // 1eVA
      if ( is1l0p0pi==1 && evis_had>30.0 ) {
        henu[k1eVA][icut]->Fill( Enu_true );
        helep[k1eVA][icut]->Fill( evis_lep );
      }

      // 1e1p
      if ( is1l1p0pi==1 ) {
        henu[k1e1p][icut]->Fill( Enu_true );
        helep[k1e1p][icut]->Fill( evis_lep );
      }

      // All
      henu[kAll][icut]->Fill( Enu_true );
      helep[kAll][icut]->Fill( evis_lep );
    }

    // nshowerhits: temporary proxy for neutrino energy
    if ( index_by_dist_v.size()>0 ) {
      // this means event had passing vertex

      int best_passing_vtx_index = index_by_dist_v.front().index;
      float num_shr_hits = (*pnu_sel_v)[best_passing_vtx_index].max_shower_nhits;
      // get number of reco showers
      std::vector<larlite::larflowcluster> nufit_shower_v = (*nufitted_v)[best_passing_vtx_index].shower_v;
      int nshower = nufit_shower_v.size();
      // std::cout<<"number of showers: "<<nshower<<std::endl;
      if (nshower>4) nshower=4;
      if(nshower==0) continue;
      // start filling histograms
      hnshower[kAll]->Fill( num_shr_hits );
      hshower_vs_enu[kAll]->Fill( num_shr_hits ,Enu_true);
      hshower_vs_evislep[kAll]->Fill(num_shr_hits, evis_lep );
      cosangle[kAll]->Fill(dir_v[best_passing_vtx_index]);
      // std::cout<<dir_v[best_passing_vtx_index]<<" degrees between reco and true"<<std::endl;
      distance_between_reco_true_start[kAll]->Fill(dist_v[best_passing_vtx_index]);
      if (ADCSum_vv[best_passing_vtx_index][2] < 200000 && ADCSum_vv[best_passing_vtx_index][2] > 150000 && evis_lep < 600){
        std::cout<<"event in lower distribution!"<<std::endl;
      }
      if (ADCSum_vv[best_passing_vtx_index][2] < 200000 && ADCSum_vv[best_passing_vtx_index][2] > 150000 && evis_lep > 600){
        std::cout<<"event in upper distribution!"<<std::endl;
      }
      if (abs(dir_v[best_passing_vtx_index]>0.95)){
        if (ADCSum_vv[best_passing_vtx_index][0] > 0) hshowerADC_vs_evislep[kAll][0]->Fill(ADCSum_vv[best_passing_vtx_index][0],evis_lep);
        if (ADCSum_vv[best_passing_vtx_index][1] > 0) hshowerADC_vs_evislep[kAll][1]->Fill(ADCSum_vv[best_passing_vtx_index][1],evis_lep);
        if (ADCSum_vv[best_passing_vtx_index][2] > 0) hshowerADC_vs_evislep[kAll][2]->Fill(ADCSum_vv[best_passing_vtx_index][2],evis_lep);
      }
      // std::cout<<"ADC SUMS: "<<ADCSum_vv[best_passing_vtx_index][0]<<" "<<ADCSum_vv[best_passing_vtx_index][1]<<" "<<ADCSum_vv[best_passing_vtx_index][2]<<std::endl;
      if ( is1l0p0pi==1 && evis_had>30.0 ) {
        hnshower[k1eVA]->Fill( num_shr_hits );
      	hshower_vs_enu[k1eVA]->Fill( num_shr_hits ,Enu_true );
      	hshower_vs_evislep[k1eVA]->Fill( num_shr_hits ,evis_lep);
        henu[k1eVA][kAllCuts]->Fill( Enu_true );
        helep[k1eVA][kAllCuts]->Fill( evis_lep );
        cosangle[k1eVA]->Fill(dir_v[best_passing_vtx_index]);
        distance_between_reco_true_start[k1eVA]->Fill(dist_v[best_passing_vtx_index]);
        if (abs(dir_v[best_passing_vtx_index]>0.95)){
          if (ADCSum_vv[best_passing_vtx_index][0] > 0) hshowerADC_vs_evislep[k1eVA][0]->Fill(ADCSum_vv[best_passing_vtx_index][0],evis_lep);
          if (ADCSum_vv[best_passing_vtx_index][1] > 0) hshowerADC_vs_evislep[k1eVA][1]->Fill(ADCSum_vv[best_passing_vtx_index][1],evis_lep);
          if (ADCSum_vv[best_passing_vtx_index][2] > 0) hshowerADC_vs_evislep[k1eVA][2]->Fill(ADCSum_vv[best_passing_vtx_index][2],evis_lep);
        }
      }

      if ( is1l1p0pi==1 ) {
        hnshower[k1e1p]->Fill( (*pnu_sel_v)[best_passing_vtx_index].max_shower_nhits );
	      hshower_vs_enu[k1e1p]->Fill( num_shr_hits ,Enu_true );
	      hshower_vs_evislep[k1e1p]->Fill(num_shr_hits , evis_lep );
        henu[k1e1p][kAllCuts]->Fill( Enu_true );
        helep[k1e1p][kAllCuts]->Fill( evis_lep );
        cosangle[k1e1p]->Fill(dir_v[best_passing_vtx_index]);
        distance_between_reco_true_start[k1e1p]->Fill(dist_v[best_passing_vtx_index]);
        if (abs(dir_v[best_passing_vtx_index]>0.95)){
          if (ADCSum_vv[best_passing_vtx_index][0] > 0) hshowerADC_vs_evislep[k1e1p][0]->Fill(ADCSum_vv[best_passing_vtx_index][0],evis_lep);
          if (ADCSum_vv[best_passing_vtx_index][1] > 0) hshowerADC_vs_evislep[k1e1p][1]->Fill(ADCSum_vv[best_passing_vtx_index][1],evis_lep);
          if (ADCSum_vv[best_passing_vtx_index][2] > 0) hshowerADC_vs_evislep[k1e1p][2]->Fill(ADCSum_vv[best_passing_vtx_index][2],evis_lep);
        }
      }


    }

    // std::cout << "[entry] to continue." << std::endl;
    // std::cin.get();

  }//end of entry loop

  // set some plot formating
  hshower_vs_evislep[k1e1p]->SetOption("COLZ");
  hshower_vs_evislep[k1e1p]->SetXTitle("Shower Hits");
  hshower_vs_evislep[k1e1p]->SetYTitle("Evis_lepton (MeV)");

  hshower_vs_evislep[k1eVA]->SetOption("COLZ");
  hshower_vs_evislep[k1eVA]->SetXTitle("Shower Hits");
  hshower_vs_evislep[k1eVA]->SetYTitle("Evis_lepton (MeV)");

  hshower_vs_evislep[kAll]->SetOption("COLZ");
  hshower_vs_evislep[kAll]->SetXTitle("Shower Hits");
  hshower_vs_evislep[kAll]->SetYTitle("Evis_lepton (MeV)");



  out->Write();
  io_larcv->finalize();
  delete in;

  return 0;
}

double GetShowerDir(larlite::pcaxis showerpc,larlite::event_mcshower* ev_mcshower){
  // get the cos of the angle between reco shower and true shower
  // initialize output
  double cosdir;

  // get reco direction shower
  std::vector<double> recoshowerdir;
  // get EigenVectors
  std::vector<std::vector<double> > showereigen_v = showerpc.getEigenVectors();
  recoshowerdir = showereigen_v[0];

  std::vector<double> trueshowerdir;
  TVector3 dir = ev_mcshower->at(0).StartDir();
  trueshowerdir ={dir.X(),dir.Y(),dir.Z()};

  // cos theta = (a dot b)/((mag a)* (mag (b)))
  double dot = recoshowerdir[0]*trueshowerdir[0]+recoshowerdir[1]*trueshowerdir[1]+recoshowerdir[2]*trueshowerdir[2];
  double magreco = sqrt(recoshowerdir[0]*recoshowerdir[0]+recoshowerdir[1]*recoshowerdir[1]+recoshowerdir[2]*recoshowerdir[2]);
  double magtrue = sqrt(trueshowerdir[0]*trueshowerdir[0]+trueshowerdir[1]*trueshowerdir[1]+trueshowerdir[2]*trueshowerdir[2]);
  cosdir = dot/(magreco*magtrue);
  double theta = acos (cosdir) * 180.0 / 3.14159;
  return theta;
}

double GetShowerStartDist(larlite::track showertrunk,larlite::event_mcshower* ev_mcshower){
  // get the distance between the true and reco shower start
  // initialize output
  double dist = 0;

  // get reco shower start
  std::vector<double> recoshowerstart;
  TVector3 reco_start_v = showertrunk.LocationAtPoint(0);
  recoshowerstart={reco_start_v.X(),reco_start_v.Y(),reco_start_v.Z()};

  // get reco shower end
  std::vector<double> recoshowerend;
  TVector3 reco_end_v = showertrunk.End();
  recoshowerend={reco_end_v.X(),reco_end_v.Y(),reco_end_v.Z()};

  // get true shower start
  std::vector<double> trueshowerstart;
  larlite::mcstep start = ev_mcshower->at(0).DetProfile();
  trueshowerstart ={start.X(),start.Y(),start.Z()};

  double xdist = recoshowerstart[0]-trueshowerstart[0];
  double ydist = recoshowerstart[1]-trueshowerstart[1];
  double zdist = recoshowerstart[2]-trueshowerstart[2];
  double dist_start = sqrt(xdist*xdist + ydist*ydist + zdist*zdist);

  xdist = recoshowerend[0]-trueshowerstart[0];
  ydist = recoshowerend[1]-trueshowerstart[1];
  zdist = recoshowerend[2]-trueshowerstart[2];
  double dist_end = sqrt(xdist*xdist + ydist*ydist + zdist*zdist);

  if (dist_start <= dist_end) dist = dist_start;
  else dist = dist_end;
  return dist;
}

std::vector<float> GetADCSum(larlite::larflowcluster shower,
              std::vector<larcv::Image2D> wire_img,int threshold ){
  //initialize output
  std::vector<float> sum_v;
  // / first turn to cluster
  //save all three planes
  std::vector<larlite::larflow3dhit>  shower_c = shower;
  // turn into 2d points (u,v,y,t)
  auto const& wireu_meta = wire_img.at(0).meta();
  auto const& wirev_meta = wire_img.at(1).meta();
  auto const& wirey_meta = wire_img.at(2).meta();
  // loop over planes
  for (int p =0;p<3;p++){
    larcv::ImageMeta meta =wire_img.at(p).meta();
    float sum =0;
    // make found list so we don't use the same point multiple times vector<(row,col)>
    // std::vector<std::vector<int>> foundpts;
    // loop over Hits
    float noutpix =0.0;
    for ( size_t ihit=0; ihit<shower_c.size(); ihit++ ) {
      int wire = shower_c[ihit].targetwire[p];
      int tick =  shower_c[ihit].tick;

      // add a check for projection failures
      if (tick < 0 || wire <0 || wire >=3456 ){
        noutpix +=1.0;
        continue;
      }

      int row = meta.row(tick);
      int col = meta.col(wire);
      int adcval =0;
      if (col < 3456 && row <1008){
        adcval = wire_img[p].pixel(row,col);
      }
      if (adcval >threshold){
        sum = sum+adcval;
      }
    }//end of loop over Hits
    if ((1.0-float(noutpix)/float(shower_c.size()))<=.98) sum_v.push_back(0.0);
    else sum_v.push_back(sum);

  }//end of loop over planes

  return sum_v;
}

int PickShower(std::vector<larlite::larflowcluster> shower_v){
  int bestshower = 0;
  int maxpix = 0;
  // loop through showers, pick largest for now
  for (int i =0; i<shower_v.size();i++){
    larflow::reco::cluster_t  shower_c = larflow::reco::cluster_from_larflowcluster( shower_v[i] );
    int numpix = shower_c.imgcoord_v.size();
    if (numpix>maxpix){
      bestshower= i;
      maxpix = numpix;
    }//end of if statement
  }//end of loop through showers
  return bestshower;
}

void TestEventDisp(larlite::larflowcluster shower, TH2D** hshower_disp,std::vector<larcv::Image2D> wire_img){
  std::cout<<"In Evt disp"<<std::endl;
  // we build a set of (tick,wire) pairs for look up reasons
  // get back to cluster_t objects


  larflow::reco::cluster_t  shower_c = larflow::reco::cluster_from_larflowcluster( shower );
  std::cout<<"number of pts "<<shower_c.imgcoord_v.size()<<std::endl;
  // turn into 2d points (u,v,y,t)
  auto const& wireu_meta = wire_img.at(0).meta();
  auto const& wirev_meta = wire_img.at(1).meta();
  auto const& wirey_meta = wire_img.at(2).meta();
  // loop over planes
  for (int p =0;p<3;p++){
    hshower_disp[p]->SetOption("COLZ");
    hshower_disp[p]->SetXTitle("Col");
    hshower_disp[p]->SetYTitle("Row");
    larcv::ImageMeta meta =wire_img.at(p).meta();

    // loop and set to one
    for (int c=1;c<3457;c++){
      for (int r=1;r<1009;r++){
          hshower_disp[p]->SetBinContent(c,r,1.0);
      }
    }

    // loop over Hits
    for ( size_t ihit=0; ihit<shower_c.imgcoord_v.size(); ihit++ ) {
      const std::vector<int>& coord = shower_c.imgcoord_v[ihit];
      int tick = coord[3];
      int wire = coord[p];
      int row = meta.row(tick);
      int col = meta.col(wire);
      int adcval = wire_img[p].pixel(row,col);
      hshower_disp[p]->SetBinContent(col,row,adcval);


    }//end of loop over Hits

  }//end of loop over planes

  return;
}

void ADCEventDisp(TH2D** h_disp,std::vector<larcv::Image2D> wire_img){
  std::cout<<"In ADC Evt disp"<<std::endl;
  // loop over planes
  for (int p =0;p<3;p++){
    h_disp[p]->SetOption("COLZ");
    h_disp[p]->SetXTitle("Col");
    h_disp[p]->SetYTitle("Row");

    // loop and set to one
    for (int c=0;c<3456;c++){
      for (int r=0;r<1008;r++){
          h_disp[p]->SetBinContent(c,r,1.0);
          if (wire_img[p].pixel(r,c)>0.0){
            h_disp[p]->SetBinContent(c,r,wire_img[p].pixel(r,c));
          }
          if (wire_img[p].pixel(r,c)>100.0){
            h_disp[p]->SetBinContent(c,r,100.0);
          }
      }
    }

  }//end of loop over planes
  return;
}
