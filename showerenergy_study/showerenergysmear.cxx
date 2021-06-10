#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
// root
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
// larcv
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventPGraph.h"
// larflow
#include "larflow/Reco/NuSelectionVariables.h"
#include "larflow/Reco/NuVertexCandidate.h"
#include "larflow/Reco/cluster_functions.h"
// larlite
#include "larlite/core/DataFormat/storage_manager.h"
#include "larlite/core/DataFormat/mcshower.h"
#include "larlite/core/DataFormat/mctruth.h"
// larutil
#include "LArUtil/LArProperties.h"
#include "LArUtil/DetectorProperties.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/ClockConstants.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"

// define helper functions
std::vector<float> GetADCSum(larlite::larflowcluster shower,
                  std::vector<larcv::Image2D> wire_img,int threshold);
std::vector<float> GetADCSumWithSmear(larlite::larflowcluster shower,
                  std::vector<larcv::Image2D> wire_img,int threshold );
double GetShowerDir(larlite::pcaxis showerpc, larlite::event_mcshower* ev_mcshower);
double GetShowerStartDist(larlite::track showertrunk,larlite::event_mcshower* ev_mcshower);
int PickShower(std::vector<larlite::larflowcluster> shower_v);
void TestEventDisp(std::vector<larlite::larflowcluster> shower_v, std::vector<float> nufit_vertex_v,
                  std::vector<larcv::Image2D> wire_img,std::vector<larcv::Image2D> thru_img,
                  int run, int subrun, int event, int low);
void ADCEventDisp(TH2D** h_disp,std::vector<larcv::Image2D> wire_img);
std::vector<int> getProjectedPixel( const std::vector<double>& pos3d,
            const larcv::ImageMeta& meta,
            const int nplanes,
            const float fracpixborder=1.5  );
void printinfo(larlite::event_mctruth* ev_mctruth);


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
  int npi0;
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
  in->SetBranchAddress( "npi0", &npi0 );
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

  // truth categories - 4 types
  const int nsamples = 4;
  enum { k1eVA=0, k1e1p, kAll, kpi0 };
  std::vector<std::string> sample_names = {"is1eVA","1e1p","all","pi0"};
  std::vector<std::string> sample_cuts  = {"is1l0p0pi==1 && evis_had>30.0",
                                           "is1l1p0pi==1",
                                           "1==1",
                                           "npi0>0"};

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

  TH2D* hshowerADC_vs_evislep[nsamples][4] = {0};
  TH1D* distance_between_reco_true_start[nsamples] = {0};
  for (int isample=0; isample<nsamples; isample++) {
    std::stringstream ss;
    ss << "distance_between_reco_true_start_" << sample_names[isample];
    distance_between_reco_true_start[isample] = new TH1D(ss.str().c_str(), "",1000,0,100);

    for (int p=0;p<4;p++){
      std::stringstream ss_vs_elep;
      ss_vs_elep << "hnshowerADC_vs_evislep_" << sample_names[isample]<<"_"<<p;
      hshowerADC_vs_evislep[isample][p] = new TH2D(ss_vs_elep.str().c_str(), "",1000,0,250000,1000,0,2000);
      hshowerADC_vs_evislep[isample][p]->SetOption("COLZ");
      hshowerADC_vs_evislep[isample][p]->SetXTitle("ADC Sum");
      hshowerADC_vs_evislep[isample][p]->SetYTitle("Evis_lepton (MeV)");

    }
  }

  // vertex position Plots
  // TH1D* hvtx_pos_low_r[3] = {0};
  // TH1D* hvtx_pos_high_r[3] = {0};
  // TH1D* hvtx_pos_all_r[3] = {0};
  // TH1D* hvtx_pos_low_c[3] = {0};
  // TH1D* hvtx_pos_high_c[3] = {0};
  // TH1D* hvtx_pos_all_c[3] = {0};
  // for (int i=0; i<3; i++) {
  //   std::stringstream ssvtxall;
  //   ssvtxall << "hvtx_pos_all" << i;
  //   std::stringstream ssvtxlow;
  //   ssvtxlow << "hvtx_pos_low" << i;
  //   std::stringstream ssvtxhigh;
  //   ssvtxhigh << "hvtx_pos_high" << i;
  //   hvtx_pos_low[i] = new TH2D(ssvtxlow.str().c_str(),"",3456,0,3456.,1008,0,1008.);
  //   hvtx_pos_high[i] = new TH2D(ssvtxhigh.str().c_str(),"",3456,0,3456.,1008,0,1008.);
  //   hvtx_pos_all[i] = new TH2D(ssvtxall.str().c_str(),"",3456,0,3456.,1008,0,1008.);
  //   for (int c=1;c<3457;c++){
  //     for (int r=1;r<1009;r++){
  //         hvtx_pos_low[i]->SetBinContent(c,r,0.000001);
  //         hvtx_pos_high[i]->SetBinContent(c,r,0.0000001);
  //         hvtx_pos_all[i]->SetBinContent(c,r,0.00000001);
  //     }
  //   }
  //
  // }


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

    larlite::event_mctruth* ev_mctruth = (larlite::event_mctruth*)ioll->get_data(larlite::data::kMCTruth,  "generator" );
    if(true){
      printinfo(ev_mctruth);
    }

    // load larcv/larlite branches
    const auto ev_img = (larcv::EventImage2D*)io_larcv->get_data( larcv::kProductImage2D, "wire" );
    const auto ev_thru = (larcv::EventImage2D*)io_larcv->get_data( larcv::kProductImage2D, "thrumu" );
    int run = ev_img->run();
		int subrun = ev_img->subrun();
		int event = ev_img->event();
    auto const& wire_img = ev_img->Image2DArray();
    auto const& thru_img = ev_thru->Image2DArray();
    larlite::event_mcshower* ev_mcshower = (larlite::event_mcshower*)ioll->get_data(larlite::data::kMCShower,"mcreco");
    larcv::EventPGraph* ev_test_pgraph  = (larcv::EventPGraph*)(io_larcv->get_data(larcv::kProductPGraph,"test"));

    // truth cuts
    bool cut_fv = vtx_dwall>10.0;

    // find best reco vertex at each cut stage, measured by closeness to true vertex
    std::vector<EventDist2True_t> index_by_dist_v;
    std::vector<bool> event_passes_cut( kNumCuts, false ); //< indicates if event has a vertex that passes this cut stage
    int nvtx = (int)(*pnu_sel_v).size();
    // std::cout<<"nvtx: "<<nvtx<<std::endl;

    std::vector<std::vector<float>> ADCSum_vv;
    std::vector<float> adcaverage_v;
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
      std::vector<float> ADCSumSmear_v;
      double dir;
      double dist;
      float adcaverage;
      if (nufit_shower_v.size()>0){
        // ADCSum_v = GetADCSum(nufit_shower_v[bestshower],wire_img,10);
        ADCSum_v = GetADCSumWithSmear(nufit_shower_v[bestshower],wire_img,10);
        // std::cout<<"Old: "<<ADCSum_v[2]<<" New: "<<ADCSumSmear_v[2]<<std::endl;
        dir = GetShowerDir(nufit_shower_pca_v[bestshower],ev_mcshower);
        dist = GetShowerStartDist(nufit_shower_trunk_v[bestshower],ev_mcshower);
        adcaverage = (ADCSum_v[0]+ADCSum_v[1]+ADCSum_v[2])/(3.0);
        // std::cout<<"ADCSums: "<<ADCSum_v[0]<<" "<<ADCSum_v[1]<<" "<<ADCSum_v[2]<<std::endl;
      }
      else {
        ADCSum_v={0.0,0.0,0.0};
        adcaverage=0.0;
      }
      ADCSum_vv.push_back(ADCSum_v);
      adcaverage_v.push_back(adcaverage);
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
        std::cout<<"is 1eVA"<<std::endl;
      }

      // 1e1p
      if ( is1l1p0pi==1 ) {
        henu[k1e1p][icut]->Fill( Enu_true );
        helep[k1e1p][icut]->Fill( evis_lep );
        std::cout<<"is 1e1p"<<std::endl;
      }

      // pi0
      if ( npi0>0 ) {
        std::cout<<"HAS PI0: "<<npi0<<std::endl;
        henu[kpi0][icut]->Fill( Enu_true );
        helep[kpi0][icut]->Fill( evis_lep );
        std::cout<<"done"<<std::endl;
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
      // get vertex
      std::vector<float> nufit_vertex_v = (*nufitted_v)[best_passing_vtx_index].pos;
      std::cout<< nufit_vertex_v[0]<< nufit_vertex_v[1]<< nufit_vertex_v[2]<<std::endl;
      if (false){
        if (ADCSum_vv[best_passing_vtx_index][2] < 200000 && ADCSum_vv[best_passing_vtx_index][2] > 150000 && evis_lep < 600  ){
          std::cout<<"event in lower distribution!"<<std::endl;
          TestEventDisp(nufit_shower_v, nufit_vertex_v ,wire_img,thru_img,run,subrun,event,1);
          std::vector<double> reco_vertex = {(double)nufit_vertex_v[0],(double)nufit_vertex_v[1], (double)nufit_vertex_v[2]};
          // now project into 2d
          std::vector<int> vtx_rc = getProjectedPixel(reco_vertex, wire_img[0].meta(), 3);
          // fill root histograms
          // float binvalu = hvtx_pos_low[0]->GetBinContent(vtx_rc[1], vtx_rc[0]);
          // float binvalv = hvtx_pos_low[1]->GetBinContent(vtx_rc[2], vtx_rc[0]);
          // float binvaly = hvtx_pos_low[2]->GetBinContent(vtx_rc[3], vtx_rc[0]);
          // hvtx_pos_low[0]->SetBinContent(vtx_rc[1], vtx_rc[0],binvalu+1.0);
          // hvtx_pos_low[1]->SetBinContent(vtx_rc[2], vtx_rc[0],binvalv+1.0);
          // hvtx_pos_low[2]->SetBinContent(vtx_rc[3], vtx_rc[0],binvaly+1.0);

        }
        else if (ADCSum_vv[best_passing_vtx_index][2] < 200000 && ADCSum_vv[best_passing_vtx_index][2] > 150000 && evis_lep > 600 ){
          std::cout<<"event in upper distribution!"<<std::endl;
          TestEventDisp(nufit_shower_v, nufit_vertex_v, wire_img,thru_img,run,subrun,event,0);
          std::vector<double> reco_vertex = {(double)nufit_vertex_v[0],(double)nufit_vertex_v[1], (double)nufit_vertex_v[2]};
          // now project into 2d
          std::vector<int> vtx_rc = getProjectedPixel(reco_vertex, wire_img[0].meta(), 3);
          // fill root histograms
          // float binvalu = hvtx_pos_high[0]->GetBinContent(vtx_rc[1], vtx_rc[0]);
          // float binvalv = hvtx_pos_high[1]->GetBinContent(vtx_rc[2], vtx_rc[0]);
          // float binvaly = hvtx_pos_high[2]->GetBinContent(vtx_rc[3], vtx_rc[0]);
          // hvtx_pos_high[0]->SetBinContent(vtx_rc[1], vtx_rc[0],binvalu+1.0);
          // hvtx_pos_high[1]->SetBinContent(vtx_rc[2], vtx_rc[0],binvalv+1.0);
          // hvtx_pos_high[2]->SetBinContent(vtx_rc[3], vtx_rc[0],binvaly+1.0);

        }
        else{
          std::vector<double> reco_vertex = {(double)nufit_vertex_v[0],(double)nufit_vertex_v[1], (double)nufit_vertex_v[2]};
          // now project into 2d
          std::vector<int> vtx_rc = getProjectedPixel(reco_vertex, wire_img[0].meta(), 3);
          // // fill root histograms
          // float binvalu = hvtx_pos_all[0]->GetBinContent(vtx_rc[1], vtx_rc[0]);
          // float binvalv = hvtx_pos_all[1]->GetBinContent(vtx_rc[2], vtx_rc[0]);
          // float binvaly = hvtx_pos_all[2]->GetBinContent(vtx_rc[3], vtx_rc[0]);
          // hvtx_pos_all[0]->SetBinContent(vtx_rc[1], vtx_rc[0],binvalu+1.0);
          // hvtx_pos_all[1]->SetBinContent(vtx_rc[2], vtx_rc[0],binvalv+1.0);
          // hvtx_pos_all[2]->SetBinContent(vtx_rc[3], vtx_rc[0],binvaly+1.0);

        }
      }
      if (abs(dir_v[best_passing_vtx_index]>0.95)){
        if (ADCSum_vv[best_passing_vtx_index][0] > 0) hshowerADC_vs_evislep[kAll][0]->Fill(ADCSum_vv[best_passing_vtx_index][0],evis_lep);
        if (ADCSum_vv[best_passing_vtx_index][1] > 0) hshowerADC_vs_evislep[kAll][1]->Fill(ADCSum_vv[best_passing_vtx_index][1],evis_lep);
        if (ADCSum_vv[best_passing_vtx_index][2] > 0) hshowerADC_vs_evislep[kAll][2]->Fill(ADCSum_vv[best_passing_vtx_index][2],evis_lep);
        if (adcaverage_v[best_passing_vtx_index] > 0) hshowerADC_vs_evislep[kAll][3]->Fill(adcaverage_v[best_passing_vtx_index],evis_lep);
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
          if (adcaverage_v[best_passing_vtx_index] > 0) hshowerADC_vs_evislep[k1eVA][3]->Fill(adcaverage_v[best_passing_vtx_index],evis_lep);
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
          if (adcaverage_v[best_passing_vtx_index] > 0) hshowerADC_vs_evislep[k1e1p][3]->Fill(adcaverage_v[best_passing_vtx_index],evis_lep);
        }
      }

      if ( npi0>0 ) {
        hnshower[kpi0]->Fill( (*pnu_sel_v)[best_passing_vtx_index].max_shower_nhits );
	      hshower_vs_enu[kpi0]->Fill( num_shr_hits ,Enu_true );
	      hshower_vs_evislep[kpi0]->Fill(num_shr_hits , evis_lep );
        henu[kpi0][kAllCuts]->Fill( Enu_true );
        helep[kpi0][kAllCuts]->Fill( evis_lep );
        cosangle[kpi0]->Fill(dir_v[best_passing_vtx_index]);
        distance_between_reco_true_start[kpi0]->Fill(dist_v[best_passing_vtx_index]);
        if (abs(dir_v[best_passing_vtx_index]>0.95)){
          if (ADCSum_vv[best_passing_vtx_index][0] > 0) hshowerADC_vs_evislep[kpi0][0]->Fill(ADCSum_vv[best_passing_vtx_index][0],evis_lep);
          if (ADCSum_vv[best_passing_vtx_index][1] > 0) hshowerADC_vs_evislep[kpi0][1]->Fill(ADCSum_vv[best_passing_vtx_index][1],evis_lep);
          if (ADCSum_vv[best_passing_vtx_index][2] > 0) hshowerADC_vs_evislep[kpi0][2]->Fill(ADCSum_vv[best_passing_vtx_index][2],evis_lep);
          if (adcaverage_v[best_passing_vtx_index] > 0) hshowerADC_vs_evislep[kpi0][3]->Fill(adcaverage_v[best_passing_vtx_index],evis_lep);
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

  hshower_vs_evislep[kpi0]->SetOption("COLZ");
  hshower_vs_evislep[kpi0]->SetXTitle("Shower Hits");
  hshower_vs_evislep[kpi0]->SetYTitle("Evis_lepton (MeV)");

  hshower_vs_evislep[kAll]->SetOption("COLZ");
  hshower_vs_evislep[kAll]->SetXTitle("Shower Hits");
  hshower_vs_evislep[kAll]->SetYTitle("Evis_lepton (MeV)");

  // for(int i=0;i<3;i++){
  //   hvtx_pos_all[i]->SetOption("COLZ");
  //   hvtx_pos_all[i]->RebinX(8);
  //   hvtx_pos_all[i]->RebinY(8);
  //   hvtx_pos_low[i]->SetOption("COLZ");
  //   hvtx_pos_low[i]->RebinX(8);
  //   hvtx_pos_low[i]->RebinY(8);
  //   hvtx_pos_high[i]->SetOption("COLZ");
  //   hvtx_pos_high[i]->RebinX(8);
  //   hvtx_pos_high[i]->RebinY(8);
  // }



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
    int npixels =0;
    for (int c=1;c<3456;c++){
      for (int r=1;r<1008;r++){
          if (wire_img[p].pixel(r,c) <threshold) wire_img[p].set_pixel(r,c,0);
      }
    }
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
      float adcval =0;
      if (col < 3456 && row <1008){
        adcval = wire_img[p].pixel(row,col);
        npixels+=1;
      }
      if (true){
        sum = sum+adcval;
      }
    }//end of loop over Hits
    if ((1.0-float(noutpix)/float(shower_c.size()))<=.98) sum_v.push_back(0.0);
    else sum_v.push_back(sum);
    std::cout<<"old: "<<p<<" "<<npixels<<std::endl;

  }//end of loop over planes

  return sum_v;
}

std::vector<float> GetADCSumWithSmear(larlite::larflowcluster shower,
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
    // initialize alredy used points so we don't double count
    larcv::Image2D usedpixels = wire_img[p];

    // loop and set to zero
    for (int c=1;c<3456;c++){
      for (int r=1;r<1008;r++){
          usedpixels.set_pixel(r,c,0);
          if (wire_img[p].pixel(r,c) <threshold) wire_img[p].set_pixel(r,c,0);
      }
    }
    int npixels =0;
    int numused =0;
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
      float adcval =0;
      // std::cout<<row<<" "<<col<<std::endl;

      // we want to get all surrounding pixels
      // updating boundary conditions - adding extra 1 pixel cut off
      if (col < 3455 && row <1007 && col > 1 && row > 1){
        // for each pixel, check if in mask
        if (usedpixels.pixel(row,col)==0){
        // if (true){
            adcval += wire_img[p].pixel(row,col);
            npixels+=1;
            usedpixels.set_pixel(row,col,1);
        }
        // else numused+=1;
        // now get surrounding pixels
        // [(r+1,c-1),(r+1,c),(r+1,c+1),(r,c-1),(r,c+1),(r-1,c-1),(r-1,c),(r-1,c+1)]
        if (usedpixels.pixel(row+1,col-1)==0){
            adcval += wire_img[p].pixel(row+1,col-1);
            usedpixels.set_pixel(row+1,col-1,1);
        }
        if (usedpixels.pixel(row+1,col)==0){
            adcval += wire_img[p].pixel(row+1,col);
            usedpixels.set_pixel(row+1,col,1);
        }
        if (usedpixels.pixel(row+1,col+1)==0){
            adcval += wire_img[p].pixel(row+1,col+1);
            usedpixels.set_pixel(row+1,col+1,1);
        }
        if (usedpixels.pixel(row,col-1)==0){
            adcval += wire_img[p].pixel(row,col-1);
            usedpixels.set_pixel(row,col-1,1);
        }
        if (usedpixels.pixel(row,col+1)==0){
            adcval += wire_img[p].pixel(row,col+1);
            usedpixels.set_pixel(row,col+1,1);
        }
        if (usedpixels.pixel(row-1,col-1)==0){
            adcval += wire_img[p].pixel(row-1,col-1);
            usedpixels.set_pixel(row-1,col-1,1);
        }
        if (usedpixels.pixel(row-1,col)==0){
            adcval += wire_img[p].pixel(row-1,col);
            usedpixels.set_pixel(row-1,col,1);
        }
        if (usedpixels.pixel(row-1,col+1)==0){
            adcval += wire_img[p].pixel(row-1,col+1);
            usedpixels.set_pixel(row-1,col+1,1);
        }
      }
      sum = sum+adcval;

    }//end of loop over Hits
    std::cout<<"new: "<<p<<" "<<npixels<<"+ "<<numused<<std::endl;

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

void TestEventDisp(std::vector<larlite::larflowcluster> shower_v,std::vector<float> nufit_vertex_v ,
      std::vector<larcv::Image2D> wire_img,std::vector<larcv::Image2D> thru_img,
      int run, int subrun, int event,int low){

  std::cout<<"In Evt disp"<<std::endl;
  // we build a set of (tick,wire) pairs for look up reasons
  // get back to cluster_t objects

  gStyle->SetOptStat(0);
  // construct histograms and save - if low ==0 ,in high dist, else in low
  TH2D* hshower_disp[3] = {0};
  TH2D* h_disp[3] = {0};
  TH2D* hvtx_disp[3] = {0};
  for (int i=0; i<3; i++) {
    std::stringstream ss;
    ss << "hshower_disp_" << i;
    std::stringstream ssadc;
    ssadc << "h_disp_" << i;
    std::stringstream ssvtx;
    ssvtx << "hvtx_disp_" << i;
    if (i==0){
      hshower_disp[i] = new TH2D(ss.str().c_str(),"",3456,0,3456.,1008,0,1008.);
      hvtx_disp[i] =new TH2D(ssvtx.str().c_str(),"",3456,0,3456.,1008,0,1008.);
      h_disp[i] = new TH2D(ssadc.str().c_str(),"",3456,0,3456.,1008,0,1008.);
    }
    if (i==1){
      hshower_disp[i] = new TH2D(ss.str().c_str(),"",3456,0,3456.,1008,0,1008.);
      h_disp[i] = new TH2D(ssadc.str().c_str(),"",3456,0,3456.,1008,0,1008.);
      hvtx_disp[i] =new TH2D(ssvtx.str().c_str(),"",3456,0,3456.,1008,0,1008.);
    }
    if (i==2){
      hshower_disp[i] = new TH2D(ss.str().c_str(),"",3456,0,3456.,1008,0,1008.);
      h_disp[i] = new TH2D(ssadc.str().c_str(),"",3456,0,3456.,1008,0,1008.);
      hvtx_disp[i] =new TH2D(ssvtx.str().c_str(),"",3456,0,3456.,1008,0,1008.);
    }

  }

  // turn into 2d points (u,v,y,t)
  auto const& wireu_meta = wire_img.at(0).meta();
  auto const& wirev_meta = wire_img.at(1).meta();
  auto const& wirey_meta = wire_img.at(2).meta();
  // loop over planes
  for (int p =0;p<3;p++){
    hshower_disp[p]->SetOption("COLZ");
    hshower_disp[p]->SetXTitle("Col");
    hshower_disp[p]->SetYTitle("Row");
    h_disp[p]->SetOption("COLZ");
    h_disp[p]->SetXTitle("Col");
    h_disp[p]->SetYTitle("Row");
    larcv::ImageMeta meta =wire_img.at(p).meta();

    // loop and set to one
    for (int c=1;c<3457;c++){
      for (int r=1;r<1009;r++){
          hshower_disp[p]->SetBinContent(c,r,1.0);
          h_disp[p]->SetBinContent(c,r,1.0);
      }
    }

    // loop over  shower Hits
    for (int i =0; i<shower_v.size();i++){
      larflow::reco::cluster_t  shower_c = larflow::reco::cluster_from_larflowcluster( shower_v[i] );
      std::cout<<"number of pts "<<shower_c.imgcoord_v.size()<<std::endl;
      for ( size_t ihit=0; ihit<shower_c.imgcoord_v.size(); ihit++ ) {
        const std::vector<int>& coord = shower_c.imgcoord_v[ihit];
        int tick = coord[3];
        int wire = coord[p];
        int row = meta.row(tick);
        int col = meta.col(wire);
        int thruval = thru_img[p].pixel(row,col);
        int adcval = wire_img[p].pixel(row,col);
        if (adcval >100 && thruval < 1 ){
          hshower_disp[p]->SetBinContent(col,row,100.0);
        }
        if(adcval <= 100 && thruval  < 1 ){
          hshower_disp[p]->SetBinContent(col,row,adcval);
        }
      }//end of loop over hits
    }//end of loop over showers

    // loop over all adc hits
    for (int c=0;c<3456;c++){
      for (int r=0;r<1008;r++){
          if (wire_img[p].pixel(r,c)>0.0 && thru_img[p].pixel(r,c) < 1){
            h_disp[p]->SetBinContent(c,r,wire_img[p].pixel(r,c));
          }
          if (wire_img[p].pixel(r,c)>100. && thru_img[p].pixel(r,c) < 1){
            h_disp[p]->SetBinContent(c,r,100.0);
          }
      }
    }

  }//end of loop over planes

  // get 3d location of reco vertex
  std::vector<double> reco_vertex = {(double)nufit_vertex_v[0],(double)nufit_vertex_v[1], (double)nufit_vertex_v[2]};
  // now project into 2d
  std::vector<int> vtx_rc = getProjectedPixel(reco_vertex, wire_img[0].meta(), 3);
  // fill root histograms
  hvtx_disp[0]->Fill(vtx_rc[1], vtx_rc[0]);
  hvtx_disp[1]->Fill(vtx_rc[2], vtx_rc[0]);
  hvtx_disp[2]->Fill(vtx_rc[3], vtx_rc[0]);

  // save canvas to png in correct folder
  // just save y plane for now
  TCanvas can("can", "histograms ", 3456, 1008);
  can.cd();
  // first draw adc images
  h_disp[2]->SetTitle(Form("Image Raw U Plane Run: %d Subrun: %d Event: %d",run,subrun,event));
  h_disp[2]->SetXTitle("Column (Wire)");
  h_disp[2]->SetYTitle("Row (6 Ticks)");
  h_disp[2]->SetOption("COLZ");
  h_disp[2]->Draw("");
  hvtx_disp[2]->SetMarkerStyle(kStar);
  hvtx_disp[2]->SetMarkerColor(2);
  hvtx_disp[2]->SetMarkerSize(5);
  hvtx_disp[2]->Draw("SAME");
  if (low==0){
    can.SaveAs(Form("/cluster/tufts/wongjiradlab/kmason03/ubdl-ana/showerenergy_study/outputshigh/adc_%d_%d_%d.png",run,subrun,event));
  }
  else{
    can.SaveAs(Form("/cluster/tufts/wongjiradlab/kmason03/ubdl-ana/showerenergy_study/outputslow/adc_%d_%d_%d.png",run,subrun,event));
  }

  // now draw shower images
  hshower_disp[2]->SetTitle(Form("Image Raw U Plane Run: %d Subrun: %d Event: %d",run,subrun,event));
  hshower_disp[2]->SetXTitle("Column (Wire)");
  hshower_disp[2]->SetYTitle("Row (6 Ticks)");
  hshower_disp[2]->SetOption("COLZ");
  hshower_disp[2]->Draw("");
  hvtx_disp[2]->SetMarkerStyle(kStar);
  hvtx_disp[2]->SetMarkerColor(2);
  hvtx_disp[2]->SetMarkerSize(5);
  hvtx_disp[2]->Draw("SAME");
  if (low==0){
    can.SaveAs(Form("/cluster/tufts/wongjiradlab/kmason03/ubdl-ana/showerenergy_study/outputshigh/shower_%d_%d_%d.png",run,subrun,event));
  }
  else{
    can.SaveAs(Form("/cluster/tufts/wongjiradlab/kmason03/ubdl-ana/showerenergy_study/outputslow/shower_%d_%d_%d.png",run,subrun,event));
  }



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

std::vector<int> getProjectedPixel( const std::vector<double>& pos3d,
            const larcv::ImageMeta& meta,
            const int nplanes,
            const float fracpixborder ) {

  // function that takes in 3d positions and returns time,wire_u,wire_v,wire_y in pixel coordinates
  std::vector<int> img_coords( nplanes+1, -1 );
  float row_border = fabs(fracpixborder)*meta.pixel_height();
  float col_border = fabs(fracpixborder)*meta.pixel_width();
  // tick/row
  float tick = pos3d[0]/(::larutil::LArProperties::GetME()->DriftVelocity()*::larutil::DetectorProperties::GetME()->SamplingRate()*1.0e-3) + 3200.0;
  if ( tick<meta.min_y() ) {
    if ( tick>meta.min_y()-row_border )
      // below min_y-border, out of image
      img_coords[0] = meta.rows()-1; // note that tick axis and row indicies are in inverse order (same order in larcv2)
    else
      // outside of image and border
      img_coords[0] = -1;
  }
  else if ( tick>meta.max_y() ) {
    if ( tick<meta.max_y()+row_border )
      // within upper border
      img_coords[0] = 0;
    else
      // outside of image and border
      img_coords[0] = -1;
  }
  else {
    // within the image
    img_coords[0] = meta.row( tick );
  }
  // Columns
  Double_t xyz[3] = { pos3d[0], pos3d[1], pos3d[2] };
  // there is a corner where the V plane wire number causes an error
  if ( (pos3d[1]>-117.0 && pos3d[1]<-116.0) && pos3d[2]<2.0 ) {
    xyz[1] = -116.0;
  }
  for (int p=0; p<nplanes; p++) {
    float wire = larutil::Geometry::GetME()->WireCoordinate( xyz, p );
    // get image coordinates
    if ( wire<meta.min_x() ) {
      if ( wire>meta.min_x()-col_border ) {
  // within lower border
  img_coords[p+1] = 0;
      }
      else
  img_coords[p+1] = -1;
    }
    else if ( wire>=meta.max_x() ) {
      if ( wire<meta.max_x()+col_border ) {
  // within border
  img_coords[p+1] = meta.cols()-1;
      }
      else
  // outside border
  img_coords[p+1] = -1;
    }
    else
      // inside image
      img_coords[p+1] = meta.col( wire );
  }//end of plane loop
  // there is a corner where the V plane wire number causes an error
  if ( pos3d[1]<-116.3 && pos3d[2]<2.0 && img_coords[1+1]==-1 ) {
    img_coords[1+1] = 0;
  }
  return img_coords;
}//end of project pixel function

void printinfo(larlite::event_mctruth* ev_mctruth){
    // get event info
    // Neutrino energy
   float _enu_true=0.0;
   _enu_true = ev_mctruth->at(0).GetNeutrino().Nu().Momentum().E()*1000.0;
   std::cout<<"Neutrino Energy: "<<_enu_true<<std::endl;

   // cc or nc event?
   bool ccevent = false;
   int ccnc = ev_mctruth->at(0).GetNeutrino().CCNC();
   if (ccnc ==0 ) ccevent=true;
   if (ccevent) std::cout<<"Is a CC Event"<<std::endl;
   else std::cout<<"Is a NC Event"<<std::endl;

   // type of Neutrino
   int nu_pdg = ev_mctruth->at(0).GetNeutrino().Nu().PdgCode();
   if (nu_pdg== 12) std::cout<<"Muon Neutrino event "<<std::endl;
   else if (nu_pdg== -12) std::cout<<"Muon Anti Neutrino event "<<std::endl;
   else if (nu_pdg== 14) std::cout<<"Electron Neutrino event "<<std::endl;
   else if (nu_pdg== -14) std::cout<<"Electon Anti Neutrino event "<<std::endl;

   // type of interction - see comments at end of script
   int int_type= ev_mctruth->at(0).GetNeutrino().InteractionType();
   if (int_type == 1001 || int_type == 1002) std::cout<<"QE Interaction "<<std::endl;
   else if (int_type >= 1003 && int_type <= 1090) std::cout<<"RES Interaction "<<std::endl;
   else if (int_type == 1092 || int_type == 1091) std::cout<<"DIS Interaction "<<std::endl;
   else std::cout<<"Other Interaction: "<<int_type<<std::endl;

   int num_protons = 0;
   int num_neutrons = 0;
   int num_pion_charged = 0;
   int num_pion_neutral = 0;
   for(int part =0;part<(int)ev_mctruth->at(0).GetParticles().size();part++){
     // pick only final state particles
     if (ev_mctruth->at(0).GetParticles().at(part).StatusCode() == 1){
       int pdg = ev_mctruth->at(0).GetParticles().at(part).PdgCode();
       if (pdg == 2212) num_protons++;
       else if (pdg == 2112) num_neutrons++;
       else if (pdg == 111 ) num_pion_neutral++;
       else if (pdg == 211 || pdg == -211) num_pion_charged++;

     }//end of if status = 1 statement
   }//end of loop over particles

   std::cout<<"Number of protons: "<<num_protons<<std::endl;
   std::cout<<"Number of neutrons: "<<num_neutrons<<std::endl;
   std::cout<<"Number of charged pions: "<<num_pion_charged<<std::endl;
   std::cout<<"Number of neutral pions: "<<num_pion_neutral<<std::endl;

  std::cout << "\n";
  return;
}
