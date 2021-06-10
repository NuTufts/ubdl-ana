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
std::vector<float> GetADCSumWithSmear(std::vector<larlite::larflowcluster> shower_v,
            int shower1,int shower2,std::vector<larcv::Image2D> wire_img,int threshold );
std::vector<float> GetShowerDir(std::vector<larlite::pcaxis> shower_v,int shower1,int shower2);
std::vector<float> GetShowerStart(std::vector<larlite::track> showertrunk,int shower1,int shower2);
std::vector<float> GetShowerEnd(std::vector<larlite::track> showertrunk,int shower1,int shower2);
int PickShower(std::vector<larlite::larflowcluster> shower_v);
std::vector<float> PickSecondShower(std::vector<larlite::larflowcluster> shower_v, int leadshower);
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
  //script to load in gen2 reco files, and get pi0 reco variables
  std::cout << "Gen-2 Pi0 Reco Test" << std::endl;

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

  // load in gen2 variables
  int npi0;
  in->SetBranchAddress( "npi0", &npi0 );
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

  // setup output variables
  std::string anatreename = "pi0reco_anatree";
  int _run;
  int _subrun;
  int _event;
  int _ispi0;
  int _best_passing_vtx_index;
  std::vector<float> _leading_shower_ADCU;
  std::vector<float> _leading_shower_ADCV;
  std::vector<float> _leading_shower_ADCY;
  std::vector<float> _subleading_shower_ADCU;
  std::vector<float> _subleading_shower_ADCV;
  std::vector<float> _subleading_shower_ADCY;
  std::vector<float> _leading_shower_Xdir;
  std::vector<float> _subleading_shower_Xdir;
  std::vector<float> _leading_shower_Ydir;
  std::vector<float> _subleading_shower_Ydir;
  std::vector<float> _leading_shower_Zdir;
  std::vector<float> _subleading_shower_Zdir;
  std::vector<float> _leading_shower_Xstart;
  std::vector<float> _subleading_shower_Xstart;
  std::vector<float> _leading_shower_Ystart;
  std::vector<float> _subleading_shower_Ystart;
  std::vector<float> _leading_shower_Zstart;
  std::vector<float> _subleading_shower_Zstart;
  std::vector<float> _leading_shower_Xend;
  std::vector<float> _subleading_shower_Xend;
  std::vector<float> _leading_shower_Yend;
  std::vector<float> _subleading_shower_Yend;
  std::vector<float> _leading_shower_Zend;
  std::vector<float> _subleading_shower_Zend;
  std::vector<float> _overlap_frac;
  TTree* _ana_tree = new TTree(anatreename.c_str(), "Pi0 reco variables");
  _ana_tree->Branch("Run",&_run);
  _ana_tree->Branch("Subrun",&_subrun);
  _ana_tree->Branch("Event",&_event);
  _ana_tree->Branch("ispi0",&_ispi0);
  _ana_tree->Branch("best_passing_vtx_index",&_best_passing_vtx_index);
  _ana_tree->Branch("overlap_frac",&_overlap_frac);
  _ana_tree->Branch("leading_shower_ADCU",&_leading_shower_ADCU);
  _ana_tree->Branch("leading_shower_ADCV",&_leading_shower_ADCV);
  _ana_tree->Branch("leading_shower_ADCY",&_leading_shower_ADCY);
  _ana_tree->Branch("subleading_shower_ADCU",&_subleading_shower_ADCU);
  _ana_tree->Branch("subleading_shower_ADCV",&_subleading_shower_ADCV);
  _ana_tree->Branch("subleading_shower_ADCY",&_subleading_shower_ADCY);
  _ana_tree->Branch("leading_shower_Xdir",&_leading_shower_Xdir);
  _ana_tree->Branch("subleading_shower_Xdir",&_subleading_shower_Xdir);
  _ana_tree->Branch("leading_shower_Ydir",&_leading_shower_Ydir);
  _ana_tree->Branch("subleading_shower_Ydir",&_subleading_shower_Ydir);
  _ana_tree->Branch("leading_shower_Zdir",&_leading_shower_Zdir);
  _ana_tree->Branch("subleading_shower_Zdir",&_subleading_shower_Zdir);
  _ana_tree->Branch("leading_shower_Xstart",&_leading_shower_Xstart);
  _ana_tree->Branch("subleading_shower_Xstart",&_subleading_shower_Xstart);
  _ana_tree->Branch("leading_shower_Ystart",&_leading_shower_Ystart);
  _ana_tree->Branch("subleading_shower_Ystart",&_subleading_shower_Ystart);
  _ana_tree->Branch("leading_shower_Zstart",&_leading_shower_Zstart);
  _ana_tree->Branch("subleading_shower_Zstart",&_subleading_shower_Zstart);
  _ana_tree->Branch("leading_shower_Xend",&_leading_shower_Xend);
  _ana_tree->Branch("subleading_shower_Xend",&_subleading_shower_Xend);
  _ana_tree->Branch("leading_shower_Yend",&_leading_shower_Yend);
  _ana_tree->Branch("subleading_shower_Yend",&_subleading_shower_Yend);
  _ana_tree->Branch("leading_shower_Zend",&_leading_shower_Zend);
  _ana_tree->Branch("subleading_shower_Zend",&_subleading_shower_Zend);

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


  //loop through entries
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
    if (npi0>0) _ispi0 =  1;
    else _ispi0 = 0;

    // load larcv/larlite branches
    const auto ev_img = (larcv::EventImage2D*)io_larcv->get_data( larcv::kProductImage2D, "wire" );
    const auto ev_thru = (larcv::EventImage2D*)io_larcv->get_data( larcv::kProductImage2D, "thrumu" );
    _run = ev_img->run();
		_subrun = ev_img->subrun();
		_event = ev_img->event();
    auto const& wire_img = ev_img->Image2DArray();
    auto const& thru_img = ev_thru->Image2DArray();
    larlite::event_mcshower* ev_mcshower = (larlite::event_mcshower*)ioll->get_data(larlite::data::kMCShower,"mcreco");
    larcv::EventPGraph* ev_test_pgraph  = (larcv::EventPGraph*)(io_larcv->get_data(larcv::kProductPGraph,"test"));
    int nvtx = (int)(*pnu_sel_v).size();

    std::vector<float> adcaverage_v;
    std::vector<double> dist_v;

    std::vector<EventDist2True_t> index_by_dist_v;
    std::cout<<nvtx<<std::endl;
    for (int ivtx=0; ivtx<nvtx; ivtx++) {
      auto const& nusel = (*pnu_sel_v)[ivtx];

      EventDist2True_t idx( nusel.dist2truevtx, ivtx );
      index_by_dist_v.push_back( idx );

      // get reco shower clusters!
      std::vector<larlite::larflowcluster> nufit_shower_v = (*nufitted_v)[ivtx].shower_v;
      std::vector<larlite::pcaxis> nufit_shower_pca_v = (*nufitted_v)[ivtx].shower_pcaxis_v;
      std::vector<larlite::track> nufit_shower_trunk_v = (*nufitted_v)[ivtx].shower_trunk_v;

      // pick largest shower first.
      // then pick largest, non overlapping
      int leadshower=PickShower(nufit_shower_v);
      std::vector<float> subshower_v=PickSecondShower(nufit_shower_v, leadshower);
      int subshower = int(subshower_v[0]);
      _overlap_frac.push_back(subshower_v[1]);
      // std::cout<<"CHOSEN SHOWERS: "<<leadshower<<" , "<<subshower<<std::endl;

      // get sum of ADC in cluster
      std::vector<float> ADCSum_v;
      std::vector<float> dir_v;
      std::vector<float> start_v;
      std::vector<float> end_v;
      float adcaverage;


      ADCSum_v = GetADCSumWithSmear(nufit_shower_v,leadshower,subshower, wire_img, 10 );
      dir_v = GetShowerDir(nufit_shower_pca_v,leadshower,subshower);
      // need shower start and end for now, later pick which is closest to the vertex
      start_v = GetShowerStart(nufit_shower_trunk_v,leadshower,subshower);
      end_v = GetShowerStart(nufit_shower_trunk_v,leadshower,subshower);
      adcaverage = (ADCSum_v[0]+ADCSum_v[1]+ADCSum_v[2])/(3.0);

      _leading_shower_ADCU.push_back(ADCSum_v[0]);
      _leading_shower_ADCV.push_back(ADCSum_v[1]);
      _leading_shower_ADCY.push_back(ADCSum_v[2]);
      _subleading_shower_ADCU.push_back(ADCSum_v[3]);
      _subleading_shower_ADCV.push_back(ADCSum_v[4]);
      _subleading_shower_ADCY.push_back(ADCSum_v[5]);
      adcaverage_v.push_back(adcaverage);
      _leading_shower_Xdir.push_back(dir_v[0]);
      _subleading_shower_Xdir.push_back(dir_v[3]);
      _leading_shower_Ydir.push_back(dir_v[1]);
      _subleading_shower_Ydir.push_back(dir_v[4]);
      _leading_shower_Zdir.push_back(dir_v[2]);
      _subleading_shower_Zdir.push_back(dir_v[5]);
      _leading_shower_Xstart.push_back(start_v[0]);
      _subleading_shower_Xstart.push_back(start_v[3]);
      _leading_shower_Ystart.push_back(start_v[1]);
      _subleading_shower_Ystart.push_back(start_v[4]);
      _leading_shower_Zstart.push_back(start_v[2]);
      _subleading_shower_Zstart.push_back(start_v[5]);
      _leading_shower_Xend.push_back(end_v[0]);
      _subleading_shower_Xend.push_back(end_v[3]);
      _leading_shower_Yend.push_back(end_v[1]);
      _subleading_shower_Yend.push_back(end_v[4]);
      _leading_shower_Zend.push_back(end_v[2]);
      _subleading_shower_Zend.push_back(end_v[5]);



    }//end of vtx loop

    std::sort( index_by_dist_v.begin(), index_by_dist_v.end() );
    if ( index_by_dist_v.size()>0 ) {
      _best_passing_vtx_index = index_by_dist_v.front().index;
    }

    _ana_tree->Fill();
    _leading_shower_ADCU.clear();
    _leading_shower_ADCV.clear();
    _leading_shower_ADCY.clear();
    _subleading_shower_ADCU.clear();
    _subleading_shower_ADCV.clear();
    _subleading_shower_ADCY.clear();
    _leading_shower_Xdir.clear();
    _subleading_shower_Xdir.clear();
    _leading_shower_Ydir.clear();
    _subleading_shower_Ydir.clear();
    _leading_shower_Zdir.clear();
    _subleading_shower_Zdir.clear();
    _leading_shower_Xstart.clear();
    _subleading_shower_Xstart.clear();
    _leading_shower_Ystart.clear();
    _subleading_shower_Ystart.clear();
    _leading_shower_Zstart.clear();
    _subleading_shower_Zstart.clear();
    _leading_shower_Xend.clear();
    _subleading_shower_Xend.clear();
    _leading_shower_Yend.clear();
    _subleading_shower_Yend.clear();
    _leading_shower_Zend.clear();
    _subleading_shower_Zend.clear();
    _overlap_frac.clear();


  }//end of entry loop

  _ana_tree->Write();
  out->Write();
  io_larcv->finalize();
  ioll->close();
  delete in;

  return 0;
}

std::vector<float> GetShowerDir(std::vector<larlite::pcaxis> shower_v,int shower1,int shower2){
  // get the direction of each the first and second shower
  // initialize output
  std::vector<float> dir_v = {-999,-999,-999,-999,-999,-999};

  // get EigenVectors
  if (shower1 >=0){
    std::vector<std::vector<double> > showereigen_v = shower_v[shower1].getEigenVectors();
    dir_v[0] = float(showereigen_v[0][0]);
    dir_v[1] = float(showereigen_v[0][1]);
    dir_v[2] = float(showereigen_v[0][2]);
  }
  if (shower1 >=0 && shower2 >=0){
    std::vector<std::vector<double> > showereigen2_v = shower_v[shower2].getEigenVectors();
    dir_v[3] = float(showereigen2_v[0][0]);
    dir_v[4] = float(showereigen2_v[0][1]);
    dir_v[5] = float(showereigen2_v[0][2]);
  }


  return dir_v;
}

std::vector<float> GetShowerStart(std::vector<larlite::track> showertrunk,int shower1,int shower2){
  // get the distance between the true and reco shower start
  // initialize output
  std::vector<float> start_v={-999,-999,-999,-999,-999,-999};
  if (shower1 >=0){
    // get  lead reco shower start
    std::vector<double> recoshowerstart;
    TVector3 reco_start_v = showertrunk[shower1].LocationAtPoint(0);
    recoshowerstart={reco_start_v.X(),reco_start_v.Y(),reco_start_v.Z()};
    start_v[0]=float(recoshowerstart[0]);
    start_v[1]=float(recoshowerstart[1]);
    start_v[2]=float(recoshowerstart[2]);

    if (shower2>=0){
      // get  lead reco shower start
      std::vector<double> recoshowerstart2;
      TVector3 reco_start_v2 = showertrunk[shower2].LocationAtPoint(0);
      recoshowerstart2={reco_start_v2.X(),reco_start_v2.Y(),reco_start_v2.Z()};
      start_v[3]=float(recoshowerstart2[0]);
      start_v[4]=float(recoshowerstart2[1]);
      start_v[5]=float(recoshowerstart2[2]);

    }
  }



  return start_v;
}

std::vector<float> GetShowerEnd(std::vector<larlite::track> showertrunk,int shower1,int shower2){
  // get the distance between the true and reco shower start
  // initialize output
  std::vector<float> start_v={-999,-999,-999,-999,-999,-999};
  if (shower1 >=0){
    // get  lead reco shower end
    std::vector<double> recoshowerstart;
    TVector3 reco_start_v = showertrunk[shower1].End();
    recoshowerstart={reco_start_v.X(),reco_start_v.Y(),reco_start_v.Z()};
    start_v[0]=float(recoshowerstart[0]);
    start_v[1]=float(recoshowerstart[1]);
    start_v[2]=float(recoshowerstart[2]);

    if(shower2 >=0){
      // get  sub reco shower end
      std::vector<double> recoshowerstart2;
      TVector3 reco_start_v2 = showertrunk[shower2].End();
      recoshowerstart2={reco_start_v2.X(),reco_start_v2.Y(),reco_start_v2.Z()};
      start_v[3]=float(recoshowerstart2[0]);
      start_v[4]=float(recoshowerstart2[1]);
      start_v[5]=float(recoshowerstart2[2]);
      return start_v;
    }
  }
}


std::vector<float> GetADCSumWithSmear(std::vector<larlite::larflowcluster> shower_v,int shower1,int shower2,
              std::vector<larcv::Image2D> wire_img,int threshold ){
  //initialize output (u1,v1,y1,u2,v2,y2)
  std::vector<float> sum_v = {-999,-999,-999,-999,-999,-999};
  if (shower1 >=0){
    bool usesec = false;
    if (shower2 >=0 && shower1 !=shower2 ) usesec = true;
    // / first turn to cluster
    //save all three planes
    std::vector<larlite::larflow3dhit>  shower_c = shower_v[shower1];
    std::vector<larlite::larflow3dhit>  secshower_c;
    if (shower2 >=0)  secshower_c=shower_v[shower2];

    // turn into 2d points (u,v,y,t)
    auto const& wireu_meta = wire_img.at(0).meta();
    auto const& wirev_meta = wire_img.at(1).meta();
    auto const& wirey_meta = wire_img.at(2).meta();
    // loop over planes
    for (int p =0;p<3;p++){
      larcv::ImageMeta meta =wire_img.at(p).meta();
      float sum =0;
      float sum2 =0;
      // make found list so we don't use the same point multiple times vector<(row,col)>
      // std::vector<std::vector<int>> foundpts;
      // loop over Hits
      float noutpix1 =0.0;
      float noutpix2 =0.0;
      // initialize alredy used points so we don't double count
      larcv::Image2D usedpixels = wire_img[p];

      // loop and set to zero
      for (int c=1;c<3456;c++){
        for (int r=1;r<1008;r++){
            usedpixels.set_pixel(r,c,0);
            if (wire_img[p].pixel(r,c) <threshold) wire_img[p].set_pixel(r,c,0);
        }
      }

      for ( size_t ihit=0; ihit<shower_c.size(); ihit++ ) {
        int wire = shower_c[ihit].targetwire[p];
        int tick =  shower_c[ihit].tick;
        // add a check for projection failures
        if (tick < 0 || wire <0 || wire >=3456 ){
          noutpix1 +=1.0;
          continue;
        }
        int row = -1;
        int col = -1;
        if (wire < 3456 && tick >0 && wire >0){
          row = meta.row(tick);
          col = meta.col(wire);
        }
        float adcval =0;
        // std::cout<<row<<" "<<col<<std::endl;
        if (col < 3455 && row <1007 && col > 1 && row > 1){
          // for each pixel, check if in mask
          if (usedpixels.pixel(row,col)==0){
              usedpixels.set_pixel(row,col,1);
              if (!usesec) adcval+=wire_img[p].pixel(row,col);
          }

          // now get surrounding pixels
          // [(r+1,c-1),(r+1,c),(r+1,c+1),(r,c-1),(r,c+1),(r-1,c-1),(r-1,c),(r-1,c+1)]
          // if (usedpixels.pixel(row+1,col-1)==0){
          //     adcval += wire_img[p].pixel(row+1,col-1);
          //     usedpixels.set_pixel(row+1,col-1,1);
          // }
          // if (usedpixels.pixel(row+1,col)==0){
          //     adcval += wire_img[p].pixel(row+1,col);
          //     usedpixels.set_pixel(row+1,col,1);
          // }
          // if (usedpixels.pixel(row+1,col+1)==0){
          //     adcval += wire_img[p].pixel(row+1,col+1);
          //     usedpixels.set_pixel(row+1,col+1,1);
          // }
          // if (usedpixels.pixel(row,col-1)==0){
          //     adcval += wire_img[p].pixel(row,col-1);
          //     usedpixels.set_pixel(row,col-1,1);
          // }
          // if (usedpixels.pixel(row,col+1)==0){
          //     adcval += wire_img[p].pixel(row,col+1);
          //     usedpixels.set_pixel(row,col+1,1);
          // }
          // if (usedpixels.pixel(row-1,col-1)==0){
          //     adcval += wire_img[p].pixel(row-1,col-1);
          //     usedpixels.set_pixel(row-1,col-1,1);
          // }
          // if (usedpixels.pixel(row-1,col)==0){
          //     adcval += wire_img[p].pixel(row-1,col);
          //     usedpixels.set_pixel(row-1,col,1);
          // }
          // if (usedpixels.pixel(row-1,col+1)==0){
          //     adcval += wire_img[p].pixel(row-1,col+1);
          //     usedpixels.set_pixel(row-1,col+1,1);
          // }
        }
        sum+=adcval;
      }//end of loop over Hits
      //loop over second shower
      if (usesec){
        for ( size_t ihit=0; ihit<secshower_c.size(); ihit++ ) {
          int wire = secshower_c[ihit].targetwire[p];
          int tick = secshower_c[ihit].tick;
          // add a check for projection failures
          if (tick < 0 || wire <0 || wire >=3456 ){
            noutpix2 +=1.0;
            continue;
          }
          int row = -1;
          int col = -1;
          if (wire < 3456 && tick >0 && wire >0){
            row = meta.row(tick);
            col = meta.col(wire);
          }
          float adcval =0;
          // std::cout<<row<<" "<<col<<std::endl;
          if (col < 3455 && row <1007 && col > 1 && row > 1){
            // for each pixel, check if in mask
            // not used in first shower
            if (usedpixels.pixel(row,col)==0){
                usedpixels.set_pixel(row,col,2);
                adcval+=wire_img[p].pixel(row,col);
            }
            // overlap with first shower
            else if (usedpixels.pixel(row,col)==1){
                usedpixels.set_pixel(row,col,3);
                adcval+=(.5 * wire_img[p].pixel(row,col));
            }
          }
          sum2+=adcval;
        }//end of loop over Hits

      // now loop back over the first shower one more time, full adc if not in second, half otherwise

        for ( size_t ihit=0; ihit<shower_c.size(); ihit++ ) {
          int wire = shower_c[ihit].targetwire[p];
          int tick = shower_c[ihit].tick;
          // add a check for projection failures
          int row = -1;
          int col = -1;
          if (wire < 3456 && tick >0 && wire >0){
            row = meta.row(tick);
            col = meta.col(wire);
          }
          float adcval =0;
          // std::cout<<row<<" "<<col<<std::endl;
          if (col < 3455 && row <1007 && col > 1 && row > 1){
            // for each pixel, check if in mask
            // not used in first shower
            if (usedpixels.pixel(row,col)==1){
                adcval+=wire_img[p].pixel(row,col);
            }
            // overlap with first shower
            else if (usedpixels.pixel(row,col)==3){
                adcval+=(.5 * wire_img[p].pixel(row,col));
            }
          }
          sum+=adcval;
        }//end of loop over Hits

      }


      if ((1.0-float(noutpix1)/float(shower_c.size())) > .98) sum_v[p] = sum;
      if (usesec && (1.0-float(noutpix2)/float(secshower_c.size())) > .98 ) sum_v[p+3]=sum2;


    }//end of loop over planes

  }
  return sum_v;
}

int PickShower(std::vector<larlite::larflowcluster> shower_v){
  int bestshower = -999;
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

std::vector<float> PickSecondShower(std::vector<larlite::larflowcluster> shower_v, int leadshower){
  // loop through showers and pick second one
  std::vector<float> secshower = {-999,-999};
  int maxpix =0;
  if (leadshower >= 0){
    // get list of pixels in lead shower
    std::vector<std::vector<int>> leadpts;
    larflow::reco::cluster_t  shower_c = larflow::reco::cluster_from_larflowcluster( shower_v[leadshower] );
    for ( size_t ihit=0; ihit<shower_c.imgcoord_v.size(); ihit++ ) {
      std::vector<int> coord = shower_c.imgcoord_v[ihit];
      leadpts.push_back(coord);
    }//end of loop over hits
    // now loop through remaining showers
    for (int i =0; i<shower_v.size();i++){
      if (i != leadshower){
        int overlappix = 0;
        int totalpix = 0;
        // loop over hits and check for amount of overlap
        larflow::reco::cluster_t  shower_2 = larflow::reco::cluster_from_larflowcluster( shower_v[i] );
        for ( size_t ihit=0; ihit < shower_2.imgcoord_v.size(); ihit++ ) {
          std::vector<int> coord2 = shower_2.imgcoord_v[ihit];
          // now loop through lead to search for overlap
          bool overlap = false;
          for (int pt =0; pt < leadpts.size();pt++){
            if(leadpts[pt][0]==coord2[0] && leadpts[pt][1]==coord2[1]){
              overlap = true;
            }
          }//end of loop over lead hits
          totalpix+=1;
          if (overlap) overlappix +=1;
        }//end of loop over hits
        float overlap_frac = float(overlappix)/float(totalpix);
        if (overlap_frac < .02 && totalpix >maxpix){
          maxpix=totalpix;
          secshower[0] = i;
          secshower[1] = overlap_frac;

        }
      }
    } // end of loop over showers

  }// end of check over first
  return secshower;
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
