#include <iostream>
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/larflow3dhit.h"
#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"
#include "larflow/Reco/EventKeypointReco.h"
#include "ublarcvapp/MCTools/MCPixelPGraph.h"
#include "ublarcvapp/MCTools/NeutrinoVertex.h"
#include "ublarcvapp/MCTools/LArbysMC.h"

#include "TTree.h"

int main( int nargs, char** argv ) {
  
  std::cout << "Analyze Larmatch ME keypoint reconstruction performance." << std::endl;

  std::string finput_larmatch = argv[1];
  std::string finput_mcinfo   = argv[2];
  std::string foutput = argv[3];

  larlite::storage_manager io( larlite::storage_manager::kBOTH );
  io.set_data_to_read(  larlite::data::kLArFlow3DHit, "larmatch" );
  io.set_data_to_read(  larlite::data::kMCTruth, "generator" );
  io.set_data_to_read(  larlite::data::kMCTrack, "mcreco" );  
  io.set_data_to_read(  larlite::data::kMCShower, "mcreco" );    
  io.set_data_to_write( larlite::data::kLArFlow3DHit, "keypoint" );
  io.set_data_to_write( larlite::data::kPCAxis, "keypoint" );
  io.add_in_filename( finput_larmatch );
  io.add_in_filename( finput_mcinfo );  
  io.set_out_filename( foutput );
  io.open();
  int nentries = io.get_entries();
  
  larflow::reco::EventKeypointReco algo;
  ublarcvapp::mctools::NeutrinoVertex nuvtx;
  larutil::SpaceChargeMicroBooNE sce;
  ublarcvapp::mctools::LArbysMC larbysmc;

  std::string outbase = foutput.substr( 0, foutput.find(".root") );
  std::cout << "outbase: " << outbase << std::endl;
  outbase += "_ana.root";
  TFile* outroot = new TFile(outbase.c_str(),"recreate");
  TTree* kprecoana = new TTree("kprecoana", "Keypoint Reco analysis per keypoint");
  TTree* kprecoana_event = new TTree("kprecoana_event", "Keypoint Reco analysis per Event");
  larbysmc.bindAnaVariables( kprecoana_event );

  float closest_nu_kp_to_truth_vtx = 1999.0;
  int num_fp = 0;
  int num_tp = 0;
  kprecoana_event->Branch("closest_to_truth_cm", &closest_nu_kp_to_truth_vtx, "closest_to_truth_cm/F" );
  kprecoana_event->Branch("num_fp", &num_fp, "num_fp/I" );
  kprecoana_event->Branch("num_tp", &num_tp, "num_tp/I" );  
  
  
  io.next_event();
  for (int ientry=0; ientry<nentries; ientry++) {
    
    std::cout << "[ENTRY " << ientry << "] ==================================" << std::endl;    
    //io.go_to(ientry);
    
    // reset variables
    closest_nu_kp_to_truth_vtx = 1999.0;
    num_fp = 0;
    num_tp = 0;    

    ublarcvapp::mctools::MCPixelPGraph mcpg;
    mcpg.buildgraphonly(io);

    larbysmc.process(io);

    // get true neutrino xyz (with space-charge applied)
    std::vector<float> vtxsce = nuvtx.getPos3DwSCE( io, &sce );
    
    auto ev_lfhits = (larlite::event_larflow3dhit*)io.get_data( larlite::data::kLArFlow3DHit, "larmatch" );
    std::cout << "num of larflow spacepoints: " << ev_lfhits->size() << std::endl;
    algo.process_larmatch_v2( io, "larmatch" );

    auto ev_kp = (larlite::event_larflow3dhit*)io.get_data( larlite::data::kLArFlow3DHit, "keypoint" );
    std::cout << "num of reco keypoints: " << ev_kp->size() << std::endl;

    
    for ( auto kp : *ev_kp ) {
      if ( int(kp[3])!=0 )
	continue;
      float dist = 0.;
      for (int i=0; i<3; i++)
	dist += (kp[i]-vtxsce[i])*(kp[i]-vtxsce[i]);
      dist = sqrt(dist);
      if ( dist < closest_nu_kp_to_truth_vtx )
	closest_nu_kp_to_truth_vtx = dist;
      if ( dist>2.0 )
	num_fp++;
    }
    std::cout << "  closest distance of nu keypoint to truth: " << closest_nu_kp_to_truth_vtx << " cm"  << std::endl;
    if ( closest_nu_kp_to_truth_vtx<2.0 )
      num_tp++;
    
    io.set_id( io.run_id(), io.subrun_id(), io.event_id() );
    io.next_event();    

    kprecoana_event->Fill();

    if ( ientry+1>=50 )
      break;
  }
  io.close();
  outroot->Write();
  
  return 0;
}
