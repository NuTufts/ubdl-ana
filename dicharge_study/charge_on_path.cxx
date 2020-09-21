#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "DataFormat/storage_manager.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larflow/Reco/TrackTruthRecoAna.h"

#include "TTree.h"

int main( int nargs, char** argv ) {

  std::string input_larcv_adc      = argv[1];
  std::string input_larlite_mcinfo = argv[2];
  std::string output_anafile       = argv[3];

  const float max_step_size = 0.1;
  
  larlite::storage_manager ioll( larlite::storage_manager::kREAD );
  ioll.add_in_filename( input_larlite_mcinfo );
  ioll.set_data_to_read( larlite::data::kMCTrack,  "mcreco" );
  ioll.open();

  larcv::IOManager iolcv( larcv::IOManager::kREAD, "iolcv", larcv::IOManager::kTickBackward );
  iolcv.add_in_file( input_larcv_adc );
  iolcv.specify_data_read( larcv::kProductImage2D, "wire" );
  iolcv.reverse_all_products();
  iolcv.initialize();

  larutil::SpaceChargeMicroBooNE* psce = new larutil::SpaceChargeMicroBooNE;

  int nentries = ioll.get_entries();

  TFile* out = new TFile( output_anafile.c_str(), "recreate" );
  TTree* anatree = new TTree("ana","ana");
  int pid;
  float phi;
  float qperpix[3];
  anatree->Branch("pid",&pid,"pid/I");
  anatree->Branch("phi",&phi,"phi/F");
  anatree->Branch("qperpix", qperpix, "qperpix[3]/F");


  for (int ientry=0; ientry<nentries; ientry++) {

    std::cout << "=== [ENTRY " << ientry << "] ===" << std::endl;
    
    ioll.go_to(ientry);
    iolcv.read_entry(ientry);

    larlite::event_mctrack* ev_mctrack
      = (larlite::event_mctrack*)ioll.get_data(larlite::data::kMCTrack,"mcreco");

    larcv::EventImage2D* evimg
      = (larcv::EventImage2D*)iolcv.get_data(larcv::kProductImage2D,"wire");
    auto& adc_v = evimg->Image2DArray();


    for (int itrack=0; itrack<(int)ev_mctrack->size(); itrack++) {

      auto const& mct = ev_mctrack->at(itrack);

      if ( mct.PdgCode()!=13 && mct.PdgCode()!=-13 && mct.PdgCode()!=2212)
        continue;

      pid = mct.PdgCode();
    
      std::vector< std::vector<float> > true_path
        = larflow::reco::TrackTruthRecoAna::getSCEtrueTrackPath( mct, psce );

      if ( true_path.size()<2 )
        continue;

      std::cout << "track[" << itrack << "] nsteps=" << true_path.size() << std::endl;
          
      int icurr_row = -1;
      std::vector<int> icurr_col(3,-1);
      
      for (int i=0; i<(int)true_path.size()-1; i++) {
        const std::vector<float>& pt1 = true_path[i];
        const std::vector<float>& pt2 = true_path[i+1];

        float steplen = 0.;
        std::vector<float> dir(3,0);
        for (int v=0; v<3; v++) {
          dir[v] = pt2[v]-pt1[v];
          steplen += dir[v]*dir[v];
        }
        steplen = sqrt(steplen);

        if ( steplen==0 )
          continue;

        for (int v=0; v<3; v++)
          dir[v] /= steplen;

        int nsubsteps = steplen/max_step_size+1;
        float substep = steplen/float(nsubsteps);

        for (int istep=0; istep<nsubsteps; istep++) {
          std::vector<double> imgpt(3,0);
          for (int v=0; v<3; v++)
            imgpt[v] = (double)pt1[v] + (double)substep*istep*dir[v];

          std::vector<float> imgcoord(4,0);
          imgcoord[3] = imgpt[0]/larutil::LArProperties::GetME()->DriftVelocity()/0.5+3200;
          if ( imgcoord[3]<2400.0 || imgcoord[3]>=adc_v.front().meta().max_y() )
            continue;

          int row = adc_v.front().meta().row( imgcoord[3] );

          bool repeat = true;
          if ( icurr_row!=row )
            repeat = false;
          for (int p=0; p<3; p++) {
            if ( icurr_col[p]!=(int)imgcoord[p] )
              repeat = false;
          }

          if ( repeat )
            continue;

          icurr_row = row;
          for (int p=0; p<3; p++)
            icurr_col[p] = (int)imgcoord[p];
          
          // sum a 3 pix window
          std::vector<float> totq(3,0);
          std::vector<float> totpix(3,0);
          for (int dr=-1; dr<=1; dr++) {
            int r = row+dr;
            if ( r<0 || r>=(int)adc_v.front().meta().rows() )
              continue;
          
            for (int p=0; p<3; p++) {
              imgcoord[p] = larutil::Geometry::GetME()->WireCoordinate( imgpt, p );
              if ( imgcoord[p]<0 || imgcoord[p]>=(float)adc_v[p].meta().max_x() )
                continue;
              
              float q = adc_v[p].pixel(r,(int)imgcoord[p]);
              if ( q>5.0 ) {
                totq[p] += q;
                totpix[p] += 1.0;
              }
            }
          }//end of dr loop
          
          phi = atan2(dir[1],dir[0]);
          for (int p=0; p<3; p++) {
            if ( totpix[p]>0 )
              qperpix[p] = totq[p]/totpix[p];
            else
              qperpix[p] = 0.;
          }
          
          anatree->Fill();
        }
      }
      
    }
      
  }

  std::cout << "[[WRITE OUTPUT]]" << std::endl;
  anatree->Write();

  std::cout << "[[CLOSE OUTPUT]]" << std::endl;
  out->Close();

  std::cout << "[[CLOSE INPUT]]" << std::endl;
  iolcv.finalize();
  ioll.close();
  
  std::cout << "[[FIN]]" << std::endl;
  
  return 0;
}
