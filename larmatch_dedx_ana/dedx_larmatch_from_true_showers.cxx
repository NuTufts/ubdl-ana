#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"

#include "DataFormat/storage_manager.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/larflow3dhit.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

#include "larflow/Reco/geofuncs.h"
#include "ublarcvapp/MCTools/crossingPointsAnaMethods.h"
#include "ublarcvapp/ubdllee/dwall.h"

int main( int nargs, char** argv )
{

  std::string dlmerged_input_file = argv[1];
  std::string larmatch_input_file = argv[2];
  std::string output_file = argv[3];

  larlite::storage_manager io( larlite::storage_manager::kREAD );
  io.add_in_filename( dlmerged_input_file );
  io.add_in_filename( larmatch_input_file );
  io.set_data_to_read( larlite::data::kLArFlow3DHit, "larmatch" );
  io.set_data_to_read( larlite::data::kMCTrack,  "mcreco" );
  io.set_data_to_read( larlite::data::kMCShower, "mcreco" );
  io.set_data_to_read( larlite::data::kMCTruth,  "generator" );  
  io.open();

  larcv::IOManager iolcv( larcv::IOManager::kREAD, "IOManager", larcv::IOManager::kTickBackward );
  iolcv.add_in_file( dlmerged_input_file );
  iolcv.specify_data_read( larcv::kProductImage2D, "wire" );
  iolcv.reverse_all_products();
  iolcv.initialize();

  int nentries = io.get_entries();  

  larutil::SpaceChargeMicroBooNE sce;
  const std::vector<Double_t> orthy = larutil::Geometry::GetME()->GetOrthVectorsY();
  const std::vector<Double_t> orthz = larutil::Geometry::GetME()->GetOrthVectorsZ();
 
  TFile* out = new TFile(output_file.c_str(),"new");
  TTree* ana = new TTree("anadedx", "Analysis of dQdx using truth tracks");

  int pid;
  float res;
  float rad;
  float s_trunk;
  float lm;  
  float pixval;
  float dqdx;
  float pixval_med;
  float dqdx_med;

  ana->Branch("pid",&pid,"pid/I");
  ana->Branch("res",&res,"res/F");
  ana->Branch("rad",&rad,"rad/F");
  ana->Branch("s_trunk",&s_trunk,"s_trunk/F");
  ana->Branch("lm",&lm,"lm/F");      
  ana->Branch("pixval",&pixval,"pixval/F");
  ana->Branch("dqdx",&dqdx,"dqdx/F");
  ana->Branch("pixval_med",&pixval_med,"pixval_med/F");
  ana->Branch("dqdx_med",&dqdx_med,"dqdx_med/F");

  TTree* shrprof = new TTree("showerprofile","Shower profiles");
  std::vector< std::vector<float> > shrprof_start_v;
  std::vector< std::vector<float> > shrprof_end_v;
  std::vector< int > shrprof_pid_v;
  std::vector< float > shrprof_avedqdx_v;
  shrprof->Branch("start_v", &shrprof_start_v );
  shrprof->Branch("end_v", &shrprof_end_v );
  shrprof->Branch("pid_v", &shrprof_pid_v );
  shrprof->Branch("ave_dedq_v",&shrprof_avedqdx_v);

  std::cout << "NUM ENTRIES: " << nentries << std::endl;
  for (int ientry=0; ientry<nentries; ientry++ ) {

    std::cout << "===[ ENTRY " << ientry << " ]===" << std::endl;
    shrprof_start_v.clear();
    shrprof_end_v.clear();
    shrprof_pid_v.clear();
    shrprof_avedqdx_v.clear();
    
    io.go_to(ientry);
    iolcv.read_entry(ientry);

    // Get wire image
    larcv::EventImage2D* ev_adc = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, "wire" );
    const std::vector<larcv::Image2D>& adc_v = ev_adc->as_vector();
    
    // Get truth showers
    std::vector< const larlite::mcshower* > primary_electrons_v;
    std::vector< const larlite::mcshower* > primary_gammas_v;

    larlite::event_mcshower* ev_mcshower = (larlite::event_mcshower*)io.get_data( larlite::data::kMCShower, "mcreco" );
    for (auto const& shr : *ev_mcshower ) {
      if ( shr.Origin()==larlite::simb::kBeamNeutrino ) {

        //if ( shr.TrackID()==shr.MotherTrackID() ) {

        std::vector<float> start(3,0);
        for (int p=0; p<3; p++) {
          start[p] = shr.DetProfile().Position()[p];
        }
        int boundary = 0;
        float dwall = ublarcvapp::dwall( start, boundary );          
        int pid = abs(shr.PdgCode());
        if ( dwall>5.0 ) {
          if ( pid==22 )
            primary_gammas_v.push_back( &shr );
          else if ( pid==11 ) 
            primary_electrons_v.push_back( &shr );
        }
          
        //}//end of if primary        
      }//end of it origin is neutrino
    }

    std::cout << "===== SHOWERS =====" << std::endl;    
    for ( auto const& pshr : primary_gammas_v ) {
      std::cout << "Shower[gamma]: start=(" << pshr->DetProfile().X() << "," << pshr->DetProfile().Y() << "," << pshr->DetProfile().Z() << ")"
                << " E=" << pshr->End().E() << " MeV"
                << std::endl;
    }
    for ( auto const& pshr : primary_electrons_v ) {
      std::cout << "Shower[electron]: start=(" << pshr->DetProfile().X() << "," << pshr->DetProfile().Y() << "," << pshr->DetProfile().Z() << ")"
                << " E=" << pshr->End().E() << " MeV"
                << std::endl;
    }
    std::cout << "--------------------------" << std::endl;

    // GET LARMATCH POINTS
    larlite::event_larflow3dhit* ev_larmatch =
      (larlite::event_larflow3dhit*)io.get_data(larlite::data::kLArFlow3DHit, "larmatch" );
    

    std::vector< std::vector< const larlite::mcshower* >* > shower_list_v;
    shower_list_v.push_back( &primary_gammas_v );
    shower_list_v.push_back( &primary_electrons_v );
    std::vector<int> shower_list_pid_v = { 22, 11 };

    std::vector< float > hit_rad_v( ev_larmatch->size(), -1.0 ); /// this vector holds

    struct TrackPt_t {
      int hitidx;
      int pid;
      float s;
      float res;
      float r;
      float q;
      float dqdx;
      float q_med;
      float dqdx_med;
      float lm;
      std::vector<float> pt;
      bool operator<( const TrackPt_t& rhs ) const
      {
        if ( s>rhs.s) return true;
        return false;
      };
    };

    struct PtQ_t {
      float q;
      float dqdx;
      bool operator<( const PtQ_t& rhs ) const
      {
        if ( q<rhs.q ) return true;
        return false;
      };
    };
    
    typedef std::vector<TrackPt_t> TrackPtList_t;
    std::vector< TrackPtList_t > trackpt_list_v;
    
    int itracklist = -1;
    for ( auto const& ptracklist : shower_list_v ) {
      itracklist++;
      pid = shower_list_pid_v[itracklist];
    
      // loop over showers, get larmatch points within some distance along trajectory of detprofile.
      // assign charge to larmatch points, assign residual range distance, plot
      for ( auto const& ptrack : *ptracklist ) {

        // index of hits filled
        std::vector<int> hit_index_v;
        hit_index_v.reserve( ev_larmatch->size() );

        // line segment defined by detprofile start and direction
        std::vector< std::vector<float> > detpath;
        
        // get start of profile, get ball around these points
        auto const& step = ptrack->DetProfile();

        std::vector<float> fpt = { step.Position().T(),
                                   step.Position().X(),
                                   step.Position().Y(),
                                   step.Position().Z() };
          
        float tick = ublarcvapp::mctools::CrossingPointsAnaMethods::getTick( fpt, 4050.0, nullptr );//&sce );
        std::vector<double> offsets = sce.GetPosOffsets( fpt[1], fpt[2], fpt[3] );
          
        std::vector<float> pathpt(3,0);
        pathpt[0] = (tick-3200.0)*larutil::LArProperties::GetME()->DriftVelocity()*0.5;
        pathpt[1] = fpt[2];// + offsets[1];
        pathpt[2] = fpt[3];// + offsets[2];
                         
        // narrow down the search set
        std::vector<int> search_index_v;
        std::vector< std::vector<float> > point_v;
        std::vector< std::vector<int> > imgcoord_v;
        
        for (int ihit=0; ihit<ev_larmatch->size(); ihit++) {
          
          auto const& hit = (*ev_larmatch)[ihit];
          float d = 0.;
          for (int i=0; i<3; i++)
            d += ( hit[i]-pathpt[i] )*( hit[i]-pathpt[i] );

          if ( d<10.0*10.0 ) {
            search_index_v.push_back( ihit );
            std::vector<float> pt = { hit[0], hit[1], hit[2] };
            point_v.push_back( pt );
            std::vector<int> imgcoord = { hit.targetwire[0], hit.targetwire[1], hit.targetwire[2], hit.tick };
            imgcoord_v.push_back( imgcoord );
          }
        }

        std::cout << "number of hits inside ball around detprofile start: " << point_v.size() << std::endl;
        
        // now collect hits along trunk, i.e. the detprofile direction

        // need to make the next point
        std::vector<float> truedir(3,0);
        const float truelen = 10.0;
        float mag = 0.;
        for (int i=0; i<3; i++) {
          truedir[i] = step.Momentum().Vect()[i];
          mag += truedir[i]*truedir[i];
        }
        if ( mag>0 ) {
          mag = sqrt(mag);
          for (int i=0; i<3; i++)
            truedir[i] /= mag;
        }
        std::vector<float> profend(4,0);
        profend[0] = fpt[0];
        for (int i=0; i<3; i++ )
          profend[i+1] = fpt[i+1] + truelen*truedir[i];        
        std::cout << "truedir: (" << truedir[0] << "," << truedir[1] << "," << truedir[2] << ")" << std::endl;
        std::cout << "define profend: (" << profend[0] << "," << profend[1] << "," << profend[2] << "," << profend[3] << ")" << std::endl;
        std::vector<double> end_offsets = sce.GetPosOffsets( profend[1], profend[2], profend[3] );
        float tickend = ublarcvapp::mctools::CrossingPointsAnaMethods::getTick( profend, 4050.0, nullptr ); //&sce );
        std::vector<float> profendpt(3,0);        
        profendpt[0] = (tickend-3200.0)*larutil::LArProperties::GetME()->DriftVelocity()*0.5;
        profendpt[1] = profend[2];//+end_offsets[1];
        profendpt[2] = profend[3];//+end_offsets[2];

        std::cout << "detprofile defined: "
                  << "(" << pathpt[0] << "," << pathpt[1] << "," << pathpt[2] << ") -> "
                  << "(" << profendpt[0] << "," << profendpt[1] << "," << profendpt[2] << ")"
                  << std::endl;
        
        detpath.push_back( pathpt );
        detpath.push_back( profendpt );

        shrprof_start_v.push_back( pathpt );
        shrprof_end_v.push_back( profendpt );
        shrprof_pid_v.push_back( pid );
        
        float current_len = 0.;
        TrackPtList_t trackpt_v;

        float ave_dqdx = 0;
        int npt_ave_dqdx = 0;
        
        for ( int istep=0; istep<(int)detpath.size()-1; istep++ ) {
          std::vector<float>& start = detpath[istep];
          std::vector<float>& end   = detpath[istep+1];
          std::vector<float> dir(3,0);
          float len = 0.;
          for (int dim=0; dim<3; dim++) {
            dir[dim] += end[dim]-start[dim];
            len += dir[dim]*dir[dim];
          }
          len = sqrt(len);
          
          for (int i=0; i<3; i++ ) {
            dir[i] /= len;
          }
          
          for (int ii=0; ii<(int)point_v.size(); ii++) {
            auto const& pt = point_v[ii];
            auto const& imgcoord = imgcoord_v[ii];
            float r = larflow::reco::pointLineDistance3f( start, end, pt );
            float s = larflow::reco::pointRayProjection3f( start, dir, pt );
            //std::cout << "  point: r=" << r << " s=" << s << std::endl;
            
            if ( r>5.0 || s<-1.0 || s>len+1.0 ) {
              continue;
            }

            // on segment
            TrackPt_t trkpt;
            trkpt.pt     = pt;
            trkpt.hitidx = search_index_v[ii];
            trkpt.pid = pid;
            trkpt.r = r;
            trkpt.s = s+current_len;
            trkpt.q = 0.;            
            trkpt.dqdx = 0.;
            trkpt.q_med = 0.;
            trkpt.dqdx_med = 0.;
            trkpt.lm = ev_larmatch->at(trkpt.hitidx).track_score;

            // get the median charge inside the image
            int row = adc_v.front().meta().row( imgcoord[3] );

            std::vector< PtQ_t > pixq_v(3);

            for ( int p=0; p<3; p++) {

              float pixsum = 0.;
              int npix = 0;
              for (int dr=-2; dr<=2; dr++ ) {
                int r = row+dr;
                if ( r<0 || r>=(int)adc_v.front().meta().rows() )
                  continue;
                pixsum += adc_v[p].pixel( r, imgcoord[p] );
                npix++;
              }
              if ( npix>0 )
                pixq_v[p].q = pixsum/float(npix);
              else
                pixq_v[p].q = 0;
              
              float dcos_yz = fabs(truedir[1]*orthy[p] + truedir[2]*orthz[p]);
              float dcos_x  = fabs(truedir[0]);
              float dx = 3.0;
              if ( dcos_yz>0.785 )
                dx = 3.0/dcos_yz;
              else
                dx = 3.0/dcos_x;
              pixq_v[p].dqdx = pixsum/dx;
            }
            // y-plane only
            trkpt.q = pixq_v[2].q;
            trkpt.dqdx = pixq_v[2].dqdx;

            // median value
            std::sort( pixq_v.begin(), pixq_v.end() );
            trkpt.q_med    = pixq_v[1].q;
            trkpt.dqdx_med = pixq_v[1].dqdx;

            if ( trkpt.dqdx_med>0 ) {
              ave_dqdx += trkpt.dqdx_med;
              npt_ave_dqdx++;
            }
            
            if ( hit_rad_v[trkpt.hitidx]<0 || trkpt.r<hit_rad_v[trkpt.hitidx] )
              hit_rad_v[trkpt.hitidx] = trkpt.r;
            
            trackpt_v.push_back( trkpt );
          }//end of point loop
          current_len += len;
        }//end of loop over detpath steps (only executes once)
        
        std::cout << "Number of hits assigned to track: " << trackpt_v.size() << std::endl;
        std::cout << "Total length of track: " << current_len << " cm" << std::endl;
        std::sort( trackpt_v.begin(), trackpt_v.end() );

        if ( npt_ave_dqdx>0 )
          ave_dqdx /= float(npt_ave_dqdx);
        shrprof_avedqdx_v.push_back( ave_dqdx );
        
        for ( auto& trkpt : trackpt_v ) {
          trkpt.res = current_len - trkpt.s;
        }
        
        trackpt_list_v.emplace_back( std::move(trackpt_v) );
        
      }//end of loop over tracks
    }//end of track list
    
    // fill track points
    int npt_notclose = 0;
    for ( auto& trackpt_v : trackpt_list_v ) {
      for (auto& trkpt : trackpt_v ) {

        // store only if hit is on closest truth track
        if ( trkpt.r<=hit_rad_v[trkpt.hitidx] ) {
          pid = trkpt.pid;
          res = trkpt.res;
          pixval = trkpt.q;
          dqdx = trkpt.dqdx;
          rad = trkpt.r;
          lm = trkpt.lm;
          pixval_med = trkpt.q_med;
          dqdx_med   = trkpt.dqdx_med;
          s_trunk = trkpt.s;
          ana->Fill();
        }
        else {
          npt_notclose++;
        }
      }
    }
    std::cout << "Number of hits in the event not the closest to track: " << npt_notclose << std::endl;
    shrprof->Fill();
  }
  
  out->Write();
  
  io.close();
  iolcv.finalize();
  
  return 0;
}
