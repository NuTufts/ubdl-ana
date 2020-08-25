#include <iostream>
#include <vector>
#include <sstream>

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

#include "ublarcvapp/MCTools/crossingPointsAnaMethods.h"
#include "ublarcvapp/ubdllee/dwall.h"

#include "larflow/Reco/geofuncs.h"
#include "larflow/Reco/NuVertexCandidate.h"

int main( int nargs, char** argv )
{

  std::string dlmerged_input_file = argv[1];
  std::string kpsreco_ana_file = argv[2];
  std::string output_file = argv[3];

  // what inputs we need
  // --------------------
  // ubdlana/FinalVertexVariables (dlfilter/dlmerged/dlana) file:
  //  - maximum BDT for the (run,subrun,event)
  //  - vtx of the one with maximum bdt. will use this to pick best NuVertexCandidate
  // fitted NuVertexCandidate in KPSRecoManagerTree tree (KPS reco ana file)
  /* Contents of the NuVertexClass
*Br   30 :nufitted_v : Int_t nufitted_v_                                     *
*Entries :        1 : Total  Size=       6426 bytes  File Size  =        107 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br   31 :nufitted_v.keypoint_producer :                                     *
*         | string keypoint_producer[nufitted_v_]                            *
*Entries :        1 : Total  Size=        821 bytes  File Size  =        152 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.01     *
*............................................................................*
*Br   32 :nufitted_v.keypoint_index : Int_t keypoint_index[nufitted_v_]      *
*Entries :        1 : Total  Size=        785 bytes  File Size  =        130 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br   33 :nufitted_v.pos : vector<float> pos[nufitted_v_]                    *
*Entries :        1 : Total  Size=        772 bytes  File Size  =        161 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br   34 :nufitted_v.cluster_v : vector<larflow::reco::NuVertexCandidate:    *
*         | :VtxCluster_t> cluster_v[nufitted_v_]                            *
*Entries :        1 : Total  Size=       1176 bytes  File Size  =        376 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.44     *
*............................................................................*
*Br   35 :nufitted_v.score : Float_t score[nufitted_v_]                      *
*Entries :        1 : Total  Size=        740 bytes  File Size  =        121 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br   36 :nufitted_v.track_v : vector<larlite::track> track_v[nufitted_v_]   *
*Entries :        1 : Total  Size=       2234 bytes  File Size  =        548 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   2.93     *
*............................................................................*
*Br   37 :nufitted_v.shower_v : vector<larlite:                              *
*         | :larflowcluster> shower_v[nufitted_v_]                           *
*Entries :        1 : Total  Size=      16088 bytes  File Size  =       1814 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   8.52     *
*............................................................................*
*Br   38 :nufitted_v.shower_trunk_v : vector<larlite:                        *
*         | :track> shower_trunk_v[nufitted_v_]                              *
*Entries :        1 : Total  Size=       1445 bytes  File Size  =        240 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   3.29     *
*............................................................................*
   */
  
  
  // larlite::storage_manager io( larlite::storage_manager::kREAD );
  // io.add_in_filename( dlmerged_input_file );
  // io.set_data_to_read( larlite::data::kLArFlow3DHit, "larmatch" );
  // io.open();

  larcv::IOManager iolcv( larcv::IOManager::kREAD, "IOManager", larcv::IOManager::kTickBackward );
  iolcv.add_in_file( dlmerged_input_file );
  iolcv.specify_data_read( larcv::kProductImage2D, "wire" );
  iolcv.reverse_all_products();
  iolcv.initialize();

  TFile fkpsreco( kpsreco_ana_file.c_str(), "open" );
  TTree* kpsreco = (TTree*)fkpsreco.Get("KPSRecoManagerTree"); // event-indexed tree
  int kpsreco_run;
  int kpsreco_subrun;
  int kpsreco_event;  
  std::vector<larflow::reco::NuVertexCandidate>*  nufitted_v = nullptr;
  kpsreco->SetBranchAddress( "run",    &kpsreco_run );
  kpsreco->SetBranchAddress( "subrun", &kpsreco_subrun );
  kpsreco->SetBranchAddress( "event",  &kpsreco_event );    
  kpsreco->SetBranchAddress( "nufitted_v", &nufitted_v );

  TFile fdlana( dlmerged_input_file.c_str(), "open" );
  TTree* fvv = (TTree*)fdlana.Get("dlana/FinalVertexVariables"); // vertex-indexed tree
  float bdtscore_1m1p;
  int dlana_run;
  int dlana_subrun;
  int dlana_event;
  int dlana_vtxid;
  float xreco;
  float yreco;
  float zreco;  
  fvv->SetBranchAddress("run",    &dlana_run);
  fvv->SetBranchAddress("subrun", &dlana_subrun);
  fvv->SetBranchAddress("event",  &dlana_event);
  fvv->SetBranchAddress("vtxid",  &dlana_vtxid);  
  fvv->SetBranchAddress("Xreco",  &xreco );
  fvv->SetBranchAddress("Yreco",  &yreco );
  fvv->SetBranchAddress("Zreco",  &zreco );
  fvv->SetBranchAddress("BDTscore_1mu1p_nu", &bdtscore_1m1p);  

  int nentries = iolcv.get_n_entries();

  larutil::SpaceChargeMicroBooNE sce;
  const std::vector<Double_t> orthy = larutil::Geometry::GetME()->GetOrthVectorsY();
  const std::vector<Double_t> orthz = larutil::Geometry::GetME()->GetOrthVectorsZ();
 
  TFile* out = new TFile(output_file.c_str(),"new");
  TTree* ana = new TTree("recodedx", "Analysis of dQdx using reco tracks");

  int pid;
  float res;
  float rad;
  float lm;  
  float pixval;
  float dqdx;
  float pixval_med;
  float dqdx_med;

  ana->Branch("pid",&pid,"pid/I");
  ana->Branch("res",&res,"res/F");
  ana->Branch("rad",&rad,"rad/F");
  ana->Branch("lm",&lm,"lm/F");      
  ana->Branch("pixval",&pixval,"pixval/F");
  ana->Branch("dqdx",&dqdx,"dqdx/F");
  ana->Branch("pixval_med",&pixval_med,"pixval_med/F");
  ana->Branch("dqdx_med",&dqdx_med,"dqdx_med/F");

  std::cout << "NUM ENTRIES: " << nentries << std::endl;
  for (int ientry=0; ientry<nentries; ientry++ ) {

    std::cout << "===[ ENTRY " << ientry << " ]===" << std::endl;
    
    //io.go_to(ientry);
    iolcv.read_entry(ientry);
    kpsreco->GetEntry(ientry);

    // Get wire image
    larcv::EventImage2D* ev_adc = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, "wire" );
    const std::vector<larcv::Image2D>& adc_v = ev_adc->as_vector();

    // get the run, subrun, event number
    int run    = ev_adc->run();
    int subrun = ev_adc->subrun();
    int event  = ev_adc->event();

    // find the vertex in the dlana. dumb complete reloop
    int matching_dlana_entry = -1;
    float min_bdt_score = 2.;
    for (int dlana_entry=0; dlana_entry<(int)fvv->GetEntries(); dlana_entry++) {
      fvv->GetEntry(dlana_entry);
      if ( dlana_run==run && dlana_subrun==subrun && dlana_event==event ) {
        if ( min_bdt_score>bdtscore_1m1p ) {
          min_bdt_score = bdtscore_1m1p;
          matching_dlana_entry = dlana_entry;
        }
      }
    }
    fvv->GetEntry( matching_dlana_entry );

    // now the (x,y,z) coordinates of the DL vertex is loaded, in principle
    if ( matching_dlana_entry<0 ) {
      std::stringstream ss;
      ss << "Did not find DL vertex in dlana tree" << std::endl
         << "  looking for (" << run << "," << subrun << "," << event << ")" << std::endl;
      throw std::runtime_error(ss.str());
    }
    else {
      std::cout << "DL vertex[" << dlana_vtxid << "]: (" << xreco << "," << yreco << "," << zreco << ")" << std::endl;
    }

    // Get the nuvertexcandidate closest to the DL vertex
    int closest_vertex_entry = -1;
    float dist_to_dlvertex = 1.0e9;
    for ( int ivtx=0; ivtx<(int)nufitted_v->size(); ivtx++ ) {
      auto const& nuvtx = nufitted_v->at(ivtx);
      
      float dist = 0.;
      dist += (nuvtx.pos[0]-xreco)*(nuvtx.pos[0]-xreco);
      dist += (nuvtx.pos[1]-yreco)*(nuvtx.pos[1]-yreco);
      dist += (nuvtx.pos[2]-zreco)*(nuvtx.pos[2]-zreco);
      dist = sqrt(dist);
      if ( dist < dist_to_dlvertex ) {
        dist_to_dlvertex = dist;
        closest_vertex_entry = ivtx;
      }
    }

    if ( closest_vertex_entry<0 )
      throw std::runtime_error("Did not indentify a closest vertex");

    std::cout << "Distance to closest DL vertex: " << dist_to_dlvertex << " cm. NuVtx Index: " << closest_vertex_entry << std::endl;
    
    auto const& nuvertex = nufitted_v->at(closest_vertex_entry);

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

    std::vector<float> track_len_v( nuvertex.track_v.size(), 0);
    typedef std::vector<TrackPt_t> TrackPtList_t;
    std::vector< TrackPtList_t > trackpt_list_v;

    std::cout << "[vertex] ntrack=" << nuvertex.track_v.size() << " nhitcluster=" << nuvertex.track_hitcluster_v.size() << std::endl;

    for (int itrack=0; itrack<nuvertex.track_v.size(); itrack++) {
      
      const larlite::track& lltrack = nuvertex.track_v.at(itrack);
      const larlite::larflowcluster& lfcluster = nuvertex.track_hitcluster_v.at(itrack);

      std::vector< float > hit_rad_v( lfcluster.size(), -1.0 ); /// this vector holds
                
      std::vector< std::vector<float> > detpath;
        
      for ( int istep=0; istep<(int)lltrack.NumberTrajectoryPoints(); istep++ )  {
        TVector3 pt = lltrack.LocationAtPoint(istep);
        std::vector<float> fpt = { (float)pt[0], (float)pt[1], (float)pt[2] };
        detpath.push_back( fpt );
      }
        
      // collect hit position and img coord
      std::vector<int> search_index_v;
      std::vector< std::vector<float> > point_v;
      std::vector< std::vector<int> > imgcoord_v;
        
      for (int ihit=0; ihit<(int)lfcluster.size(); ihit++) {
        auto const& hit = lfcluster[ihit];
        search_index_v.push_back( ihit );
        std::vector<float> pt = { hit[0], hit[1], hit[2] };
        point_v.push_back( pt );
        std::vector<int> imgcoord = { hit.targetwire[0], hit.targetwire[1], hit.targetwire[2], hit.tick };
        imgcoord_v.push_back( imgcoord );
      }
      
      // Make TrackPtList 
      float current_len = 0.;
      TrackPtList_t trackpt_v;
        
      for ( int istep=0; istep<(int)detpath.size()-1; istep++ ) {
        std::vector<float>& start = detpath[istep];
        std::vector<float>& end   = detpath[istep+1];
        std::vector<float> dir(3,0);
        std::vector<float> truedir(3,0);          
        float len = 0.;
        float truelen = 0.;
        for (int dim=0; dim<3; dim++) {
          dir[dim] += end[dim]-start[dim];
          len += dir[dim]*dir[dim];

          truedir[dim] = end[dim]-start[dim];
          truelen += truedir[dim]*truedir[dim];
        }
        len = sqrt(len);
        truelen = sqrt(truelen);
          
        if (len<=0.1 ) {
          current_len += len;
          continue;
        }
          
        for (int i=0; i<3; i++ ) {
          dir[i] /= len;
          truedir[i] /= truelen;
        }
          
        for (int ii=0; ii<(int)point_v.size(); ii++) {
          auto const& pt = point_v[ii];
          auto const& imgcoord = imgcoord_v[ii];
          float r = larflow::reco::pointLineDistance3f( start, end, pt );
          float s = larflow::reco::pointRayProjection3f( start, dir, pt );
          //std::cout << "  point: r=" << r << " s=" << s << std::endl;
          
          if ( r>5.0 || s<0 || s>len ) {
            continue;
          }

          // on segment
          TrackPt_t trkpt;
          trkpt.pt     = pt;
          trkpt.hitidx = search_index_v[ii];
          trkpt.pid = itrack;
          trkpt.r = r;
          trkpt.s = s+current_len;
          trkpt.q = 0.;            
          trkpt.dqdx = 0.;
          trkpt.q_med = 0.;
          trkpt.dqdx_med = 0.;
          trkpt.lm = lfcluster.at(trkpt.hitidx).track_score;

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
          
          if ( hit_rad_v[trkpt.hitidx]<0 || trkpt.r<hit_rad_v[trkpt.hitidx] )
            hit_rad_v[trkpt.hitidx] = trkpt.r;
          
          trackpt_v.push_back( trkpt );

        }//end of point loop
        
        current_len += len;
      }//end of loop over detpath steps
      
      std::cout << "Number of hits assigned to track: " << trackpt_v.size() << std::endl;
      std::cout << "Total length of track: " << current_len << " cm" << std::endl;
      std::sort( trackpt_v.begin(), trackpt_v.end() );

      for ( auto& trkpt : trackpt_v ) {
        trkpt.res = current_len - trkpt.s;
      }

      track_len_v[itrack] = current_len;
        
      trackpt_list_v.emplace_back( std::move(trackpt_v) );
        
    }//end of track list

    // we have no PID of course, but we can label by longest and shortest track
    std::vector<int> pid_by_len_v(nuvertex.track_v.size(),0);
    float longest_len = 0;
    int longest_len_id = 0;
    int shortest_len_id = 0;
    float shortest_len = 1e9;

    for (int itrack=0; itrack<(int)nuvertex.track_v.size(); itrack++) {
      if ( longest_len<track_len_v[itrack] ) {
        longest_len = track_len_v[itrack];
        longest_len_id = itrack;
      }
      if ( shortest_len>track_len_v[itrack] ) {
	shortest_len = track_len_v[itrack];
	shortest_len_id = itrack;
      }
    }
    
    
    // fill track points
    for ( auto& trackpt_v : trackpt_list_v ) {
      for (auto& trkpt : trackpt_v ) {

        if ( trkpt.pid==longest_len_id )
          pid = 13;
	else if ( trkpt.pid==shortest_len_id && shortest_len_id!=longest_len_id )
          pid = 2212;
	else if ( shortest_len_id!=longest_len_id )
	  pid = 111;
        
        res = trkpt.res;
        pixval = trkpt.q;
        dqdx = trkpt.dqdx;
        rad = trkpt.r;
        lm = trkpt.lm;
        pixval_med = trkpt.q_med;
        dqdx_med   = trkpt.dqdx_med;
        ana->Fill();
      }
    }
  }
  
  out->Write();
  
  //io.close();
  iolcv.finalize();
  
  return 0;
}
