#include <iostream>
#include <vector>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TSpline.h"

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
#include "larflow/Reco/TrackTruthRecoAna.h"

int main( int nargs, char** argv )
{

  std::string dlmerged_1m1p_file  = argv[1];
  std::string larmatch_input_file = argv[2];
  std::string kpsreco_ana_file    = argv[3];
  std::string output_file         = argv[4];

  std::string splinefile = argv[5];

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
  // io.add_in_filename( larmatch_input_file );
  // io.set_data_to_read( larlite::data::kLArFlow3DHit, "larmatch" );
  // io.open();

  // filtered file
  larcv::IOManager iolcv( larcv::IOManager::kREAD, "IOManager", larcv::IOManager::kTickBackward );
  iolcv.add_in_file( dlmerged_1m1p_file );
  iolcv.specify_data_read( larcv::kProductImage2D, "wire" );
  iolcv.reverse_all_products();
  iolcv.initialize();

  // full event file
  TFile fkpsreco( kpsreco_ana_file.c_str(), "open" );
  TTree* kpsreco = (TTree*)fkpsreco.Get("KPSRecoManagerTree"); // event-indexed tree
  int kpsreco_run;
  int kpsreco_subrun;
  int kpsreco_event;  
  std::vector<larflow::reco::NuVertexCandidate>*  nufitted_v = nullptr;
  std::vector<larflow::reco::VertexTrackTruthRecoInfo>* truthmatch_vtxinfo_v = nullptr;  
  kpsreco->SetBranchAddress( "run",    &kpsreco_run );
  kpsreco->SetBranchAddress( "subrun", &kpsreco_subrun );
  kpsreco->SetBranchAddress( "event",  &kpsreco_event );    
  kpsreco->SetBranchAddress( "nufitted_v", &nufitted_v );
  kpsreco->SetBranchAddress( "track_truthreco_vtxinfo_v", &truthmatch_vtxinfo_v );

  /**
   * DLMERGED 1M1P FILE
   */
  TFile fdlana( dlmerged_1m1p_file.c_str(), "open" );

  // Vertex-indexed tree with analysis variables
  TTree* fvv = (TTree*)fdlana.Get("dlana/FinalVertexVariables"); // vertex-indexed tree
  float bdtscore_1m1p;
  int fvv_run;
  int fvv_subrun;
  int fvv_event;
  int fvv_vtxid;
  float xreco;
  float yreco;
  float zreco;  
  fvv->SetBranchAddress("run",    &fvv_run);
  fvv->SetBranchAddress("subrun", &fvv_subrun);
  fvv->SetBranchAddress("event",  &fvv_event);
  fvv->SetBranchAddress("vtxid",  &fvv_vtxid);  
  fvv->SetBranchAddress("Xreco",  &xreco );
  fvv->SetBranchAddress("Yreco",  &yreco );
  fvv->SetBranchAddress("Zreco",  &zreco );
  fvv->SetBranchAddress("BDTscore_1mu1p_nu", &bdtscore_1m1p);

  // event-indexed trees with run,subrun,event
  int ubdlana_run;
  int ubdlana_subrun;
  int ubdlana_event;
  TTree* ubdlana_evtree = (TTree*)fdlana.Get("dlana/ubdlana_id_tree"); // vertex-indexed tree
  ubdlana_evtree->SetBranchAddress("run",    &ubdlana_run);
  ubdlana_evtree->SetBranchAddress("subrun", &ubdlana_subrun);
  ubdlana_evtree->SetBranchAddress("event",  &ubdlana_event);  

  // spline file
  float q2adc = 93.0/2.2;

  TFile* splinefile_rootfile = new TFile( splinefile.c_str(), "open" );
  TSpline3* sMuonRange2dEdx = (TSpline3*)splinefile_rootfile->Get("sMuonRange2dEdx");
  TSpline3* sProtonRange2dEdx = (TSpline3*)splinefile_rootfile->Get("sProtonRange2dEdx");  

  // WE ARE GOING TO LOOP OVER THE FILTERED ENTRIES
  int nentries = ubdlana_evtree->GetEntries();

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
  float llpid_pt;
  float truth_mse_track;
  int truth_pid_track;

  ana->Branch("pid",&pid,"pid/I");
  ana->Branch("res",&res,"res/F");
  ana->Branch("rad",&rad,"rad/F");
  ana->Branch("lm",&lm,"lm/F");      
  ana->Branch("pixval",&pixval,"pixval/F");
  ana->Branch("dqdx",&dqdx,"dqdx/F");
  ana->Branch("pixval_med",&pixval_med,"pixval_med/F");
  ana->Branch("dqdx_med",&dqdx_med,"dqdx_med/F");
  ana->Branch("llpid_pt",&llpid_pt,"llpid_pt");
  ana->Branch("truth_mse",&truth_mse_track,"truth_mse/F");
  ana->Branch("truth_pid",&truth_pid_track,"truth_pid/I");  

  TTree* llana = new TTree("llana","log-likelihood score per track");
  int pid_track;
  float llpid_track;
  float len_track;
  float ddlvertex;
  llana->Branch("pid",&pid_track,"pid/I");
  llana->Branch("llpid",&llpid_track,"llpid/F");
  llana->Branch("len",&len_track,"len/F");
  llana->Branch("ddlvertex",&ddlvertex,"ddlvertex/F");
  llana->Branch("truth_mse",&truth_mse_track,"truth_mse/F");
  llana->Branch("truth_pid",&truth_pid_track,"truth_pid/I");  

  std::cout << "NUM ENTRIES: " << nentries << std::endl;
  int current_fvv_entry = 0;
  int current_fvv_run = 0;
  int current_fvv_subrun = 0;
  int current_fvv_event = 0;
  float min_bdt_score = 1e9;
  int min_bdt_entry = -1;

  // FVV LOOP WILL BE WEIRD
  for (int ientry=0; ientry<nentries; ientry++ ) {

    std::cout << "===[ ENTRY " << ientry << " ]===" << std::endl;
    ubdlana_evtree->GetEntry(ientry);
    iolcv.read_entry(ientry);
    
    larcv::EventImage2D* ev_adc = nullptr;

    ev_adc = (larcv::EventImage2D*)iolcv.get_data( larcv::kProductImage2D, "wire" );
    
    // get the run, subrun, event number
    int run    = ev_adc->run();
    int subrun = ev_adc->subrun();
    int event  = ev_adc->event();
    
    // SEARCH FOR THE (RUN,SUBRUN,ENTRY) in the larcv,larlite,and kpsreco trees
    bool found_rse = false;
    for (int iunfilter=0; iunfilter<(int)kpsreco->GetEntries(); iunfilter++) {
    
      //io.go_to(iunfilter);
      kpsreco->GetEntry(iunfilter);
      
      if ( ubdlana_run==kpsreco_run && ubdlana_subrun==kpsreco_subrun && ubdlana_event==kpsreco_event
           && ubdlana_run==run && ubdlana_subrun==subrun && ubdlana_event==event ) {
        found_rse = true;
        break;
      }
    }

    if ( !found_rse ) {
      throw std::runtime_error("COULD NOT FIND RSE IN FILTERED FILE IN THE FULL FILES");
    }

    // Get wire image
    const std::vector<larcv::Image2D>& adc_v = ev_adc->as_vector();

    // find the min-bdt-score vertex in the dlana. dumb complete reloop
    int matching_fvv_entry = -1;
    float min_bdt_score = 2.;
    for (int fvv_entry=0; fvv_entry<(int)fvv->GetEntries(); fvv_entry++) {
      fvv->GetEntry(fvv_entry);
      if ( fvv_run==run && fvv_subrun==subrun && fvv_event==event ) {
        if ( min_bdt_score>bdtscore_1m1p ) {
          min_bdt_score = bdtscore_1m1p;
          matching_fvv_entry = fvv_entry;
        }
      }
    }
    fvv->GetEntry( matching_fvv_entry );
    
    // now the (x,y,z) coordinates of the DL vertex is loaded, in principle
    if ( matching_fvv_entry<0 ) {
      std::stringstream ss;
      ss << "Did not find DL vertex in dlana tree" << std::endl
         << "  looking for (" << run << "," << subrun << "," << event << ")" << std::endl;
      throw std::runtime_error(ss.str());
    }
    else {
      std::cout << "BEST DL vertex[" << fvv_vtxid << "]: (" << xreco << "," << yreco << "," << zreco << ") bdt-score=" << min_bdt_score << std::endl;
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
    auto const& truthmatch_vtxinfo = truthmatch_vtxinfo_v->at(closest_vertex_entry);

    struct TrackPt_t {
      int hitidx;
      int pid;
      int truthmatch_pid;
      float s;
      float res;
      float r;
      float q;
      float dqdx;
      float q_med;
      float dqdx_med;
      float lm;
      float ll;
      float llw;
      float truthmatch_mse;
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
    std::vector<float> track_ll_v( nuvertex.track_v.size(), 0 );
    typedef std::vector<TrackPt_t> TrackPtList_t;
    std::vector< TrackPtList_t > trackpt_list_v;

    std::cout << "[vertex] ntrack=" << nuvertex.track_v.size() << " nhitcluster=" << nuvertex.track_hitcluster_v.size() << std::endl;

    for (int itrack=0; itrack<nuvertex.track_v.size(); itrack++) {
      
      const larlite::track& lltrack = nuvertex.track_v.at(itrack);
      const larlite::larflowcluster& lfcluster = nuvertex.track_hitcluster_v.at(itrack);
      const larflow::reco::TrackTruthRecoInfo& truthmatch_trackinfo = truthmatch_vtxinfo.trackinfo_v.at(itrack);

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
          trkpt.truthmatch_pid = truthmatch_trackinfo.matched_true_pid;
          trkpt.truthmatch_mse = truthmatch_trackinfo.matched_mse;

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

      // calculate residual range
      // calculate likelihood
      float totw = 0.;
      float totll = 0.;
      for ( auto& trkpt : trackpt_v ) {
        trkpt.res = current_len - trkpt.s;

        float mu_dedx = sMuonRange2dEdx->Eval(trkpt.res);
        float mu_dedx_birks = q2adc*mu_dedx/(1+mu_dedx*0.0486/0.273/1.38);
        float p_dedx = sProtonRange2dEdx->Eval(trkpt.res);
        float p_dedx_birks = q2adc*p_dedx/(1+p_dedx*0.0486/0.273/1.38);
        
        float dmu = trkpt.dqdx_med-mu_dedx_birks;
        float dp  = trkpt.dqdx_med-p_dedx_birks;

        float llpt = -0.5*dmu*dmu/100.0 + 0.5*dp*dp/100.0;
        float w_dedx = (mu_dedx_birks-p_dedx_birks)*(mu_dedx_birks-p_dedx_birks);
        trkpt.ll = llpt;
        trkpt.llw = w_dedx;
	if ( trkpt.dqdx_med>10.0 ) {
	  totll += llpt*w_dedx;
	  totw  += w_dedx;
	}
      }
      if ( totw>0 )
        totll /= totw;

      track_len_v[itrack] = current_len;
      track_ll_v[itrack] = totll;
      
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

    for (int itrack=0; itrack<(int)nuvertex.track_v.size(); itrack++) {
      if ( itrack==longest_len_id )
        pid_track = 13;
      else if ( itrack==shortest_len_id && shortest_len_id!=longest_len_id )
        pid_track = 2212;
      else if ( shortest_len_id!=longest_len_id )
        pid_track = 111;
      llpid_track = track_ll_v[itrack];
      len_track   = track_len_v[itrack];
      ddlvertex   = dist_to_dlvertex;

      truth_mse_track = truthmatch_vtxinfo.trackinfo_v.at(itrack).matched_mse;
      truth_pid_track = truthmatch_vtxinfo.trackinfo_v.at(itrack).matched_true_pid;
      
      llana->Fill();
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
        truth_pid_track = trkpt.truthmatch_pid;
        truth_mse_track = trkpt.truthmatch_mse;
        ana->Fill();
      }
    }
  }
  
  out->Write();
  
  //io.close();
  iolcv.finalize();
  
  return 0;
}
