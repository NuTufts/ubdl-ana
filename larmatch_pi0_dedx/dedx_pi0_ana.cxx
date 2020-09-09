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

  // per space-point tree
  TTree* ana = new TTree("recodedx", "Analysis of dQdx using reco tracks");

  float dist2start;  // distance of spacepoint to shower start
  float rad;         // distance from radius
  float lm;          // larmatch score
  float pixval;      // pixval sum on Y plane
  float dqdx;        // dqdx on Y plane
  float pixval_med;  // median plane pixval 
  float dqdx_med;    // median plane dqdx
  float llpid_pt;    // -log-likelihood ratio

  // shower level variables to cut on
  float dist2vertex_shower; // distance of shower to vertex  
  float llpid_shower;
  float qtot_shower;
  float ddlvertex;
  
  ana->Branch("dist2start", &dist2start, "dist2start/F");
  ana->Branch("rad",&rad,"rad/F");
  ana->Branch("lm",&lm,"lm/F");      
  ana->Branch("pixval",&pixval,"pixval/F");
  ana->Branch("dqdx",&dqdx,"dqdx/F");
  ana->Branch("pixval_med",&pixval_med,"pixval_med/F");
  ana->Branch("dqdx_med",&dqdx_med,"dqdx_med/F");
  ana->Branch("llpid_pt",&llpid_pt,"llpid_pt");
  ana->Branch("dist2vertex_shower",&dist2vertex_shower,"dist2vertex_shower/F");
  ana->Branch("qtot_shower",&qtot_shower,"qtot_shower/F");  

  // per shower tree
  TTree* llana = new TTree("llana","log-likelihood score per track");
  ana->Branch("dist2vertex",&dist2vertex_shower,"dist2vertex/F");  
  llana->Branch("llpid",&llpid_shower,"llpid/F");
  llana->Branch("qtot",&qtot_shower,"qtot/F");
  llana->Branch("ddlvertex",&ddlvertex,"ddlvertex/F");  

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

    // get the vertex!
    auto const& nuvertex = nufitted_v->at(closest_vertex_entry);

    struct ShowerPt_t {
      int hitidx;
      int pid;
      float s;
      float r;
      float q;
      float dqdx;
      float q_med;
      float dqdx_med;
      float lm;
      float ll;
      float llw;
      std::vector<float> pt;
      bool operator<( const ShowerPt_t& rhs ) const
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

    std::vector<float> shower_d2vtx_v( nuvertex.shower_v.size(), 0);
    std::vector<float> shower_ll_v(    nuvertex.shower_v.size(), 0 );
    std::vector<float> shower_qtot_v(  nuvertex.shower_v.size(), 0 );    
    typedef std::vector<ShowerPt_t> ShowerPtList_t;
    std::vector< ShowerPtList_t > showerpt_list_v;

    std::cout << "[vertex] nshowers=" << nuvertex.shower_v.size() << " ntrunks=" << nuvertex.shower_trunk_v.size() << std::endl;

    for (int ishower=0; ishower<nuvertex.shower_v.size(); ishower++) {
      
      const larlite::track& lltrunk = nuvertex.shower_trunk_v.at(ishower);
      const larlite::larflowcluster& lfcluster = nuvertex.shower_v.at(ishower);

      std::vector< float > hit_rad_v( lfcluster.size(), -1.0 ); /// this vector holds
                
      std::vector< std::vector<float> > detpath; // will hold start pt then end point
      TVector3 pt_0 = lltrunk.LocationAtPoint(0);
      TVector3 pt_1 = lltrunk.LocationAtPoint(1);

      float dist0 = 0.;
      float dist1 = 0.;
      
      for (int i=0; i<3; i++) {
        dist0 += (pt_0[i]-nuvertex.pos[i])*(pt_0[i]-nuvertex.pos[i]);
        dist1 += (pt_1[i]-nuvertex.pos[i])*(pt_1[i]-nuvertex.pos[i]);
      }
      dist0 = sqrt(dist0);
      dist1 = sqrt(dist1);
      if ( dist0<dist1 ) {
        detpath.push_back( std::vector<float>( {(float)pt_0[0], (float)pt_0[1], (float)pt_0[2] } ) );
        detpath.push_back( std::vector<float>( {(float)pt_1[0], (float)pt_1[1], (float)pt_1[2] } ) );
        shower_d2vtx_v[ishower] = dist0;
      }
      else {
        detpath.push_back( std::vector<float>( {(float)pt_1[0], (float)pt_1[1], (float)pt_1[2] } ) );
        detpath.push_back( std::vector<float>( {(float)pt_0[0], (float)pt_0[1], (float)pt_0[2] } ) );
        shower_d2vtx_v[ishower] = dist1;
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
      
      // Make ShowerPtList 
      float current_len = 0.;
      ShowerPtList_t showerpt_v;
        
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
          ShowerPt_t shwrpt;
          shwrpt.pt     = pt;
          shwrpt.hitidx = search_index_v[ii];
          shwrpt.pid = ishower;
          shwrpt.r = r;
          shwrpt.s = s+current_len;
          shwrpt.q = 0.;            
          shwrpt.dqdx = 0.;
          shwrpt.q_med = 0.;
          shwrpt.dqdx_med = 0.;
          shwrpt.lm = lfcluster.at(shwrpt.hitidx).track_score;

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
          shwrpt.q = pixq_v[2].q;
          shwrpt.dqdx = pixq_v[2].dqdx;
          
          // median value
          std::sort( pixq_v.begin(), pixq_v.end() );
          shwrpt.q_med    = pixq_v[1].q;
          shwrpt.dqdx_med = pixq_v[1].dqdx;
          
          if ( hit_rad_v[shwrpt.hitidx]<0 || shwrpt.r<hit_rad_v[shwrpt.hitidx] )
            hit_rad_v[shwrpt.hitidx] = shwrpt.r;
          
          showerpt_v.push_back( shwrpt );

        }//end of point loop
        
        current_len += len;
      }//end of loop over detpath steps
      
      std::cout << "Number of hits assigned to track: " << showerpt_v.size() << std::endl;
      std::cout << "Total length of track: " << current_len << " cm" << std::endl;
      std::sort( showerpt_v.begin(), showerpt_v.end() );

      // calculate residual range
      // calculate likelihood
      float totw = 0.;
      float totll = 0.;
      float qsum_med = 0.;      
      for ( auto& shwrpt : showerpt_v ) {

        float dgamma     = shwrpt.dqdx_med-140.0;
        float delectron  = shwrpt.dqdx_med-75.0;

        float llpt = -0.5*dgamma*dgamma/(25.0*25.0) + 0.5*delectron*delectron/(15.0*15.0);
        float w_dedx = 1.0;
        shwrpt.ll = llpt;
        shwrpt.llw = w_dedx;
	if ( shwrpt.dqdx_med>10.0 && shwrpt.s>0 && shwrpt.s<3.0 ) {
	  totll += llpt*w_dedx;
	  totw  += w_dedx;
	}
        qsum_med += shwrpt.q_med;
      }
      if ( totw>0 )
        totll /= totw;

      shower_ll_v[ishower] = totll;
      shower_qtot_v[ishower] = qsum_med;
      
      showerpt_list_v.emplace_back( std::move(showerpt_v) );

        
    }//end of track list

    // Fill Shower Info
    for (int ishower=0; ishower<(int)nuvertex.shower_v.size(); ishower++) {

      dist2vertex_shower = shower_d2vtx_v[ishower];
      llpid_shower = shower_ll_v[ishower];
      qtot_shower  = shower_qtot_v[ishower];
      ddlvertex   = dist_to_dlvertex;
      llana->Fill();

      // fill track points
      auto& showerpt_v = showerpt_list_v[ishower];
      
      for (auto& shwrpt : showerpt_v ) {
        
        dist2start = shwrpt.s;
        rad = shwrpt.r;
        lm = shwrpt.lm;
        pixval = shwrpt.q;
        dqdx = shwrpt.dqdx;
        pixval_med = shwrpt.q_med;
        dqdx_med   = shwrpt.dqdx_med;
        llpid_pt   = shwrpt.ll;
        ana->Fill();
      }
    }
  }
  
  out->Write();
  
  //io.close();
  iolcv.finalize();
  
  return 0;
}
