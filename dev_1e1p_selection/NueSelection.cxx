#include "NueSelection.h"

#include "ublarcvapp/ubdllee/dwall.h"

std::vector<bool>
NueSelection::runVertexSelection(  const larflow::reco::NuSelectionVariables& nusel,
                                   const larflow::reco::NuVertexCandidate& nuvtx,
                                   bool cut_fv )
{
  
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

  // selection cuts
  std::vector<bool> vtx_pass( kNumCuts, false );
  vtx_pass[kFV] = cut_fv; // [0]
  vtx_pass[kVertexCand3cm] = nusel.dist2truevtx<3.0; // [1]
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
  
  vtx_pass[kAllCuts]       = true;
  
  // reco variable cuts only
  for ( int i=kMinShowerSize; i<kAllCuts; i++)
    vtx_pass[kAllCuts] = vtx_pass[kAllCuts] && vtx_pass[i];
  
  // bool vtx_seq = true;
  // for (int icut=0; icut<kNumCuts; icut++) {
  //   vtx_seq = vtx_seq && vtx_pass[icut]; // follows sequence
  //   // if still true mark as passing
  //   // or if previously passed, another vertex had passed this stage
  //   event_passes_cut[ icut ] = event_passes_cut[icut] || vtx_seq;
  // }

  return vtx_pass;
  
}


bool NueSelection::runEventSelection( const std::vector< larflow::reco::NuSelectionVariables>& nusel_v,
                                      const std::vector< larflow::reco::NuVertexCandidate>& nuvtx_v,
                                      const float vtx_dwall,
                                      const int ientry )
{

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
  
  std::vector<EventDist2True_t> index_by_dist_v;
  //std::vector<bool> event_passes_cut( kNumCuts, false ); //< indicates if event has a vertex that passes this cut stage
  event_passes_cut.clear();
  event_passes_cut.resize( kNumCuts, false ); //< indicates if event has a vertex that passes this cut stage
  event_vtx_pass.clear();
  
  int nvtx = (int)(nusel_v.size());

  bool cut_fv = (vtx_dwall>10.0);

  for (int ivtx=0; ivtx<nvtx; ivtx++) {

    auto const& nusel = nusel_v[ivtx];
    auto const& nuvtx = nuvtx_v[ivtx];
    
    EventDist2True_t idx( nusel.dist2truevtx, ivtx );

    std::vector<bool> vtx_pass = runVertexSelection( nusel, nuvtx, cut_fv );

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

    event_vtx_pass.push_back( vtx_pass );
  }
  std::sort( index_by_dist_v.begin(), index_by_dist_v.end() );

  // for debug
  if ( DEBUG) {
    std::cout << "[entry results] ------------" <<  std::endl;
    for (int i=0; i<=kAllCuts; i++) {
      std::cout << "  " << selcut_names[i] << ": " << event_passes_cut[i] << std::endl;
    }
    std::cout << "----------------------------" << std::endl;
  }
  
  if ( event_passes_cut[kAllCuts] ) {
    event_vtx_index = index_by_dist_v.front().index;
    event_vtx_dist  = index_by_dist_v.front().dist;
  }
  else {
    event_vtx_index = -1;
    event_vtx_dist  = -1;
  }

  return event_passes_cut[kAllCuts];

}
