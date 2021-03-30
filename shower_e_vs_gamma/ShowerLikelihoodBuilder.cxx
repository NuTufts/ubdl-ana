#include "ShowerLikelihoodBuilder.h"

#include "larflow/Reco/geofuncs.h"

namespace larflow {
namespace reco {

  ShowerLikelihoodBuilder::ShowerLikelihoodBuilder()
    : larcv::larcv_base("ShowerLikelihoodBuilder")
  {
    _make_profile_hists();
  }
  
  void ShowerLikelihoodBuilder::fill( const std::vector<larlite::larflowcluster>& shower_v,
                                      const std::vector<larlite::track>& shower_trunk_v,                                      
                                      const std::vector<larcv::Image2D>& adc_v ) {

    for (int ishower=0; ishower<(int)shower_v.size(); ishower++ ) {
      
      auto const& shower = shower_v[ishower];
      auto const& trunk  = shower_trunk_v[ishower];

      float adcsum = _calcSum( shower, adc_v, 2, 10.0 );

      float MeV_reco = 400.0/100.0*adcsum;

      TH2F* hprofile = getProfileHist( MeV_reco );

      std::vector<float> shower_vtx(3,0);
      std::vector<float> shower_dir(3,0);

      const TVector3 start = trunk.LocationAtPoint(0);
      const TVector3 end   = trunk.LocationAtPoint( trunk.NumberTrajectoryPoints()-1 );
      TVector3 dir = end-start;

      for (int i=0; i<3; i++) {
        shower_vtx[i] = start[i];
        shower_dir[i] = (end[i]-start[i])/dir.Mag();
      }

      _fillProfileHist( shower, shower_dir, shower_vtx, *hprofile );
      
    }

  }
  
  void ShowerLikelihoodBuilder::_make_profile_hists()
  {
    float min_energy = 0;
    float max_energy = 500;
    float step_energy = 25.0;
    int nsteps = 20;

    for (int istep=0; istep<nsteps; istep++) {
      float E1 = min_energy + istep*step_energy;
      float E2 = min_energy + (istep+1)*step_energy;

      EnergyRange_t erange;
      erange.min = E1;
      erange.max = E2;
      char histname[50];
      sprintf( histname, "hprofile_%02d",istep);
      char hist_title[100];
      sprintf( hist_title, "Profile %d for %0.2f to %0.2f_Me", istep, E1, E2 );
      erange.hist = new TH2F( histname, hist_title, 2000, -10, 190, 1000, 0, 100 );

      _profile_maps[ erange ] = erange.hist;

    }
      
  }

  void ShowerLikelihoodBuilder::clearProfileHists()
  {
    for ( auto it=_profile_maps.begin(); it!=_profile_maps.end(); it++ ) {
      it->second->Reset();
    }
  }

  float ShowerLikelihoodBuilder::_calcSum( const larlite::larflowcluster& shower,
                                           const std::vector< larcv::Image2D >& adc_v,
                                           const int dpix,
                                           const float threshold )
  {
    
    auto const& meta = adc_v.front().meta();
    
    std::vector<float> pixsum_v( adc_v.size(), 0 );
    
    for ( auto const& hit : shower ) {
      
      float tick = hit.tick;
      
      if ( tick<=meta.min_y() || tick>=meta.max_y() ) {
        continue;
      }

      int row = meta.row(tick);
      std::vector<int> imgcoord(adc_v.size()+1,0);
      imgcoord[adc_v.size()] = row;

      for (int dr=-abs(dpix); dr<=abs(dpix); dr++) {
        int r = row+dr;
        for (int p=0; p<(int)adc_v.size(); p++) {

          int col = hit.targetwire[p];
          
          for (int dc=-abs(dpix); dc<=abs(dpix); dc++) {

            int c = col + dc;

            if ( c<0 || c>=(int)meta.cols() ) {

              float pixval = adc_v[p].pixel( r, c );
              if ( pixval>threshold ) {
                pixsum_v[p] += pixval;
              }
            }
          }//dc loop
        }//plane loop
      }//dr loop
      
    }//end of shower loop
    
    return pixsum_v[2];
  }
  
  void ShowerLikelihoodBuilder::_fillProfileHist( const larlite::larflowcluster& hit_v,
                                                  const std::vector<float>& shower_dir,
                                                  const std::vector<float>& shower_vtx,
                                                  TH2F& hprofile )
  {

    // get distance of point from pca-axis
    // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    
    LARCV_DEBUG() << "Fill hits for shower[ "
                  << "dir=(" << shower_dir[0] << "," << shower_dir[1] << "," << shower_dir[2] << ") "
                  << "vtx=(" << shower_vtx[0] << "," << shower_vtx[1] << "," << shower_vtx[2] << ") "
                  << "] with nhits=" << hit_v.size()
                  << std::endl;
    
    std::vector<float> endpt(3,0);
    for (int i=0; i<3; i++) {
      endpt[i] = shower_vtx[i] + 10*shower_dir[i];
    }
    
    for ( auto const& hit : hit_v ) {
      std::vector<float> pt = { hit[0], hit[1], hit[2] };
      float rad  = larflow::reco::pointLineDistance3f( shower_vtx, endpt, pt );
      float proj = larflow::reco::pointRayProjection3f( shower_vtx, shower_dir, pt );
      hprofile.Fill( proj, rad );
    }
    
  }
  

}
}
