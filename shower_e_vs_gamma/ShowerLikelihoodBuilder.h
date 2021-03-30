#ifndef __LARFLOW_RECO_SHOWERLIKELIHOODBUILDER_H__
#define __LARFLOW_RECO_SHOWERLIKELIHOODBUILDER_H__

/**
 * This class contains functions to build distibutions needed 
 * to make shower likelihood functions to be used in the larfow shower reco code
 *
 * This includes:
 *
 *  1) shower profile likelihood in 3D. its the location of charge deposited as a function of 
 *       the distance along the trunk line and the perpendicular dist from the trunk line
 *  2) brem segment impact param, distance to trunk line, cosine of pca between trunk lines
 *
 * We fill a tree to later use to make distributions.
 * Can feed it, single shower MC  (best) or low energy neutrino. OK.
 *
 */

#include <map>
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "DataFormat/storage_manager.h"
#include "DataFormat/larflow3dhit.h"
#include "DataFormat/larflowcluster.h"
#include "DataFormat/track.h"

#include "larflow/PrepFlowMatchData/PrepMatchTriplets.h"

#include "TFile.h"
#include "TH2F.h"

namespace larflow {
namespace reco {

  class ShowerLikelihoodBuilder : public larcv::larcv_base {

  public:

    ShowerLikelihoodBuilder();
    virtual ~ShowerLikelihoodBuilder() {};

    void process( larcv::IOManager& iolcv, larlite::storage_manager& ioll );

    TH2F* getProfileHist( float energy ) {
      for ( auto it=_profile_maps.begin(); it!=_profile_maps.end(); it++ ) {
        if ( it->first.isinside( energy ) )
          return it->second;
      }
      return nullptr;
    };

    void fill( const std::vector<larlite::larflowcluster>& shower_v,
               const std::vector<larlite::track>& shower_trunk_v,                                      
               const std::vector<larcv::Image2D>& adc_v );

    void clearProfileHists();    
    
  protected:

    std::vector< TH2F* > _hprofile;
    struct EnergyRange_t {
      float min;
      float max;
      TH2F* hist;
      bool isinside( float energy ) const {
        if ( min<=energy && energy<=max )
          return true;
        return false;
      };
      bool operator<( const EnergyRange_t& rhs ) const
      {
        if ( rhs.min < min )
          return true;
        return false;
      }
    };

    std::map< EnergyRange_t, TH2F* > _profile_maps;
    
    void _make_profile_hists();

    void _fillProfileHist( const larlite::larflowcluster& hit_v,
                           const std::vector<float>& shower_dir,
                           const std::vector<float>& shower_vtx,
                           TH2F& hprofile );

    float _calcSum( const larlite::larflowcluster& shower,
                    const std::vector< larcv::Image2D >& adc_v,
                    const int dpix, const float threshold );
    

    
    
  };
  
}
}

#endif
