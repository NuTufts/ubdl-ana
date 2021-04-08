#ifndef __NUESELECTION_H__
#define __NUESELECTION_H__

#include <string>
#include <vector>

#include "larflow/Reco/NuSelectionVariables.h"
#include "larflow/Reco/NuVertexCandidate.h"

class NueSelection {
 public:

  // cut stages
  enum { kFV=0,           // [0] true vertex in FV (10 cm from TPC boundary), sets baseline for efficiency study
         kVertexCand3cm,  // [1] reco candidate formed within 3 cm of vertex
         kMinShowerSize,  // [2] min shower size cut (might want to loosen)
         kNShowerProngs,  // [3] number of shower prongs
         kNTrackProngs,   // [4] number of track prongs         
         kShowerGap,      // [5] shower gap
         kTrackGap,       // [6] track gap
         kMaxTrackLen,    // [7] max track len
         kSecondShower,   // [8] second shower size
         kVertexAct,      // [9] vertex activity cut
         kRecoFV,         // [10] reco fv cut
         kShowerLLCut,    // [11] shower likelihood cut         
         kWCPixel,        // [12] Wire-Cell pixel cut
         kHadronic,       // [13] see hadronic particles (proton or vertex activity)
         kUnrecoQ,        // [14] how much of the in-time pixels have been used
         kAllCuts,        // [15] All cuts applied except FV -- represents reco pass rate
         kNumCuts };      // [16] Number in enum
  
  std::vector<std::string> selcut_names;
  
  NueSelection() {
    selcut_names
      = { "fv",             // [0]
          "vertexcand",     // [1]
          "minshower",      // [2]
          "nshowerprongs",  // [3]
          "ntrackprongs",   // [4]
          "showergap",      // [5]
          "trackgap",       // [6]
          "maxtracklen",    // [7]
          "secondshower",   // [8]
          "vertexact",      // [9]
          "showerll",       // [10]        
          "recofv",         // [11]
          "wcpixel",        // [12]        
          "hadronic",       // [13]
          "unrecoq",        // [14]
          "allreco",        // [15]          
          "numcuts"};       // [16]

    
  };
  virtual ~NueSelection() {};

  std::vector< std::vector<bool> > event_vtx_pass;
  std::vector<bool> event_passes_cut;
  int event_vtx_index;
  float event_vtx_dist;
  
  std::vector<bool>
    runVertexSelection(  const larflow::reco::NuSelectionVariables& nusel,
                         const larflow::reco::NuVertexCandidate& nuvtx,
                         bool cut_fv );

  bool runEventSelection( const std::vector< larflow::reco::NuSelectionVariables>& nusel_v,
                          const std::vector< larflow::reco::NuVertexCandidate>& nuvtx_v,
                          const float vtx_dwall,
                          const int ientry );                          
  
  bool DEBUG;
  
};

#endif
