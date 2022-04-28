import os,sys
import ROOT as rt
from math import sqrt

"""
Compare the distance from the true vertex for reco vertices before and after optimization using prong charge.
"""

# input file is a KPSRecoManager file
inputfile = "test_bnbnu_fitted_kpsrecomanagerana.root"
rfile = rt.TFile(inputfile,"open")
ttree = rfile.Get("KPSRecoManagerTree")
nentries = ttree.GetEntries()

output = rt.TFile("output_compare_orig_vs_fitted_vtx.root","recreate")
hdist2vtx_prefit = rt.TH1D("hdist2vtx_prefit","",15,0,3.0)
hdist2vtx_fit = rt.TH1D("hdist2vtx_fit","",15,0,3.0)
himprovement = rt.TH1D("himprovement","",31,-3.0,3.0)

for ientry in range(nentries):


    print("==============")
    print("ENTRY ",ientry)
    print("===============")
    ttree.GetEntry(ientry)
    
    # true vtx
    true_vtx = [ttree.vtx_sce_x, ttree.vtx_sce_y, ttree.vtx_sce_z]

    nmerged = ttree.numerged_v.size()
    nvetoed = ttree.nuvetoed_v.size()

    # look for the closest kptype=0 (neutrino) vertex
    best_merged_dist = 999999.;
    best_merged_pos  = None
    for i in range(nmerged):
        nuvtx = ttree.numerged_v.at(i)
        if nuvtx.keypoint_type!=0:
            continue
        dist = 0.
        for v in range(3):
            dist = (nuvtx.pos[v]-true_vtx[v])*(nuvtx.pos[v]-true_vtx[v])
        dist = sqrt(dist)
        if dist<best_merged_dist:
            best_merged_dist = dist
            best_merged_pos = [nuvtx.pos[v] for v in range(3)]

    best_vetoed_dist = 999999.;
    best_vetoed_pos  = None
    for i in range(nvetoed):
        nuvtx = ttree.nuvetoed_v.at(i)
        if nuvtx.keypoint_type!=0:
            continue
        dist = 0.
        for v in range(3):
            dist = (nuvtx.pos[v]-true_vtx[v])*(nuvtx.pos[v]-true_vtx[v])
        dist = sqrt(dist)
        if dist<best_vetoed_dist:
            best_vetoed_dist = dist
            best_vetoed_pos = [nuvtx.pos[v] for v in range(3)]

    if best_vetoed_pos is not None and best_merged_pos is not None:
        hdist2vtx_prefit.Fill( best_merged_dist )
        hdist2vtx_fit.Fill( best_vetoed_dist )
        himprovement.Fill( best_merged_dist-best_vetoed_dist )
            
    print("best merged: ",best_merged_dist)
    print("best vetoed: ",best_vetoed_dist)

output.Write()
