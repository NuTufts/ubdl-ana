import os,sys
import ROOT as rt
rt.gStyle.SetOptStat(0)

rt.gStyle.SetHistTopMargin(0)

pot = 4.5e19
#pot = 1e20

ftypes = {"bnbnu":"plots_1e1p_sel_bnbnu_run3_v1reco.root",
          "intrinsicnue":"plots_1e1p_sel_intrinsic_v1reco.root"}
fcolors = {"intrinsicnue":rt.kRed-3,
           "bnbnu":rt.kBlue-3}
fpot = {"bnbnu":5.409417753e+19,
        "intrinsicnue":4.4204107333999975e+22}
#fill_order = ["intrinsicnue","bnbnu"]
#fill_order = ["bnbnu"]
fill_order = ["intrinsicnue"]

sample_names = [ "is1eVA","1e1p","all" ]

reco_state_names = ["good","bad"]
reco_state_styles = {"good":1001,
                     "bad":3002}

cutvar_names  = [ "dwall",
                  "dist2true",
                  "nmaxshowerhits",
                  "nshowerprongs",
                  "ntrackprongs",
                  "llpid",
                  "hipfraction",
                  "minshowergap",
                  "maxshowergap",
                  "maxtracklen",
                  "vertexcharge",
                  "largestshowerll",
                  "closestshowerll",
                  "largestshoweravedqdx",
                  "closestshoweravedqdx",
                  "unrecoqmedfrac"]


mode_names = ["ccqe",
              "ccres",
              "ccother",
              "ncqe",
              "ncres",
              "ncother",
              "all"]

mode_styles = {"ccqe":1001,
               "ccres":3002,
               "ccother":3144,
               "ncqe":3006,
               "ncres":3007,
               "ncother":3025,
               "all":3003}


rfile = {}
for nutype in ftypes:
    rfile[nutype] = rt.TFile( ftypes[nutype], 'open' )

canvas = {}
hists = {}
hstack_v = {}
tlen_v = []
for cutvar in cutvar_names:
    print "========[ ",cutvar," ]==========="
    canvas[cutvar] = rt.TCanvas("c%s"%(cutvar),cutvar,1200,600)

    # build this stacked histogram
    hstack = rt.THStack("hstack_"+cutvar,cutvar)
    stack_tot = 0.
    
    # load hists

    for recostat in reco_state_names:
        hname = "hCutVar_"+cutvar+"_"+recostat

        for nutype in fill_order:
            h = rfile[nutype].Get( hname )
            h.SetFillColor(fcolors[nutype])
            h.SetFillStyle( reco_state_styles[recostat] )
            scale = pot/fpot[nutype]
            h.Scale(scale)
            print nutype," ",hname,": ",h.Integral()
            hstack.Add(h)
            stack_tot += h.Integral()
            hists[(nutype,cutvar,recostat)] = h
    

    ymax = hstack.GetMaximum()
    print ymax
    hstack.SetMaximum(ymax*10 )
    hstack.SetMinimum( 0.01 )
    hstack.Draw("hist")    
    canvas[cutvar].Draw()
    print hname,": ",stack_tot
    hstack_v[cutvar] = hstack
    canvas[cutvar].SetLogy(1)
    canvas[cutvar].Update()
    
print "enter to end"
raw_input()


