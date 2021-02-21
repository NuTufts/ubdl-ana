import os,sys
import ROOT as rt
rt.gStyle.SetOptStat(0)

pot = 4.5e19
#pot = 2e20

ftypes = {"bnbnu":"plots_1e1p_sel_bnbnu_run3_merged.root",
          "intrinsicnue":"plots_1e1p_sel_intrinsic_merged.root"}
fcolors = {"intrinsicnue":rt.kRed-3,
           "bnbnu":rt.kBlue-3}
fpot = {"bnbnu":4.159516125e+19,
        "intrinsicnue":4.613240272e+22}
fill_order = ["intrinsicnue","bnbnu"]        

sample_names = [ "is1eVA","1e1p","all" ]

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
for s in sample_names:
    canvas[s] = rt.TCanvas("c%s"%(s),s,1200,600)
    # load hists
    hname = "hnshower_%s"%(s)
    hstack = rt.THStack(hname+"_stack","")    
    for nutype in fill_order:
        h = rfile[nutype].Get( hname )
        h.SetFillColor(fcolors[nutype])
        h.SetFillStyle(3003)
        scale = pot/fpot[nutype]
        h.Scale(scale)
        print nutype," ",s,": ",h.Integral()
        hstack.Add(h)        
        hists[(nutype,s)] = h
    
    hstack.Draw("hist")
    pred = h.Integral()
    canvas[s].Draw()
    print hname,": ",pred
    hstack_v[s] = hstack
    
print "enter to end"
raw_input()


