import os,sys
import ROOT as rt
rt.gStyle.SetOptStat(0)

pot = 4.5e19
#pot = 2e20

ftypes = {"bnbnu":"plots_1e1p_sel_bnbnu_run3_merged.root",
          "intrinsicnue":"plots_1e1p_sel_intrinsic_merged.root"}
fcolors = {"intrinsicnue":rt.kRed-3,
           "bnbnu":rt.kBlue-3}
fpot = {"bnbnu":1.2495750974e+20,#5.112648436e+19,
        "intrinsicnue":4.613240272e+22}
fill_order = ["intrinsicnue","bnbnu"]
#fill_order = ["bnbnu"]

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

    # building up this stack
    hstack = rt.THStack("h"+s+"_stack","")
    all_tot = 0. # a check
    stack_tot = 0.
    
    # load hists
    for mode in mode_names:
        #hname = "hEnu_%s_vertexactcut_%s"%(s,mode)
        hname = "hEnu_%s_allrecocut_%s"%(s,mode)        

        for nutype in fill_order:
            h = rfile[nutype].Get( hname )
            h.SetFillColor(fcolors[nutype])
            h.SetFillStyle(mode_styles[mode])
            scale = pot/fpot[nutype]
            h.Scale(scale)
            print nutype," ",s," ",mode,": ",h.Integral()
            if mode != "all":
                hstack.Add(h)
                stack_tot += h.Integral()
            else:
                all_tot += h.Integral()
                print "  all <=500 Mev: ",h.Integral(1,h.GetXaxis().FindBin(500))
            hists[(nutype,s,mode)] = h
            
    hstack.Draw("hist")
    print "stack tot=", stack_tot," vs. all tot=",all_tot        
    canvas[s].Draw()
    hstack_v[s] = hstack
    
print "enter to end"
raw_input()


