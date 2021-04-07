import os,sys
import ROOT as rt
rt.gStyle.SetOptStat(0)

pot = 5.0e19
#pot = 1e20

#ftypes = {"bnbnu":"plots_1e1p_sel_bnbnu_run3_v1reco.root",
#          "intrinsicnue":"plots_1e1p_sel_intrinsic_v1reco.root"}
ftypes = {"bnbnu":"plots_1e1p_sel_bnbnu_run3_v1reco_1trackmin.root",
          "intrinsicnue":"plots_1e1p_sel_intrinsic_v1reco_1trackmin.root"}

fcolors = {"intrinsicnue":rt.kRed-3,
           "bnbnu":rt.kBlue-3}
fpot = {"bnbnu":1.9611765276e+20,
        "intrinsicnue":4.4204107333999975e+22}
fill_order = ["intrinsicnue","bnbnu"]
#fill_order = ["bnbnu"]
#fill_order = ["intrinsicnue"]

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
nu_symbols = {"intrinsicnue":"#nu_{e}",
              "bnbnu":"#nu_{#mu}",
              "nc":"#nu"}

rfile = {}
for nutype in ftypes:
    rfile[nutype] = rt.TFile( ftypes[nutype], 'open' )

canvas = {}
hists = {}
hstack_v = {}
tlen_v = []
for s in sample_names:
    print "====== SAMPLE ",s," =========="
    canvas[s] = rt.TCanvas("c%s"%(s),s,1200,600)

    # building up this stack
    hstack = rt.THStack("h"+s+"_stack","")
    all_tot = 0. # a check
    stack_tot = 0.
    lowe_tot = 0.;

    tlen = rt.TLegend(0.75,0.35, 0.9,0.9)
    
    # load hists
    for nutype in fill_order:
        print "--",nutype,"-------"
        for mode in mode_names:
            #hname = "hEnu_%s_vertexactcut_%s"%(s,mode)
            hname = "hEnu_%s_allrecocut_%s"%(s,mode)        

            h = rfile[nutype].Get( hname )
            h.SetFillColor(fcolors[nutype])
            h.SetFillStyle(mode_styles[mode])
            scale = pot/fpot[nutype]
            h.Scale(scale)
            #print nutype," ",s," ",mode,": ",h.Integral()
            if mode != "all":
                hstack.Add(h)
                stack_tot += h.Integral()
                lowe_tot +=  h.Integral(1,h.GetXaxis().FindBin(500))
            else:
                all_tot += h.Integral()
                #print "  all <=500 Mev: ",h.Integral(1,h.GetXaxis().FindBin(500))
            #if nutype in ["intrinsicnue"]:
            print nutype," ",s," ",mode,": ",h.Integral()," low-E (<500 MeV): ",h.Integral(1,h.GetXaxis().FindBin(450))
            if h.Integral()>0 and mode!="all":
                if "nc" in mode:
                    tlen.AddEntry(h,"%s:%s"%( nu_symbols["nc"],mode),"F")
                else:
                    tlen.AddEntry(h,"%s:%s"%( nu_symbols[nutype],mode),"F")
            hists[(nutype,s,mode)] = h
            
    hstack.Draw("hist")
    hstack.SetTitle(";true E_{#nu};counts/100 MeV/1E20 POT")
    print "stack tot=", stack_tot," vs. all tot=",all_tot
    print "low E, <500 MeV tot=",lowe_tot
    tlen.Draw()
    canvas[s].Draw()
    hstack_v[s] = hstack
    canvas[s].Update()
    canvas[s].SaveAs("cplotenu_%s.png"%(s))
    tlen_v.append(tlen)
    raw_input()
    
print "POT: ",pot
print "enter to end"



