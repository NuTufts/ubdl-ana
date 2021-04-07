import os,sys
import numpy as np
import ROOT as rt
rt.gStyle.SetOptStat(0)

#pot = 4.5e19
pot = 1e20

ftypes = {"bnbnu":"plots_final_states_bnbnu_run3_merged.root",
          "intrinsicnue":"plots_final_states_intrinsic_merged.root"}
fcolors = {"intrinsicnue":rt.kRed-3,
           "bnbnu":rt.kBlue-3}
fpot = {"bnbnu":2.2610676018999955e+20,
        "intrinsicnue":4.597582955e+22}
fill_order = ["intrinsicnue","bnbnu"]
#fill_order = ["bnbnu"]

finalstate_names = [ "1eVA",
                     "1e_1p_0pi_nopi0",  "1e_0p_1pi_nopi0", "1e_0p_0pi_wpi0",
                     "1e_Np_0pi_nopi0",  "1e_1p_Npi_nopi0", "1e_1p_0pi_wpi0",
                     "1e_Np_Npi_wpi0",
                     "1mVA",
                     "1m_1p_0pi_nopi0",  "1m_0p_1pi_nopi0", "1m_0p_0pi_wpi0",
                     "1m_Np_0pi_nopi0",  "1m_1p_Npi_nopi0", "1m_1p_0pi_wpi0",
                     "1m_Np_Npi_wpi0",
                     # NC
                     "nc_1p_0pi_nopi0",  "nc_0p_1pi_nopi0", "nc_0p_0pi_wpi0",
                     "nc_Np_0pi_nopi0",  "nc_1p_Npi_nopi0", "nc_1p_0pi_wpi0",
                     "nc_Np_Npi_wpi0",
                     # cosmic
                     "offnu" ]
finalstate_colors = np.arange(1,50)
print finalstate_colors
np.random.shuffle( finalstate_colors )
print finalstate_colors

nu_symbols = {"intrinsicnue":"#nu_{e}",
              "bnbnu":"#nu_{#mu}",
              "nc":"#nu"}
nu_styles = {"intrinsicnue":1001,
             "bnbnu":3001,
             "nc":3144}

mode_names = ["all"]
mode_styles = {"ccqe":1001,
               "ccres":3002,
               "ccother":3144,
               "ncqe":3006,
               "ncres":3007,
               "ncother":3025,
               "all":1001}

rfile = {}
for nutype in ftypes:
    rfile[nutype] = rt.TFile( ftypes[nutype], 'open' )

canvas = {}
hists = {}
hstack_v = {}
tlen_v = []

for mode in mode_names:
    print "====== MODE: ",mode," =========="

    canvas[mode] = rt.TCanvas("c%s"%(mode),mode,1200,600)

    # building up this stack
    hstack = rt.THStack("h"+mode+"_stack","")
    all_tot = 0. # a check
    stack_tot = 0.
    lowe_tot = 0.;

    tlen = rt.TLegend(0.75,0.35, 0.9,0.9)
    
    # load hists
    for nutype in fill_order:
        print "--",nutype,"-------"
        for ifs,fs in enumerate(finalstate_names):
            #hname = "hEnu_%s_vertexactcut_%s"%(s,mode)
            hname = "hEnu_%s_allrecocut_%s"%(fs,mode)        

            h = rfile[nutype].Get( hname )
            #h.SetFillColor(fcolors[nutype])
            h.SetFillColor( finalstate_colors[ifs] )
            h.SetFillStyle( nu_styles[nutype] )
            scale = pot/fpot[nutype]
            h.Scale(scale)
            #print nutype," ",s," ",mode,": ",h.Integral()
            if mode == "all":
                hstack.Add(h)
                stack_tot += h.Integral()
                lowe_tot +=  h.Integral(1,h.GetXaxis().FindBin(500))
            else:
                all_tot += h.Integral()
                #print "  all <=500 Mev: ",h.Integral(1,h.GetXaxis().FindBin(500))
            #if nutype in ["intrinsicnue"]:
            print nutype," ",fs," ",mode," allreco: ",h.Integral(),": ",h.Integral(1,h.GetXaxis().FindBin(500))
            if h.Integral()>0 and mode=="all":
                if "nc" in mode:
                    tlen.AddEntry(h,"%s:%s (%.1f)"%( nu_symbols["nc"],fs,h.Integral()),"F")
                else:
                    tlen.AddEntry(h,"%s:%s (%.1f)"%( nu_symbols[nutype],fs,h.Integral()),"F")
            hists[(nutype,fs,mode)] = h
            
    hstack.Draw("hist")
    hstack.SetTitle(";true E_{#nu};counts/100 MeV/1E20 POT")
    print "stack tot=", stack_tot," vs. all tot=",all_tot
    print "low E, <500 MeV tot=",lowe_tot
    tlen.Draw()
    canvas[mode].Draw()
    hstack_v[mode] = hstack
    canvas[mode].Update()
    canvas[mode].SaveAs("cplotenu_%s.png"%(mode))
    tlen_v.append(tlen)
    raw_input()
    
print "POT: ",pot
print "enter to end"



