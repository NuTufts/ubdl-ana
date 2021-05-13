import os,sys
import ROOT as rt

rt.gStyle.SetOptStat(0)

input_nue  = rt.TFile("output_eff_plots_intrinsic_nue.root")
input_numu = rt.TFile("output_eff_plots_bnbnu.root")


h1e1p = input_nue.Get("heff_enu_dl_1e1p")
h1m1p = input_numu.Get("heff_enu_dl_1e1p")

hspec_1e1p  = input_nue.Get("heff_enu_all_1e1p")
hspec_1e1p.SetFillColor(rt.kRed)
hspec_1e1p.SetFillStyle(3003)
tot_1e1p = hspec_1e1p.Integral()
hspec_1e1p.Scale( 1.0/tot_1e1p )

hspec_1m1p = input_numu.Get("heff_enu_all_1e1p")
hspec_1m1p.SetFillColor(rt.kBlue+2)
hspec_1m1p.SetFillStyle(3003)
tot_1m1p = hspec_1m1p.Integral()
hspec_1m1p.Scale( 1.0/tot_1m1p )


c = rt.TCanvas("c","c",800,400)
t = rt.TLegend(0.2,0.75,0.5,0.85)
h1e1p.Draw("histE1")
h1e1p.SetTitle(";true neutrino energy (MeV);efficiency")
h1e1p.GetYaxis().SetRangeUser(0,1.0)
h1m1p.SetLineColor(rt.kBlue+2)
h1m1p.Draw("histE1same")
#hspec_1m1p.Draw("histsame")
#hspec_1e1p.Draw("histsame")
t.AddEntry(h1e1p,"true 1e1p events")
t.AddEntry(h1m1p,"true 1#mu1p events")
t.SetBorderSize(0)
t.Draw()
c.Draw()
c.Update()
input()
