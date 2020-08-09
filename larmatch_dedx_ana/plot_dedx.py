import os,sys
import ROOT as rt

rt.gStyle.SetOptStat(0)

inputfile = "dedxana_mcc9_v29e_bnb_nu_overlay.root"


infile = rt.TFile( inputfile, "open" )
anadedx = infile.Get( "anadedx" )

c = rt.TCanvas("c","c",1400,600)
c.Divide(2,1)

c.cd(1)
anadedx.Draw("dqdx_med:res>>hresrange_mu_med(50,0,50,175,0,350)","pid==13 && rad<3.0 && lm>0.8","colz")

c.cd(2)
anadedx.Draw("dqdx_med:res>>hresrange_p_med(50,0,50,175,0,350)","pid==2212 && rad<3.0 && lm>0.8","colz")


raw_input()
