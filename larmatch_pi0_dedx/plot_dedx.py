import os,sys
import ROOT as rt

rt.gStyle.SetOptStat(0)

inputfile = "pi0_dedxana.root"

q2adc = 93.0/2.2

infile = rt.TFile( inputfile, "open" )
anadedx = infile.Get( "recodedx" )

outfile = rt.TFile( "temp.root", "recreate" )
nbinsy = 50
xend   = 10
hresrange_el_med = rt.TH2D("hresrange_el_med",";residual range (cm); #frac{dQ}{dx} (pixel sum per cm)",100,0,xend,nbinsy,0,350)
hresrange_ga_med = rt.TH2D("hresrange_ga_med",";residual range (cm); #frac{dQ}{dx} (pixel sum per cm)",100,0,xend,nbinsy,0,350)

c = rt.TCanvas("c","c",1400,600)
c.Divide(2,1)

c.cd(1).SetLeftMargin(0.12)
c.cd(1).SetTopMargin(0.05)
    
hresrange_el_med.GetYaxis().SetTitleOffset(1.7)
hresrange_el_med.GetXaxis().SetTitleOffset(1.15)
anadedx.Draw("dqdx_med:dist2start>>hresrange_el_med","rad<3.0 && lm>0.8 && dist2vertex_shower<10.0","colz")
#el_curve.Draw("L")
#ga_curve.Draw("L")

c.cd(2).SetLeftMargin(0.12)
c.cd(2).SetTopMargin(0.05)

hresrange_ga_med.GetYaxis().SetTitleOffset(1.7)
hresrange_ga_med.GetXaxis().SetTitleOffset(1.15)    
anadedx.Draw("dqdx_med:dist2start>>hresrange_ga_med","rad<3.0 && lm>0.8 && qtot_shower>3000.0 && dist2vertex_shower>10.0","colz")
#el_curve.Draw("L")
#ga_curve.Draw("L")

c.Update()

raw_input()
