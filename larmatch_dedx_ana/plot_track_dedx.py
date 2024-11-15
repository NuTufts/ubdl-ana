import os,sys
import ROOT as rt

rt.gStyle.SetOptStat(0)

inputfile = "dedxana_mcc9_v29e_bnb_nu_overlay.root"

"""
TFile**		Proton_Muon_Range_dEdx_LAr_TSplines.root	
 TFile*		Proton_Muon_Range_dEdx_LAr_TSplines.root	
  KEY: TSpline3	sMuonRange2T;1	sMuonRange2T;CSDA Range [cm];T [MeV]
  KEY: TSpline3	sMuonT2dEdx;1	sMuonT2dEdx;T [MeV];dEdX [Mev/cm]
  KEY: TSpline3	sMuonRange2dEdx;1	sMuonRange2dEdx;T [MeV];dEdX [Mev/cm]
  KEY: TSpline3	sProtonRange2T;1	sProtonRange2T;CSDA Range [cm];T [MeV]
  KEY: TSpline3	sProtonT2dEdx;1	sProtonT2dEdx;T [MeV];dEdX [Mev/cm]
  KEY: TSpline3	sProtonRange2dEdx;1	sProtonRange2dEdx;T [MeV];dEdX [Mev/cm]
"""

q2adc = 93.0/2.2

splinefile = rt.TFile( "Proton_Muon_Range_dEdx_LAr_TSplines.root" )
sMuonRange2dEdx = splinefile.Get("sMuonRange2dEdx")
sProtonRange2dEdx = splinefile.Get("sProtonRange2dEdx")

infile = rt.TFile( inputfile, "open" )
anadedx = infile.Get( "anadedx" )

outfile = rt.TFile( "temp.root", "recreate" )
hresrange_mu_med = rt.TH2D("hresrange_mu_med",";residual range (cm); #frac{dQ}{dx} (pixel sum per cm)",50,0,50,175,0,350)
hresrange_p_med  = rt.TH2D("hresrange_p_med",";residual range (cm); #frac{dQ}{dx} (pixel sum per cm)",50,0,50,175,0,350)

mu_curve = rt.TGraph( 50 )
for i in xrange(50):
    x = i*1.0+0.5
    y = sMuonRange2dEdx.Eval(x)
    y2 = q2adc*y/(1+y*0.0486/0.273/1.38)    
    mu_curve.SetPoint(i,x,y2)
mu_curve.SetLineColor(rt.kCyan)
mu_curve.SetLineWidth(2)
    
p_curve = rt.TGraph( 50 )
for i in xrange(50):
    x = i*1.0+0.5
    y = sProtonRange2dEdx.Eval(x)
    y2 = q2adc*y/(1+y*0.0486/0.273/1.38)
    p_curve.SetPoint(i,x,y2)
p_curve.SetLineColor(rt.kMagenta)
p_curve.SetLineWidth(2)

c = rt.TCanvas("c","c",1400,600)
c.Divide(2,1)

c.cd(1).SetLeftMargin(0.12)
c.cd(1).SetTopMargin(0.05)
    
hresrange_mu_med.GetYaxis().SetTitleOffset(1.7)
hresrange_mu_med.GetXaxis().SetTitleOffset(1.15)
anadedx.Draw("dqdx_med:res>>hresrange_mu_med","pid==13 && rad<3.0 && lm>0.8","colz")
mu_curve.Draw("L")
p_curve.Draw("L")

c.cd(2).SetLeftMargin(0.12)
c.cd(2).SetTopMargin(0.05)

hresrange_p_med.GetYaxis().SetTitleOffset(1.7)
hresrange_p_med.GetXaxis().SetTitleOffset(1.15)    
anadedx.Draw("dqdx_med:res>>hresrange_p_med","pid==2212 && rad<3.0 && lm>0.8","colz")
mu_curve.Draw("L")
p_curve.Draw("L")

c.Update()

raw_input()
