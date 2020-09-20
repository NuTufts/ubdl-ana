import os,sys
import ROOT as rt

rt.gStyle.SetOptStat(0)

samples = ["data","mc"]
particles = ["all","p","mu"]

inputfile = { "data":"dedxana_1m1p.root",
              "mc":"output_dedxanamc_bnb_nu_overlay_1m1p.root" }


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
nbins = 101

infile = {}
ttree = {}
for sample in samples:
    infile[sample] = rt.TFile( inputfile[sample], "open" )
    ttree[sample] = infile[sample].Get("llana")

outfile = rt.TFile( "temp.root", "recreate" )
hllr = {}
cut = {"all":"pid!=111",
       "p":"pid!=111 && abs(truth_pid)==2212 && truth_mse<1.0",
       "mu":"pid!=111 && abs(truth_pid)==13 && truth_mse<1.0"}
color = {"all":rt.kBlack,
         "p":rt.kBlue,
         "mu":rt.kRed}

ctemp = rt.TCanvas("ctemp","ctemp",800,600)
for pid in particles:
    for sample in samples:
        if sample=="data" and pid!="all":
            continue
        
        hname = "hllr_%s_%s"%(pid,sample)
        hllr[(pid,sample)] = rt.TH1D(hname,";log-likelihood ratio; events",nbins,-100,100)
        ttree[sample].Draw("llpid>>%s"%(hname),cut[pid])
        hllr[(pid,sample)].SetLineColor(color[pid])
        if pid!="all":
            hllr[(pid,sample)].SetFillColor(color[pid])
            hllr[(pid,sample)].SetFillStyle(3003)

mcscale = hllr[("all","data")].Integral()/hllr[("all","mc")].Integral()
for pid in particles:
    hllr[(pid,"mc")].Scale(mcscale)

#mu_curve = rt.TGraph( int(xend) )
#for i in xrange(int(xend)):
#    x = i*1.0+0.5
#    y = sMuonRange2dEdx.Eval(x)
#    y2 = q2adc*y/(1+y*0.0486/0.273/1.38)    
#    mu_curve.SetPoint(i,x,y2)
#mu_curve.SetLineColor(rt.kCyan)
#mu_curve.SetLineWidth(2)
    
#p_curve = rt.TGraph( int(xend) )
#for i in xrange(int(xend)):
#    x = i*1.0+0.5
#    y = sProtonRange2dEdx.Eval(x)
#    y2 = q2adc*y/(1+y*0.0486/0.273/1.38)
#    p_curve.SetPoint(i,x,y2)
#p_curve.SetLineColor(rt.kMagenta)
#p_curve.SetLineWidth(2)

c = rt.TCanvas("c","c",800,600)

c.cd(1).SetLeftMargin(0.12)
c.cd(1).SetTopMargin(0.05)
    
#hresrange_mu_med.GetYaxis().SetTitleOffset(1.7)
#hresrange_mu_med.GetXaxis().SetTitleOffset(1.15)
hllr[("all","mc")].Draw("hist")
hllr[("p","mc")].Draw("histsame")
hllr[("mu","mc")].Draw("histsame")
hllr[("all","data")].Draw("E1same")

cutbin = hllr[("p","mc")].GetXaxis().FindBin(0.0)
npass_p  = hllr[("p","mc")].Integral(0,cutbin)
npass_mu = hllr[("mu","mc")].Integral(cutbin+1,nbins+1)

nfail_p  = hllr[("p","mc")].Integral(cutbin+1,nbins+1)
nfail_mu = hllr[("mu","mc")].Integral(0,cutbin)

p_eff = npass_p/hllr[("p","mc")].Integral()
p_pur = npass_p/(npass_p+nfail_mu)

mu_eff = npass_mu/hllr[("mu","mc")].Integral()
mu_pur = npass_mu/(npass_mu+nfail_p)

print "proton eff: ",p_eff
print "muon eff: ",mu_eff
print "proton purity: ",p_pur
print "muon purity: ",mu_pur


c.Update()

raw_input()
