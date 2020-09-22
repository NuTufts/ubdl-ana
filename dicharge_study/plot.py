import os,sys
import ROOT as rt

rt.gStyle.SetOptStat(0)

samples = ["highstat","detvarcv"]
color = {"highstat":rt.kBlue,
         "detvarcv":rt.kRed}

anafile = {}
ana = {}
h2d = {}
h1d = {}

c = rt.TCanvas("c","c",2000,1800)
c.Divide(3,3)
c.cd(1)



for s  in samples:
    filepath = "out_%s.root"%(s)
    anafile[s] = rt.TFile(filepath)
    ana[s] = anafile[s].Get("ana")
    for p in xrange(3):
        h2d[(s,p)] = rt.TH2D("hqperpix_%s_p%d"%(s,p),"Q/pixel: %s,plane %d;phi: |atan2(dir_y,dir_x)| (rad);Q/pix"%(s,p),101,-3.14169,3.14160,50,0,100)
        h1d[(s,p)] = rt.TH1D("hqperpix1d_%s_p%d"%(s,p),"Q/pixel for |#phi|<0.1, plane %d"%(p),50,0,100)
        ana[s].Draw("qperpix[%d]:phi>>hqperpix_%s_p%d"%(p,s,p),"qperpix[%d]<100 && abs(pid)==13"%(p),"colz")
        #ana[s].Draw("qperpix[%d]>>hqperpix1d_%s"%(s),"fabs(phi)<0.5 || fabs(3.14159-phi)<0.5","colz")
        ana[s].Draw("qperpix[%d]>>hqperpix1d_%s_p%d"%(p,s,p),"fabs(phi)<0.1 && abs(pid)==13","colz")    
        h1d[(s,p)].SetLineColor(color[s])
        tot = h1d[(s,p)].Integral()
        h1d[(s,p)].Scale(1.0/tot)

t_v  = []
for p in xrange(3):
    c.cd(3*p+0+1)
    h2d[("highstat",p)].Draw("colz")
    c.cd(3*p+1+1)
    h2d[("detvarcv",p)].Draw("colz")
    c.cd(3*p+2+1)
    if h1d[("highstat",p)].GetMaximum()>h1d[("detvarcv",p)].GetMaximum():
        h1d[("highstat",p)].Draw("histE1")
        h1d[("detvarcv",p)].Draw("histE1same")        
    else:
        h1d[("detvarcv",p)].Draw("histE1")        
        h1d[("highstat",p)].Draw("histE1same")
        
    t = rt.TLegend(0.5,0.5,0.9,0.9)        
    t.AddEntry(h1d[("detvarcv",p)],"DetVar CV","L")
    t.AddEntry(h1d[("highstat",p)],"HighStats CV","L")
    t.Draw()    
    t_v.append(t)


raw_input()
