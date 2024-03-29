import os,sys
import ROOT as rt
from math import sqrt

rt.gStyle.SetOptStat(0)

# input ana file
infile = rt.TFile( sys.argv[1], "open" )
ev_ana = infile.Get("nuvertexana_event")
vtx_ana = infile.Get("vertexana_recovertex")

# CUTS
fvcut = "(vtx_x>10 && vtx_x<246 && vtx_y>-106 && vtx_y<106 && vtx_z>10.0 && vtx_z<1026.0)"
selection = [("all","1==1"),
             #("1#mu1p","is1l1p0pi==1"),
             ("1e1p","is1l1p0pi==1"),
             ("1shower","nlepton_35mev==1 && nproton_60mev==0 && nmeson_35mev==0 && npi0==0")]

# Eff. versus Enu
ceff = {}
heff = {}
teff = {}
eff_nbins = 12
eff_range = (0,1200)
good_rad = 3.0
eff_plot = [("eff_enu_all","Enu_true","1==1",rt.kBlue,"All",True),
            ("eff_enu_kp","Enu_true","min_dist_to_vtx<%.2f"%(good_rad),rt.kBlue,"Keypoint Net",False),
            #("eff_enu_wc","Enu_true","min_dist_to_vtx_wct<5.0",rt.kCyan,"Keypoint Net+WC filter",False),
            ("eff_enu_dl","Enu_true","min_dist_to_vtx_dl<%.2f"%(good_rad),rt.kRed,"DLLEE Vertexer",True)]

out = rt.TFile("output_eff_plots.root","recreate")

for (cutname,select_cut) in selection:
    ceff[cutname] = rt.TCanvas("ceff_"+cutname,"Eff vs. true Enu, %s"%(cutname),800,400)
    for (varname,var,varcut,color,label,plotme) in eff_plot:
        hname = "h"+varname+"_"+cutname
        heff[hname] = rt.TH1F(hname,";true neutrino energy (MeV)",eff_nbins, eff_range[0], eff_range[1])
        cut = "(" + fvcut + ") && (" + select_cut + ") && (" + varcut +")"
        print var,":: ",cut
        ev_ana.Draw("%s>>%s"%(var,hname),cut)
        heff[hname].SetLineColor(color)
    hallname = "h"+eff_plot[0][0]+"_"+cutname
    teff[cutname] = rt.TLegend(0.1,0.7,0.4,0.9)
    teff[cutname].SetBorderSize(0)
    for (varname,var,varcut,color,label,plotme) in eff_plot[1:]:
        hname = "h"+varname+"_"+cutname
        print "Divide ",hname," by ",hallname
        heff[hname].Divide( heff[hallname] )
        for ibin in range(heff[hname].GetXaxis().GetNbins()):
            binval = heff[hname].GetBinContent(ibin+1)
            binN   = heff[hallname].GetBinContent(ibin+1)
            k = binval*binN
            if binN>0:
                deff = (1.0/float(binN))*sqrt( k*(1-binval) )
                heff[hname].SetBinError(ibin+1,deff)
            else:                
                heff[hname].SetBinError(ibin+1,0.0)

        if plotme:
            teff[cutname].AddEntry( heff[hname], label, "L" )
        

for (cutname,select_cut) in selection:
    ceff[cutname].Clear()
    c = ceff[cutname]
    c.SetBottomMargin(0.15)
    c.SetLeftMargin(0.1)

    #h = heff["heff_enu_kp_"+cutname]
    h = heff["heff_enu_dl_"+cutname]
    h.SetTitle(";true neutrino energy (MeV); efficiency for %s events"%(cutname))
    h.GetYaxis().SetRangeUser(0,1.0)
    h.GetXaxis().SetTitleOffset(1.0)
    h.GetYaxis().SetTitleOffset(0.75)    
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetYaxis().SetLabelSize(0.05)

    
    #heff["heff_enu_kp_"+cutname].Draw()
    #heff["heff_enu_wc_"+cutname].Draw("same")
    heff["heff_enu_dl_"+cutname].Draw("histE1")

    teff[cutname].Draw()
    
    ceff[cutname].Update()
    ceff[cutname].SaveAs("ceff_vs_trueenu_%s.png"%(cutname))
    
out.Write()
raw_input()
