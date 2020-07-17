import os,sys
import ROOT as rt

rt.gStyle.SetOptStat(0)

# input ana file
infile = rt.TFile( sys.argv[1], "open" )
ev_ana = infile.Get("nuvertexana_event")
vtx_ana = infile.Get("vertexana_recovertex")

# CUTS
fvcut = "(vtx_x>10 && vtx_x<246 && vtx_y>-106 && vtx_y<106 && vtx_z>10.0 && vtx_z<1026.0)"
selection = [("all","1==1"),
             ("1l1p","is1l1p0pi==1"),
             ("1shower","nlepton_35mev==1 && nproton_60mev==0 && nmeson_35mev==0 && npi0==0")]

# Eff. versus Enu
ceff = {}
heff = {}
teff = {}
eff_nbins = 40
eff_range = (0,400)
eff_plot = [("eff_enu_all","Enu_true","1==1",rt.kBlue,"All"),
            ("eff_enu_kp","Enu_true","min_dist_to_vtx<5.0",rt.kBlue,"Keypoint Net"),
            ("eff_enu_wc","Enu_true","min_dist_to_vtx_wct<5.0",rt.kCyan,"Keypoint Net+WC filter"),
            ("eff_enu_dl","Enu_true","min_dist_to_vtx_dl<5.0",rt.kRed,"DLLEE Vertexer")]

for (cutname,select_cut) in selection:
    ceff[cutname] = rt.TCanvas("ceff_"+cutname,"Eff vs. true Enu, %s"%(cutname),800,400)
    for (varname,var,varcut,color,label) in eff_plot:
        hname = "h"+varname+"_"+cutname
        heff[hname] = rt.TH1F(hname,";true neutrino energy (MeV)",eff_nbins, eff_range[0], eff_range[1])
        cut = "(" + fvcut + ") && (" + select_cut + ") && (" + varcut +")"
        print var,":: ",cut
        ev_ana.Draw("%s>>%s"%(var,hname),cut)
        heff[hname].SetLineColor(color)
    hallname = "h"+eff_plot[0][0]+"_"+cutname
    teff[cutname] = rt.TLegend(0.1,0.7,0.4,0.9)    
    for (varname,var,varcut,color,label) in eff_plot[1:]:
        hname = "h"+varname+"_"+cutname
        print "Divide ",hname," by ",hallname
        heff[hname].Divide( heff[hallname] )
        teff[cutname].AddEntry( heff[hname], label, "L" )
        

for (cutname,select_cut) in selection:
    ceff[cutname].Clear()
    c = ceff[cutname]
    c.SetBottomMargin(0.15)
    c.SetLeftMargin(0.1)

    h = heff["heff_enu_kp_"+cutname]
    h.SetTitle(";true neutrino energy (MeV); efficiency for %s events"%(cutname))
    h.GetYaxis().SetRangeUser(0,1)
    h.GetXaxis().SetTitleOffset(1.0)
    h.GetYaxis().SetTitleOffset(0.75)    
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetYaxis().SetLabelSize(0.05)

    
    heff["heff_enu_kp_"+cutname].Draw()
    heff["heff_enu_wc_"+cutname].Draw("same")
    heff["heff_enu_dl_"+cutname].Draw("same")

    teff[cutname].Draw()
    
    ceff[cutname].Update()
    ceff[cutname].SaveAs("ceff_vs_trueenu_%s.png"%(cutname))
    

raw_input()
