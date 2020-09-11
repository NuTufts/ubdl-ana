import os,sys
import ROOT as rt

rt.gStyle.SetOptStat(0)

infile_nu = sys.argv[1]
infile_cosmic  = sys.argv[2]

# compare false vertex rate for neutrino MC file with vertexing rate of EXTBNB sample
nu_file = rt.TFile( infile_nu, "open" )
bg_file = rt.TFile( infile_cosmic, "open" )

nu_tree = nu_file.Get( "nuvertexana_event" )
bg_tree = bg_file.Get( "nuvertexana_event" )


# SAMPLE
samples = [("nu",nu_tree,rt.kRed),
           ("bg",bg_tree,rt.kBlack)]


# Eff. versus Enu
cbadrate = {}
hbadrate = {}
tbadrate = {}
badrate_nbins = 15
badrate_range = (0,15)
badrate_plot = [("badrate_kp","n_false","n_reco","Keypoint Net"),
                ("badrate_wc","n_false_wct","n_reco_wct","Keypoint Net+WC filter"),
                ("badrate_dl","n_false_dl","n_reco_dl","DLLEE Vertexer")]


for (varname,varnu,varbg,label) in badrate_plot:
    cbadrate[varname] = rt.TCanvas("cbadrate_"+varname,"Badrate vs. true Enu, %s"%(varname),800,400)    
    for (samplename,tree,color) in samples:
        hname = "h"+varname+"_"+samplename
        hbadrate[hname] = rt.TH1F(hname,";number of non-nu vertices",badrate_nbins, badrate_range[0], badrate_range[1])
        if samplename=="nu":
            tree.Draw("%s>>%s"%(varnu,hname))
        elif samplename=="bg":
            tree.Draw("%s>>%s"%(varbg,hname))
        hbadrate[hname].SetLineColor(color)
        s =  hbadrate[hname].Integral()
        hbadrate[hname].Scale( 1.0/s )
    tbadrate[samplename] = rt.TLegend(0.1,0.7,0.4,0.9)

    hname_nu = "h"+varname+"_nu"
    hname_bg = "h"+varname+"_bg"
    maxnu = hbadrate[hname_nu].GetMaximum()
    maxbg = hbadrate[hname_bg].GetMaximum()

    if maxnu:
        hbadrate[hname_nu].GetYaxis().SetRangeUser(0,1.5*maxnu)        
        hbadrate[hname_nu].Draw("hist")
    else:
        hbadrate[hname_bg].GetYaxis().SetRangeUser(0,1.5*maxbg)        
        hbadrate[hname_bg].Draw("hist")
            
    hname_nu = "h"+varname+"_nu"
    hname_bg = "h"+varname+"_bg"
    hbadrate[hname_nu].SetLineWidth(2)
    hbadrate[hname_bg].SetLineWidth(2)
    hbadrate[hname_nu].Draw("histsame")
    hbadrate[hname_bg].Draw("histsame")
    cbadrate[varname].Update()
        

# for (cutname,select_cut) in selection:
#     cbadrate[cutname].Clear()
#     c = cbadrate[cutname]
#     c.SetBottomMargin(0.15)
#     c.SetLeftMargin(0.1)

#     h = hbadrate["hbadrate_enu_kp_"+cutname]
#     h.SetTitle(";true neutrino energy (MeV); badrateiciency for %s events"%(cutname))
#     h.GetYaxis().SetRangeUser(0,1)
#     h.GetXaxis().SetTitleOffset(1.0)
#     h.GetYaxis().SetTitleOffset(0.75)    
#     h.GetXaxis().SetTitleSize(0.05)
#     h.GetYaxis().SetTitleSize(0.05)
#     h.GetXaxis().SetLabelSize(0.05)
#     h.GetYaxis().SetLabelSize(0.05)

    
#     hbadrate["hbadrate_enu_kp_"+cutname].Draw()
#     hbadrate["hbadrate_enu_wc_"+cutname].Draw("same")
#     hbadrate["hbadrate_enu_dl_"+cutname].Draw("same")

#     tbadrate[cutname].Draw()
    
#     cbadrate[cutname].Update()
#     cbadrate[cutname].SaveAs("cbadrate_%s.png"%(cutname))
    

raw_input()

