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

# Dist. versus Enu
cdist = {}
hdist = {}
tdist = {}
dist_nbins = 50
dist_range = (0.0,10.0)
dist_plot = [("dist_all","min_dist_to_vtx","1==1",rt.kBlue,"All"),
             ("dist_kp","min_dist_to_vtx","1==1",rt.kBlue,"Keypoint Net"),
             #("dist_wc","min_dist_to_vtx_fit","1==1",rt.kCyan,"Keypoint Net+WC filter"),
             ("dist_dl","min_dist_to_vtx_dl","1==1",rt.kRed,"DLLEE Vertexer")]

hmax = {}

for (cutname,select_cut) in selection:
    cdist[cutname] = rt.TCanvas("cdist_"+cutname,"Dist vs. true Enu, %s"%(cutname),800,400)
    for (varname,var,varcut,color,label) in dist_plot:
        hname = "h"+varname+"_"+cutname
        hdist[hname] = rt.TH1F(hname,";distance to true vertex",dist_nbins, dist_range[0], dist_range[1])
        cut = "(" + fvcut + ") && (" + select_cut + ") && (" + varcut +")"
        print var,":: ",cut
        ev_ana.Draw("%s>>%s"%(var,hname),cut)
        hdist[hname].SetLineColor(color)
    #hallname = "h"+dist_plot[0][0]+"_"+cutname
    tdist[cutname] = rt.TLegend(0.6,0.7,0.9,0.9)
    histmax = 0.0
    for (varname,var,varcut,color,label) in dist_plot[1:]:
        hname = "h"+varname+"_"+cutname
        #print "Divide ",hname," by ",hallname
        s = hdist[hname].Integral()
        hdist[hname].Scale(1.0/s)
        if histmax<hdist[hname].GetMaximum():
            histmax = hdist[hname].GetMaximum()
        #hdist[hname].Divide( hdist[hallname] )
        tdist[cutname].AddEntry( hdist[hname], label, "L" )
    hmax[(cutname,select_cut)] = histmax
        

for (cutname,select_cut) in selection:
    cdist[cutname].Clear()
    c = cdist[cutname]
    c.SetBottomMargin(0.15)
    c.SetLeftMargin(0.15)

    h = hdist["hdist_kp_"+cutname]
    h.SetTitle(";distance to true vertex (cm); fraction of %s events"%(cutname))
    #h.GetYaxis().SetRangeUser(0,1)
    h.GetXaxis().SetTitleOffset(1.0)
    h.GetYaxis().SetTitleOffset(1.0)    
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetYaxis().SetLabelSize(0.05)

    hdist["hdist_kp_"+cutname].GetYaxis().SetRangeUser(0,1.2*hmax[(cutname,select_cut)])
    hdist["hdist_kp_"+cutname].Draw("hist")
    #hdist["hdist_wc_"+cutname].Draw("histsame")
    hdist["hdist_dl_"+cutname].Draw("histsame")

    tdist[cutname].Draw()
    
    cdist[cutname].Update()
    cdist[cutname].SaveAs("cdist_%s.png"%(cutname))
    

raw_input()
