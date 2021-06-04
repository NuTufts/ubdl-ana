import os,sys
import ROOT as rt
rt.gStyle.SetOptStat(0)
rt.gStyle.SetPadBottomMargin(0.15)
rt.gStyle.SetPadLeftMargin(0.15)
rt.gStyle.SetPadTopMargin(0.10)
rt.gStyle.SetPadRightMargin(0.10)

infile = sys.argv[1]

f = rt.TFile(infile)
tree = f.Get("showerbilinear")

pid_v = ["electron","gamma"]
pid_cut = {"electron":"abs(true_match_pdg)==11",
           "gamma":"abs(true_match_pdg)==22"}
goodreco_cut = "true_dir_cos>0.90 && abs(true_vertex_err_dist)<3.0"
hists = {}

## ==============================================================
## dq/dx e vs. gamma
## ==============================================================
cdqdx = rt.TCanvas("cdqdx","dqdx",600,400)
hmax = None
for pid in pid_v:
    hname = "hdqdx_%s"%(pid)
    h = rt.TH1D(hname,";dq/dx over within first 3 cm (pixsum/cm);counts",50,0,1000)
    tree.Draw("best_pixsum_dqdx>>%s"%(hname),goodreco_cut + " && " + pid_cut[pid])
    if hmax is None or h.GetMaximum()>hmax.GetMaximum():
        hmax = h
    if pid in ("gamma"):
        h.SetLineColor(rt.kRed)
    h.SetLineWidth(2)
    hists[hname] = h

hmax.Draw()
hmax.GetYaxis().SetLabelSize(0.05)
hmax.GetXaxis().SetLabelSize(0.05)
hmax.GetXaxis().SetTitleSize(0.05)
hmax.GetXaxis().SetTitleOffset(1.2)
hmax.GetYaxis().SetTitleSize(0.05)
for pid in pid_v:
    hists["hdqdx_%s"%(pid)].Draw("samehist")
#cdqdx.SetPadBottomMargin(0.15)
#cdqdx.SetPadLeftMargin(0.15)
#cdqdx.SetPadTopMargin(0.05)
#cdqdx.SetPadRightMargin(0.05)
cdqdx.Update()
cdqdx.SaveAs("cdqdx.png")

## ==============================================================
## truth-reco matching variables
## ==============================================================
ctruth = rt.TCanvas("ctruth","ctruth",1200,400)
ctruth.Divide(2,1)
ctruth.cd(1)
hmax = 0.0
for pid in pid_v:
    hname = "htruth_"+pid
    h = rt.TH2D(hname,"best match to true %s;distance from true vertex projected on shower trunk (cm); cos of true and reco shower trunk"%(pid),
                50,-5,10,51,-1.01,1.01)
    tree.Draw("true_dir_cos:true_vertex_err_dist>>%s"%(hname),pid_cut[pid])
    if hmax<h.GetMaximum():
        hmax = h.GetMaximum()
    h.GetYaxis().SetLabelSize(0.05)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetTitleSize(0.05)        
    hists[hname] = h
ctruth.cd(1).SetGridy(1)
ctruth.cd(1).SetGridx(1)
hists["htruth_electron"].GetZaxis().SetRangeUser(0,hmax)
hists["htruth_electron"].Draw("colz")
ctruth.cd(2)
ctruth.cd(2).SetGridy(1)
ctruth.cd(2).SetGridx(1)
hists["htruth_gamma"].GetZaxis().SetRangeUser(0,hmax)
hists["htruth_gamma"].Draw("colz")
ctruth.Update()
ctruth.SaveAs("ctruth.png")

## ==============================================================
## error vs. dq/dx value
## ==============================================================
cerr_v_dqdx = rt.TCanvas("cerr_v_dqdx","cerr_v_dqdx",1200,400)
cerr_v_dqdx.Divide(2,1)
cerr_v_dqdx.cd(1)
hmax = 0.0
for pid in pid_v:
    hname = "herr_v_dqdx_"+pid
    h = rt.TH2D(hname,"best match to true %s;dq/dx within first 3 cm (pixsum/cm); fractional dq/dx error, (reco-true)/true"%(pid),
                25,0,1000,50,-1.01,2.01)
    tree.Draw("(best_pixsum_dqdx-true_matched_best_pixsum)/true_matched_best_pixsum:best_pixsum_dqdx>>%s"%(hname),
              goodreco_cut + " && " + pid_cut[pid] + " && best_pixsum_dqdx>0 && true_matched_best_pixsum>0")
    if hmax<h.GetMaximum():
        hmax = h.GetMaximum()
    h.GetYaxis().SetLabelSize(0.05)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetTitleSize(0.05)        
    hists[hname] = h
cerr_v_dqdx.cd(1).SetGridy(1)
cerr_v_dqdx.cd(1).SetGridx(1)
hists["herr_v_dqdx_electron"].GetZaxis().SetRangeUser(0,hmax)
hists["herr_v_dqdx_electron"].Draw("colz")
cerr_v_dqdx.cd(2)
cerr_v_dqdx.cd(2).SetGridy(1)
cerr_v_dqdx.cd(2).SetGridx(1)
hists["herr_v_dqdx_gamma"].GetZaxis().SetRangeUser(0,hmax)
hists["herr_v_dqdx_gamma"].Draw("colz")
cerr_v_dqdx.Update()
cerr_v_dqdx.SaveAs("cerr_v_dqdx.png")

## ==============================================================
## cos of other particle vs. dq/dx value
## ==============================================================
cdqdx_v_pcos = rt.TCanvas("cdqdx_v_pcos","cdqdx_v_pcos",1200,400)
cdqdx_v_pcos.Divide(2,1)
cdqdx_v_pcos.cd(1)
hmax = 0.0
for pid in pid_v:
    hname = "hdqdx_v_pcos_"+pid
    h = rt.TH2D(hname,"best match to true %s;dq/dx within first 3 cm (pixsum/cm);cos between shower and another particle"%(pid),
                25,0,1000,50,-1.01,1.01)
    tree.Draw("true_max_primary_cos:best_pixsum_dqdx>>%s"%(hname),
              goodreco_cut + " && " + pid_cut[pid])
    if hmax<h.GetMaximum():
        hmax = h.GetMaximum()
    h.GetYaxis().SetLabelSize(0.05)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetTitleSize(0.05)        
    hists[hname] = h
cdqdx_v_pcos.cd(1).SetGridy(1)
cdqdx_v_pcos.cd(1).SetGridx(1)
hists["hdqdx_v_pcos_electron"].GetZaxis().SetRangeUser(0,hmax)
hists["hdqdx_v_pcos_electron"].Draw("colz")
cdqdx_v_pcos.cd(2)
cdqdx_v_pcos.cd(2).SetGridy(1)
cdqdx_v_pcos.cd(2).SetGridx(1)
hists["hdqdx_v_pcos_gamma"].GetZaxis().SetRangeUser(0,hmax)
hists["hdqdx_v_pcos_gamma"].Draw("colz")
cdqdx_v_pcos.SaveAs("cdqdx_v_pcos.png")

raw_input()

