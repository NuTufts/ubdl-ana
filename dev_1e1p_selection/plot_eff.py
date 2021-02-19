import os,sys
import ROOT as rt
rt.gStyle.SetOptStat(0)

sample_names = [ "is1eVA","1e1p","all" ]
    
selcut_names  = [ "fv",
                  "vertexcand",
                  "minshower",
                  "nshowerprongs",
                  "ntrackprongs",        
                  "hadronic",
                  "showergap",
                  "vertexact",
                  "allreco",
                  "numcuts"]
selcut_colors = { "fv":rt.kBlack,
                  "vertexcand":rt.kBlue+1,
                  "minshower":rt.kCyan+2,
                  "nshowerprongs":rt.kGreen+3,
                  "ntrackprongs":rt.kOrange-3,
                  "hadronic":rt.kMagenta-2,
                  "showergap":rt.kOrange+4,
                  "vertexact":rt.kRed,
                  "allreco":rt.kBlack,
                  "numcuts":rt.kBlack}

rfile = rt.TFile(sys.argv[1])

canvas = {}
hists = []
tlen_v = []
for s in sample_names:
    canvas[s] = rt.TCanvas("c%s"%(s),s,1200,600)
    # load hists
    hs = {}
    for cut in selcut_names:
        h = rfile.Get( "hEff_%s_%scut"%(s,cut) )
        print h
        if h:
            h.SetLineColor( selcut_colors[cut] )
            if cut in ["vertexcand","vertexact","allreco"]:
                h.SetLineWidth(2)
            hs[cut] = h
    canvas[s].Draw()
    hs["fv"].SetTitle("%s Events"%(s))
    hs["fv"].Draw()
    tlen = rt.TLegend(0.8,0.2,0.95,0.8)
    for cut in selcut_names[0:-2]:
        tlen.AddEntry(hs[cut],cut,"L")
        hs[cut].Draw("same")
    tlen.Draw()
    canvas[s].SetGridx(1)
    canvas[s].SetGridy(1)     
    canvas[s].Update()
    canvas[s].SaveAs("c%s_eff_vs_enu.png"%(s))
    tlen_v.append(tlen)
    hists.append(hs)
print "enter to end"
raw_input()


