from __future__ import print_function
import os,sys,argparse,json
from ctypes import c_int
from math import log

parser = argparse.ArgumentParser("Keypoint and Triplet truth")
parser.add_argument("input_file",type=str,help="file produced by 'dedx_...'")
parser.add_argument("-lm","--larmatch",type=str,default=None,required=True,help="output of larmatch deploy, contains larmatch hits")
parser.add_argument("-mc","--has-mc",action='store_true',default=False,help="If given, will try and plot MC tracks")
parser.add_argument("--draw-flash",action='store_true',default=False,help="If true, draw in-time flash PMT data [default: false]")
args = parser.parse_args()

import numpy as np
import ROOT as rt
from larlite import larlite
from larcv import larcv
from ublarcvapp import ublarcvapp
from larflow import larflow
larcv.SetPyUtil()

import plotly.graph_objects as go
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate


import lardly


color_by_options = ["larmatch","ssn-bg","ssn-track","ssn-shower","ssn-class","keypoint-nu","keypoint-track","keypoint-shower","flow-field"]
colorscale = "Viridis"
option_dict = []
for opt in color_by_options:
    option_dict.append( {"label":opt,"value":opt} )

# LOAD TREES
io = larlite.storage_manager(larlite.storage_manager.kREAD)
io.add_in_filename( args.larmatch )
io.open()

showerprof_rootfile = rt.TFile( args.input_file )
shrproftree = showerprof_rootfile.Get("showerprofile")

nentries = io.get_entries()
print("NENTRIES: ",nentries)
    

from larlite import larutil
dv = larutil.LArProperties.GetME().DriftVelocity()

def make_figures(entry,plotby="larmatch",minprob=0.0):

    print("making figures for entry={} plot-by={}".format(entry,plotby))
    global io

    io.go_to(entry)
    shrproftree.GetEntry(entry)
    
    ev_lfhits = io.get_data(larlite.data.kLArFlow3DHit,"larmatch")
    npoints = ev_lfhits.size()
    print("num larflow hits: ",npoints)

    hitindex = 0
    if plotby=="ssn-bg":
        hitindex = 10
    elif plotby=="ssn-track":
        hitindex = 11
    elif plotby=="ssn-shower":
        hitindex = 12
    elif plotby=="keypoint-nu":        
        hitindex = 13
    elif plotby=="keypoint-track":
        hitindex = 14
    elif plotby=="keypoint-shower":
        hitindex = 15

    detdata = lardly.DetectorOutline()

    traces_v = []

    if args.draw_flash:
        ev_flash = io.get_data(larlite.data.kOpFlash,"simpleFlashBeam")
        nflashes = 0
        for iflash in range(ev_flash.size()):
            flash = ev_flash.at(iflash)
            if flash.Time()>2.94 and flash.Time()<4.86:            
                flash_trace_v = lardly.data.visualize_larlite_opflash_3d( flash )
                traces_v += flash_trace_v
                nflashes += 1
                break
        if nflashes==0:
            traces_v += lardly.data.visualize_empty_opflash()
    

    if plotby=="larmatch":
        lfhits_v =  [ lardly.data.visualize_larlite_larflowhits( ev_lfhits, "larmatch", score_threshold=minprob, max_hits=30000) ]
        traces_v += lfhits_v + detdata.getlines()
    elif plotby in ["ssn-bg","ssn-track","ssn-shower","ssn-class","keypoint-nu","keypoint-track","keypoint-shower"]:
        xyz = np.zeros( (npoints,4 ) )
        ptsused = 0
        for ipt in xrange(npoints):
            hit = ev_lfhits.at(ipt)

            if hit.track_score<minprob:
                continue

            xyz[ptsused,0] = hit[0]
            xyz[ptsused,1] = hit[1]
            xyz[ptsused,2] = hit[2]
            if plotby=="ssn-class":
                idx = np.argmax( np.array( (hit[10],hit[11],hit[12]) ) )
                xyz[ptsused,3] = float(idx)/2.0
            else:
                xyz[ptsused,3] = hit[hitindex]
            #print(xyz[ptsused,3])
            ptsused += 1

        print("make hit data[",plotby,"] npts=",npoints," abovethreshold(plotted)=",ptsused)
        larflowhits = {
            "type":"scatter3d",
            "x": xyz[:ptsused,0],
            "y": xyz[:ptsused,1],
            "z": xyz[:ptsused,2],
            "mode":"markers",
            "name":plotby,
            "marker":{"color":xyz[:ptsused,3],"size":1,"opacity":0.8,"colorscale":'Viridis'},
        }
        #print(xyz[:ptsused,3])
        traces_v += [larflowhits]+detdata.getlines()
    elif plotby in ["flow-field"]:
        # must sample, if trying to draw triangles
        ptsused = 0
        index = np.arange(npoints)
        np.random.shuffle(index)     
        xyz = np.zeros( (npoints,7) )

        for ipt in index:

            if ptsused>=10000:
                break
            
            hit = ev_lfhits.at(ipt)
            if hit.track_score<minprob:
                continue

            paflen = np.sqrt( hit[16]*hit[16]+hit[17]*hit[17] + hit[18]*hit[18] )
            if paflen==0:
                paflen = 1.0
            
            xyz[ptsused,0] = hit[0]
            xyz[ptsused,1] = hit[1]
            xyz[ptsused,2] = hit[2]
            xyz[ptsused,3] = hit.track_score
            xyz[ptsused,4] = hit[16]/paflen
            xyz[ptsused,5] = hit[17]/paflen
            xyz[ptsused,6] = hit[18]/paflen
            ptsused += 1

        print("make hit data[",plotby,"] npts=",npoints," abovethreshold(plotted)=",ptsused)

        fig = go.Cone(
            x=xyz[:ptsused,0],
            y=xyz[:ptsused,1],
            z=xyz[:ptsused,2],
            u=xyz[:ptsused,4],
            v=xyz[:ptsused,5],
            w=xyz[:ptsused,6],
            colorscale='Blues',
            sizeref=10,
            sizemode="absolute")
        traces_v += [fig]+detdata.getlines()


    # Print profiles
    nprofs = shrproftree.start_v.size()
    for ishr in range(nprofs):
        profpid = shrproftree.pid_v[ishr]
        pos = np.zeros( (2,3) )
        for i in range(3):
            pos[0,i] = shrproftree.start_v[ishr][i]
            pos[1,i] = shrproftree.end_v[ishr][i]

        if profpid==11:
            profcolor = "rgb(255,0,255)"
        else:
            profcolor = "rgb(0,255,255)"
            
        profiletrace = {
            "type":"scatter3d",
            "x":pos[:,0],
            "y":pos[:,1],
            "z":pos[:,2],
            "mode":"lines",
            "name":"shr[%d]"%(ishr),
            "line":{"color":profcolor,"width":3}
            }
        traces_v.append( profiletrace )
        
    if args.has_mc:
        mctrack_v = lardly.data.visualize_larlite_event_mctrack( io.get_data(larlite.data.kMCTrack, "mcreco"), origin=1)
        traces_v += mctrack_v

        mcshower_v = lardly.data.visualize_larlite_event_mcshower( io.get_data(larlite.data.kMCShower, "mcreco"), return_dirplot=True )
        traces_v += mcshower_v[1:]
        
    return traces_v

def test():
    pass
    
app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)

server = app.server

axis_template = {
    "showbackground": True,
    #"backgroundcolor": "#141414", # black
    #"gridcolor": "rgba(255, 255, 255)",
    #"zerolinecolor": "rgba(255, 255, 255)",    
    "backgroundcolor": "rgba(100, 100, 100,0.5)",    
    "gridcolor": "rgb(50, 50, 50)",
    "zerolinecolor": "rgb(0, 0, 0)",
}

plot_layout = {
    "title": "",
    "height":800,
    "margin": {"t": 0, "b": 0, "l": 0, "r": 0},
    "font": {"size": 12, "color": "black"},
    "showlegend": False,
    #"plot_bgcolor": "#141414",
    #"paper_bgcolor": "#141414",
    "plot_bgcolor": "#ffffff",
    "paper_bgcolor": "#ffffff",
    "scene": {
        "xaxis": axis_template,
        "yaxis": axis_template,
        "zaxis": axis_template,
        "aspectratio": {"x": 1, "y": 1, "z": 3},
        "camera": {"eye": {"x": 1, "y": 1, "z": 1},
                   "up":dict(x=0, y=1, z=0)},
        "annotations": [],
    },
}

eventinput = dcc.Input(
    id="input_event",
    type="number",
    placeholder="Input Event")

minprob_input = dcc.Input(
    id="min_prob",
    type="text",
    placeholder="0.0")

plotopt = dcc.Dropdown(
    options=option_dict,
    value='truthmatch',
    id='plotbyopt',
    )
        

app.layout = html.Div( [
    html.Div( [ eventinput,
                minprob_input,
                plotopt,
                html.Button("Plot",id="plot")
    ] ),
    html.Hr(),
    html.Div( [
        dcc.Graph(
            id="det3d",
            figure={
                "data": [],
                "layout": plot_layout,
            },
            config={"editable": True, "scrollZoom": False},
        )],
              className="graph__container"),
    html.Div(id="out")
] )

                       
@app.callback(
    [Output("det3d","figure"),
     Output("out","children")],
    [Input("plot","n_clicks")],
    [State("input_event","value"),
     State("min_prob","value"),
     State("plotbyopt","value"),
     State("det3d","figure")],
    )
def cb_render(*vals):
    if vals[1] is None:
        print("Input event is none")
        raise PreventUpdate
    if vals[1]>=nentries or vals[1]<0:
        print("Input event is out of range")
        raise PreventUpdate
    if vals[3] is None:
        print("Plot-by option is None")
        raise PreventUpdate
    try:
        minprob = float(vals[2])
    except:
        print("min prob cannot be turned into float")
        raise PreventUpdate
    if minprob<0:
        minprob = 0.0
    if minprob>1.0:
        minprob = 1.0
    

    cluster_traces_v = make_figures(int(vals[1]),plotby=vals[3],minprob=minprob)
    #print(cluster_traces_v)
    vals[-1]["data"] = cluster_traces_v
    return vals[-1],"event requested: {} {} {}".format(vals[1],vals[2],vals[3])

if __name__ == "__main__":
    app.run_server(debug=True)
