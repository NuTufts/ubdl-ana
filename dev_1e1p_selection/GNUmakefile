CXX=g++

CXXFLAGS=-g -fPIC
CXXFLAGS+= `root-config --cflags`
CXXFLAGS+= `larlite-config --includes`
CXXFLAGS+= -I$(LARLITE_BASEDIR)/../
CXXFLAGS+= `larcv-config --includes`
CXXFLAGS+= `ublarcvapp-config --includes`
CXXFLAGS+= -I$(LARFLOW_INCDIR)

LDFLAGS=-L$(LARFLOW_LIBDIR) -lLArFlow_LArFlowConstants -lLArFlow_PrepFlowMatchData -lLArFlow_KeyPoints  -lLArFlow_Reco -lLArFlow_PrepFlowMatchData -lLArFlow_SCBoundary 
LDFLAGS+=`ublarcvapp-config --libs` -lLArCVApp_MCTools -lLArCVApp_ubdllee -lLArCVApp_UBWireTool -lLArCVApp_LArliteHandler -lLArCVApp_ContourTools -lLArCVApp_UBWireTool
LDFLAGS+=`larcv-config --libs` -lLArCVCorePyUtil
LDFLAGS+=`larlite-config --libs`
LDFLAGS+=`root-config --libs`

EXECS = make_1e1p_sel_plots make_finalstate_plots #make_shower_reco_acc_plots
LIBS = libSelection.so

all: $(LIBS) $(EXECS)

%: %.cxx
	@echo "<< compile $@ >>"
	$(CXX) $(CXXFLAGS)  $< -o $@ $(LDFLAGS)

libSelection.so: NueSelection.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^

%.o: %.cxx %.h
	@echo "<< compile $@ >>"
	$(CXX) $(CXXFLAGS)  -c $< -o $@

clean:
	@rm $(EXECS)
