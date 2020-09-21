CXX=g++

CXXFLAGS=-g -fPIC
CXXFLAGS+= `root-config --cflags`
CXXFLAGS+= `larlite-config --includes`
CXXFLAGS+= -I$(LARLITE_BASEDIR)/../
CXXFLAGS+= `larcv-config --includes`
CXXFLAGS+= `ublarcvapp-config --includes`
CXXFLAGS+= -I$(LARFLOW_INCDIR)

LDFLAGS=
LDFLAGS+=-L$(LARFLOW_LIBDIR) -lLArFlow_LArFlowConstants -lLArFlow_PrepFlowMatchData -lLArFlow_KeyPoints  -lLArFlow_Reco -lLArFlow_PrepFlowMatchData -lLArFlow_SCBoundary
LDFLAGS+=`ublarcvapp-config --libs` -lLArCVApp_Reco3D -lLArCVApp_LArliteHandler -lLArCVApp_ContourTools -lLArCVApp_UBWireTool -lLArCVApp_MCTools -lLArCVApp_ubdllee 
LDFLAGS+=`larcv-config --libs` -lLArCVCorePyUtil
LDFLAGS+=`larlite-config --libs`
LDFLAGS+=`root-config --libs`


all: charge_on_path

%: %.cxx
	@echo "<< compile $@ >>"
	$(CXX) $(CXXFLAGS)  $< -o $@ $(LDFLAGS)

clean:
	@rm dedx_1m1p_ana
