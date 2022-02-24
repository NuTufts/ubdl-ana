add_executable( ana_kpreco ana_kpreco.cxx )

set(LARCV_LIBS_USED LArCVCoreBase LArCVCoreDataFormat LArCVCoreProcessor  )
set(UBLARCVAPP_LIBS_USED LArCVApp_ContourTools LArCVApp_dbscan LArCVApp_UBImageMod LArCVApp_UBWireTool LArCVApp_ubdllee LArCVApp_MCTools LArCVApp_LLCVProcessor )
set(LARLITE_LIBS_USED larlite_DataFormat larlite_LArUtil larlite_Base )

target_include_directories( ana_kpreco PUBLIC ROOT )
target_include_directories( ana_kpreco PUBLIC ${LARFLOW_INCLUDE_DIR} )
target_include_directories( ana_kpreco PUBLIC ${LARLITE_INCLUDE_DIR} )
target_include_directories( ana_kpreco PUBLIC ${LARCV_INCLUDE_DIR} )
target_link_libraries( ana_kpreco PUBLIC LArFlow_Reco LArFlow_SCBoundary )
target_link_libraries( ana_kpreco PUBLIC ${UBLARCVAPP_LIBS_USED} )
target_link_libraries( ana_kpreco PUBLIC ${LARLITE_LIBS_USED} )
target_link_libraries( ana_kpreco PUBLIC ${LARCV_LIBS_USED} )
target_link_libraries( ana_kpreco PUBLIC ${ROOT_LIBRARIES} )
target_link_libraries( ana_kpreco PUBLIC ${OpenCV_LIBRARIES} )



# ublarcvapp libs used                                                                                                                                                                        
