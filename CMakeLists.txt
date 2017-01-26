PROJECT(JetCounting)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  #set(VTK_LIBS_ORIG ${VTK_LIBRARIES})  # store the libs in an extra variable
  include(${VTK_USE_FILE})
else()
  find_package(ItkVtkGlue REQUIRED)
  include(${ItkVtkGlue_USE_FILE})
  set(Glue ItkVtkGlue)
endif()

# now restore the variables back = switch required to use ItkVtkGlue
#set(VTK_LIBRARIES ${VTK_LIBS_ORIG})

add_executable(JetCounting MACOSX_BUNDLE JetCounting.cxx JetCountingMeshToolbox.cxx)

target_link_libraries(JetCounting ${Glue} ${VTK_LIBRARIES} ${ITK_LIBRARIES})
