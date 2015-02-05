#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(voxelSValues)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleVoxelSValues exampleVoxelSValues.cc ${sources} ${headers})
target_link_libraries(exampleVoxelSValues ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build VoxelSValues. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEVoxelSValues_SCRIPTS
  exampleVoxelSValues.in
  exampleVoxelSValues.out
  init.mac
  init_vis.mac
  run1.mac
  run2.mac
  vis.mac
  1MeV_gamma_run.mac
  100keV_gamma_run.mac
  10keV_gamma_run.mac
  1MeV_electron_run.mac
  100keV_electron_run.mac
  10keV_electron_run.mac
  )

foreach(_script ${EXAMPLEVoxelSValues_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(VoxelSValues DEPENDS exampleVoxelSValues)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exampleVoxelSValues DESTINATION bin)


