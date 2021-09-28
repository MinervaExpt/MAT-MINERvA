#A CMakeLists.txt that automatically sets up all dependencies for this package.

cmake_minimum_required(VERSION 3.2) #Requiring version 3.2

include(ExternalProject)

set(MY_GIT_TAG "main" CACHE STRING "Check out the same git tag for all projects")
set(GIT_CHECKOUT_PREFIX "https://github.com/" CACHE STRING "Select HTTPS or SSH checkout for external projects") #git@github.com: to use ssh checkout
set(FLUX_FILE_DIR "/cvmfs/minerva.opensciencegrid.org/minerva/CentralizedFluxAndReweightFiles" CACHE STRING "Look for centralized flux files here")
set(SOURCE_CODE_DIR ${CMAKE_SOURCE_DIR}/../..)

#Try to get flux files from a central location
set(DOWNLOAD_FLUX_FILES FALSE)
if(EXISTS ${FLUX_FILE_DIR}/MATFluxAndReweightFiles AND EXISTS ${FLUX_FILE_DIR}/MParamFiles)
  message("Found centralized flux and reweight files at ${FLUX_FILE_DIR}.  I will not make another copy so that you save disk space.")
else()
  message("There are no flux files at ${FLUX_FILE_DIR}, so queueing a fresh download.  You'll get TWO copies of these: one to develop and another to go to the grid.  This will eat >4GB of disk space and be very slow!")
  set(DOWNLOAD_FLUX_FILES TRUE)
  set(FLUX_FILE_DIR ../etc)
endif()

#Stop the build if any source directories already exist.  Otherwise, CMake would blindly overwrite any changes!
#Contact aolivier@ur.rochester.edu for manual install instructions if this is stopping you from doing something.
set(DISABLE_GIT_DOWNLOAD)
set(DISABLE_CVS_DOWNLOAD)

if(NOT CMAKE_NOT_FIRST_RUN)
  if(EXISTS ${SOURCE_CODE_DIR}/MAT OR EXISTS ${SOURCE_CODE_DIR}/UnfoldUtils)
    message("One or more source code directories already exist.  I'm not going to download anything so I don't overwrite your hard work.")
    set(DISABLE_GIT_DOWNLOAD DOWNLOAD_COMMAND "")
  endif()
  if(EXISTS ${CMAKE_INSTALL_PREFIX}/etc/MParamFiles OR EXISTS ${CMAKE_INSTALL_PREFIX}/etc/MATFluxAndReweightFiles)
    message("Flux and/or reweight files already exist.  I'm not going to download anything so I don't overwrite your hard work.")
    set(DISABLE_CVS_DOWNLOAD DOWNLOAD_COMMAND "")
  endif()
endif()

ExternalProject_Add(MAT
                    "${DISABLE_GIT_DOWNLOAD}"
                    GIT_REPOSITORY ${GIT_CHECKOUT_PREFIX}MinervaExpt/MAT.git
                    GIT_TAG ${MY_GIT_TAG}
                    SOURCE_DIR "${SOURCE_CODE_DIR}/MAT"
                    STEP_TARGETS update
                    UPDATE_DISCONNECTED true
                    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE})

ExternalProject_Add(MAT-MINERvA
                    GIT_REPOSITORY ${GIT_CHECKOUT_PREFIX}MinervaExpt/MAT-MINERvA.git
                    GIT_TAG ${MY_GIT_TAG}
                    DOWNLOAD_COMMAND "" #Let git updates work, but not downloads.  Otherwise, you'd overwrite your MAT-MINERvA area every time you tried to install it.  In reality, I think the install just stops and fails.
                    SOURCE_DIR "${SOURCE_CODE_DIR}/MAT-MINERvA"
                    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
                    STEP_TARGETS update
                    UPDATE_DISCONNECTED true
                    DEPENDS MAT)

ExternalProject_Add(UnfoldUtils
                    "${DISABLE_GIT_DOWNLOAD}"
                    GIT_REPOSITORY ${GIT_CHECKOUT_PREFIX}MinervaExpt/UnfoldUtils.git
                    GIT_TAG ${MY_GIT_TAG}
                    SOURCE_DIR "${SOURCE_CODE_DIR}/UnfoldUtils"
                    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
                    STEP_TARGETS update
                    UPDATE_DISCONNECTED true
                    DEPENDS MAT MAT-MINERvA)

#MParamFiles doesn't need to be built or installed.  It's a bunch of text and .root files that programs need to be able to find.
#Nothing depends on MParamFiles to build, but it does need to be installed before much of libPlotUtils can be used.
if(DOWNLOAD_FLUX_FILES)
  ExternalProject_Add(MParamFiles
                      "${DISABLE_CVS_DOWNLOAD}"
                      CVS_REPOSITORY minervacvs@cdcvs.fnal.gov:/cvs/mnvsoft
                      CVS_MODULE MParamFiles
                      SOURCE_DIR "${SOURCE_CODE_DIR}/MParamFiles"
                      CONFIGURE_COMMAND ""
                      BUILD_COMMAND ""
                      INSTALL_COMMAND ""
                      UPDATE_DISCONNECTED true
                      STEP_TARGETS update
                      )
  install(DIRECTORY ${SOURCE_CODE_DIR}/MParamFiles DESTINATION etc)
  
  #The flux files and many of our reweights now live in this package.
  ExternalProject_Add(MATFluxAndReweightFiles
                      "${DISABLE_CVS_DOWNLOAD}"
                      CVS_REPOSITORY minervacvs@cdcvs.fnal.gov:/cvs/mnvsoft
                      CVS_MODULE AnalysisFramework/Ana/MATFluxAndReweightFiles
                      SOURCE_DIR "${SOURCE_CODE_DIR}/MATFluxAndReweightFiles"
                      CONFIGURE_COMMAND ""
                      BUILD_COMMAND ""
                      INSTALL_COMMAND ""
                      UPDATE_DISCONNECTED true
                      STEP_TARGETS update
                      )
  install(DIRECTORY ${SOURCE_CODE_DIR}/MATFluxAndReweightFiles DESTINATION etc)

  add_custom_target(update
                    COMMENT "Updating all projects from github/CVS")
  add_dependencies(update MAT-MINERvA-update MAT-update UnfoldUtils-update MParamFiles-update MATFluxAndReweightFiles-update)
else()
  add_custom_target(update
                    COMMENT "Updating all projects EXCEPT FLUX FILES from github")
  add_dependencies(update MAT-MINERvA-update MAT-update UnfoldUtils-update)
endif()

#Make a symbolic link to MATFluxAndReweightFiles because Andrew didn't think about
#naming it "data" before he created a bunch of CVS directories that only Rob can remove :(
#Nota Bene: This may not work well on Windows.  If you need to do that one day, read https://stackoverflow.com/questions/35765106/symbolic-links-cmake/41037224
#Also make a symbolic link for MParamFiles because I can't be bothered to figure out the right way to detect FLUX_FILE_DIR in setup.sh.in
install(DIRECTORY DESTINATION etc)

install(CODE "EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E create_symlink ${FLUX_FILE_DIR}/MATFluxAndReweightFiles ${CMAKE_INSTALL_PREFIX}/lib/data )" )
if(NOT EXISTS ${CMAKE_INSTALL_PREFIX}/etc/MParamFiles)
  install(CODE "EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E create_symlink ${FLUX_FILE_DIR}/MParamFiles ${CMAKE_INSTALL_PREFIX}/etc/MParamFiles )" )
endif()

#Tell setup.sh about the installation directory of this package
configure_file(setup.sh.in setup.sh @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.sh DESTINATION bin)
install(FILES setupROOT6OnGPVMs.sh DESTINATION bin)
