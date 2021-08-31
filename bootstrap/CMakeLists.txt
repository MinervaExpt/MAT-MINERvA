#A CMakeLists.txt that automatically sets up all dependencies for this package.

cmake_minimum_required(VERSION 3.2) #Requiring version 3.2

include(ExternalProject)

set(MY_GIT_TAG "main" CACHE STRING "Check out the same git tag for all projects") #TODO: Make this come from the command line/cache
set(GIT_CHECKOUT_PREFIX "https://github.com/" CACHE STRING "Select HTTPS or SSH checkout for external projects") #git@github.com: to use ssh checkout
set(SOURCE_CODE_DIR ${CMAKE_SOURCE_DIR}/../..)

#Stop the build if any source directories already exist.  Otherwise, CMake would blindly overwrite any changes!
#Contact aolivier@ur.rochester.edu for manual install instructions if this is stopping you from doing something.
if(EXISTS ${SOURCE_CODE_DIR}/MAT OR EXISTS ${SOURCE_CODE_DIR}/UnfoldUtils)
  message(FATAL_ERROR "One or more source code directories already exist.  I'm not going to download anything so I don't overwrite your hard work.")
endif()
if(EXISTS ${CMAKE_INSTALL_PREFIX}/MParamFiles OR EXISTS ${CMAKE_INSTALL_PREFIX}/MATFluxAndReweightFiles)
  message(FATAL_ERROR "Flux and/or reweight files are exist.  I'm not going to download anything so I don't overwrite your hard work.")
endif()

ExternalProject_Add(MAT
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
                    GIT_REPOSITORY ${GIT_CHECKOUT_PREFIX}MinervaExpt/UnfoldUtils.git
                    GIT_TAG ${MY_GIT_TAG}
                    SOURCE_DIR "${SOURCE_CODE_DIR}/UnfoldUtils"
                    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
                    STEP_TARGETS update
                    UPDATE_DISCONNECTED true
                    DEPENDS MAT MAT-MINERvA)

#MParamFiles doesn't need to be built or installed.  It's a bunch of text and .root files that programs need to be able to find.
#Nothing depends on MParamFiles to build, but it does need to be installed before much of libPlotUtils can be used.
ExternalProject_Add(MParamFiles
                    CVS_REPOSITORY minervacvs@cdcvs.fnal.gov:/cvs/mnvsoft
                    CVS_MODULE MParamFiles
                    ${CVS_DOWNLOAD_COMMAND}
                    SOURCE_DIR "${CMAKE_INSTALL_PREFIX}/etc/MParamFiles" #Check out directly into install prefix
                    CONFIGURE_COMMAND ""
                    BUILD_COMMAND ""
                    INSTALL_COMMAND ""
                    UPDATE_DISCONNECTED true
                    STEP_TARGETS update
                    )

#The flux files and many of our reweights now live in this package.
ExternalProject_Add(MATFluxAndReweightFiles
                    CVS_REPOSITORY minervacvs@cdcvs.fnal.gov:/cvs/mnvsoft
                    CVS_MODULE AnalysisFramework/Ana/MATFluxAndReweightFiles
                    ${CVS_DOWNLOAD_COMMAND}
                    SOURCE_DIR "${CMAKE_INSTALL_PREFIX}/etc/MATFluxAndReweightFiles" #Check out directly into install prefix
                    CONFIGURE_COMMAND ""
                    BUILD_COMMAND ""
                    INSTALL_COMMAND ""
                    UPDATE_DISCONNECTED true
                    STEP_TARGETS update
                    )

add_custom_target(update
                  COMMENT "Updating all projects from github/CVS")
add_dependencies(update MAT-MINERvA-update MAT-update UnfoldUtils-update MParamFiles-update MATFluxAndReweightFiles-update)

#Make a symbolic link to MATFluxAndReweightFiles because Andrew didn't think about
#naming it "data" before he created a bunch of CVS directories that only Rob can remove :(
#Nota Bene: This may not work well on Windows.  If you need to do that one day, read https://stackoverflow.com/questions/35765106/symbolic-links-cmake/41037224
#TODO: If /cvmfs exists, symlink to central MParamFiles and MATFluxAndReweightFiles instead of downloading them.
#      Not entirely sure symlinks work between /cvmfs, BlueArc, and CONDIR_INPUT_DIR to be fair.
install(CODE "EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E create_symlink ../etc/MATFluxAndReweightFiles ${CMAKE_INSTALL_PREFIX}/lib/data )" )

#Tell setup.sh about the installation directory of this package
configure_file(setup.sh.in setup.sh @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.sh DESTINATION bin)
install(FILES setupROOT6OnGPVMs.sh DESTINATION bin)
