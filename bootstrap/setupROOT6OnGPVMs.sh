#!/bin/bash
if [ -n "$ROOTSYS" ]; then
  echo "ROOT is already set up.  You must set up the MINERvA 101 tutorial from scratch.  Check your .bash_profile for \"root\""
fi

source /cvmfs/minerva.opensciencegrid.org/minerva/hep_hpc_products/setups
setup root v6_10_04d -q e14:prof
setup cmake v3_7_1
