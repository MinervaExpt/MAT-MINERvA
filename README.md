# MAT-MINERvA
MINERvA-specific plugins to the MAT like FluxReweighter, MuonFunctions, and our CCInclusiveCuts.  Extends the [MAT](https://github.com/MinervaExpt/MAT)'s systematic uncertainty infrastructure with standard implementations of systematics across all MINERvA analyses.

# Installation (Installs MAT too)
- git checkout TODO
- mkdir -p opt/build && cd opt/build
- cmake ../../MAT-MINERvA/bootstrap
- kinit #Fermilab Kerberos ticket for getting flux and reweight files
- make install #-j 4
