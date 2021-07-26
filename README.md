# MAT-MINERvA
MINERvA-specific plugins to the MAT like FluxReweighter, MuonFunctions, and our CCInclusiveCuts.  Extends the [MAT](https://github.com/MinervaExpt/MAT)'s systematic uncertainty infrastructure with standard implementations of systematics across all MINERvA analyses.  Requires MAT, and best used with UnfoldUtils and MATFluxAndReweightFiles.

# Installation (Installs MAT too)
```
git checkout TODO
mkdir -p opt/build && cd opt/build
cmake ../../MAT-MINERvA/bootstrap -DCMAKE_INSTALL_PREFIX=`pwd`/.. -DCMAKE_BUILD_TYPE=Release
kinit #Fermilab Kerberos ticket for getting flux and reweight files
make install #-j 4
```

# How to Use MAT-MINERvA
MAT-MINERvA is compiled and installed using CMake 3 or later.

## Learn to use MINERvA's standard systematics
Follow the [tutorial](https://github.com/MinervaExpt/MINERvA-101-Cross-Section) to work through a simple inclusive analysis.  It's a good blueprint/fork point for building a new analysis.

## Use it in your own analysis
- CMake: Works automatically if you install in the same `opt` prefix.  Just create a separate `buildYourPackage` area under `opt`.
- Makefile: Your package can link against libMAT-MINERvA in `opt/lib`.  Headers to include are in `opt/include`.
- PyROOT: MAT-MINERvA automatically generates python bindings thanks to ROOT.  Use the `.rootlogon.C` in the interpreter instructions below, and `from ROOT import PlotUtils`.
- ROOT's CLING Interpreter: Install this `.rootlogon.C` in either your home area or your working area
```
{
  if( gSystem->Getenv("PLOTUTILSROOT") )
  {
    string newpath = string(gROOT->GetMacroPath()) + ":" + string("${PLOTUTILSROOT}/../bin" );
    gROOT->SetMacroPath( newpath.c_str() );
    gInterpreter->AddIncludePath( "${PLOTUTILSROOT}/../include" );
    gInterpreter->AddIncludePath( "${PLOTUTILSROOT}/../include/PlotUtils" );
    std::vector<std::string> packages = { "MAT", "MAT-MINERvA" };
    for(const std::string& package: packages)
    {
      gSystem->Load( gSystem->ExpandPathName(("$PLOTUTILSROOT/lib" + package + ".so").c_str()) );
    }
  }
}
```  

# Contributing
- Commit directly to the main branch for now, effectively the way CVS used to work.  If we get more breaking changes than we can handle, the respository maintainers will move to a pull-request-powered contribution workflow.
- Use `git merge` to fold in changes to your working area and other branches.  We're keeping the commit procedure as simple as possible at the cost of a little commit history bloat.
- For big changes that you want to collaborate on, create a new branch with "feature/" at the beginning of its name.
