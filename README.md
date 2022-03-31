# MAT-MINERvA
MINERvA-specific plugins to the MAT like FluxReweighter, MuonFunctions, and our CCInclusiveCuts.  Extends the [MAT](https://github.com/MinervaExpt/MAT)'s systematic uncertainty infrastructure with standard implementations of systematics across all MINERvA analyses.  Requires MAT, and best used with UnfoldUtils and MATFluxAndReweightFiles.

## Installation (Installs MAT and UnfoldUtils too)
```
#If you're on a MINERvA GPVM hosted by Fermilab, you need to set up newer versions of ROOT and CMake first
source /cvmfs/minerva.opensciencegrid.org/minerva/hep_hpc_products/setups
setup root v6_10_04d -q e14:prof
setup cmake v3_7_1
kinit #Fermilab Kerberos ticket for getting flux and reweight files
unset SSH_ASKPASS  # this will stop it putting up an xwindow to ask for your token for https git access.

#Have you set up a github personal access token yet?  If you don't know what I'm talking about:
#First, follow these instructions: https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token
#When you get to setting permissions, tick only the "repo" box.
#Next, write down the key it generated.
#Then, run this before proceeding: git config --global credential.helper store #Saves your access token as a plaintext file in your home directory
#Later, when you're asked for a username and a password, use your github username and your access key.  This will save the access key forever.

git clone https://github.com/MinervaExpt/MAT-MINERvA.git #Use a github personal access token for checkout
#git clone git@github.com:MinervaExpt/MAT-MINERvA.git #Use an ssh key instead of a github personal access token
mkdir -p opt/build && cd opt/build
cmake ../../MAT-MINERvA/bootstrap -DCMAKE_INSTALL_PREFIX=`pwd`/.. -DCMAKE_BUILD_TYPE=Release

make install #-j 4
```

## How to Use MAT-MINERvA
MAT-MINERvA is compiled and installed using CMake 3 or later.  Once it has been installed, **`source opt/bin/setup.sh` to tell the operating system to use it**.

### Learn to use MINERvA's standard systematics
Follow the [tutorial](https://github.com/MinervaExpt/MINERvA-101-Cross-Section) to work through a simple inclusive analysis.  It's a good blueprint/fork point for building a new analysis.

### Use it in your own analysis
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
## Testing older python binding that allows `from PlotUtils import *`

Try this to get PlotUtils into your PYTHONPATH - this likely needs to go into opt/build/setup.sh

```export PYTHONPATH=<path to MAT>/MAT-MINERvA/python:<path to MAT>/MAT-MINERvA/python/PlotUtils ```


## More Expert Options
- Developing the flux files: During the installation procedure, when you run `cmake ...`, add this flag on the command line: `-DFLUX_FILE_DIR=none`.
- Keeping only 1 copy of the flux files without CVMFS: During the installation procedure, when you run `cmake ...`, add this flag on the command line: `-DFLUX_FILE_DIR=/path/to/flux/files`.  This is a good way to save several GB per MAT installation if you have multiple installations.

## Contributing
- Commit directly to the main branch for now, effectively the way CVS used to work.  If we get more breaking changes than we can handle, the respository maintainers will move to a pull-request-powered contribution workflow.
- Use `git merge` to fold in changes to your working area and other branches.  We're keeping the commit procedure as simple as possible at the cost of a little commit history bloat.
- For big changes that you want to collaborate on, create a new branch with "feature/" at the beginning of its name.
