source /cvmfs/sft.cern.ch/lcg/contrib/gcc/6.3/x86_64-slc6/setup.sh
source /cvmfs/geant4.cern.ch/geant4/10.5/x86_64-slc6-gcc63-opt-MT/CMake-setup.sh
export CXX=`which g++`
export CC=`which gcc`
export PATH=$PATH:/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.11.1/Linux-x86_64/bin
source /cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc62-opt/setup.sh 
source $ROOTSYS/bin/thisroot.sh
export LIBRARY_PATH=/home/eno/dualReadout/fakelib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/home/eno/dualReadout/fakelib:$LD_LIBRARY_PATH





