# dualCrysta
for a crystal like the one used in https://www.sciencedirect.com/science/article/pii/S0168900212005888

source ./g4env.sh
cmake -DGeant4_DIR=/cvmfs/geant4.cern.ch/geant4/10.5/x86_64-slc6-gcc63-opt/lib64/GEANT4-10.5.0

to do display
./DualCrystal -c template.cfg -u Xm 

to run without display
./DualCrystal -c template.cfg -m run.mac -o haha

make histograms in ROOT with CC_Analyzer.cc(1)
