# L1MuGlobalNtuples

### To install a complete CMSSW area (including this package)

The code works inside a CMSSW environment. For Setting up the necessary CMSSW environment follow the recipe reported here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#Recipe

Then check out the code from this github repository and compile it on top of the CMSSW release:

```
cd $CMSSW_BASE/src/L1Trigger/
git clone https://github.com/pellicci/L1MuGlobalNtuples.git
cd $CMSSW_BASE/src/
scram b -j4 
```

### Code documentation

Please see the Twiki page https://twiki.cern.ch/twiki/bin/view/CMS/L1Phase2MuonNtuples for more documentation. 
