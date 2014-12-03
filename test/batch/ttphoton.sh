#!/bin/bash
export DPM_HOST=se1.particles.ipm.ac.ir
export DPNS_HOST=se1.particles.ipm.ac.ir
#lcg-ls "srm://se1.particles.ipm.ac.ir:8446/srm/managerv2?SFN=/dpm/particles.ipm.ac.ir/home/cms/store/user/goldouzian"
cd /home/goldouzian/CMSSW_5_3_7/src/myanalysis/Atq/test/
eval `scramv1 runtime -sh`
export LD_LIBRARY_PATH=/home/goldouzian:${LD_LIBRARY_PATH}
cmsRun ttphoton_cfg.py

