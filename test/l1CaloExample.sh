#!/bin/bash
cd $LS_SUBCWD
eval `scramv1 runtime -sh`
export STAGE_SVCCLASS=default
cmsRun l1CaloExample.cfg
#END
