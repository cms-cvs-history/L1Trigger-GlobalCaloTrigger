#!/bin/bash
cd $LS_SUBCWD
eval `scramv1 runtime -sh`
cmsRun l1CaloExample.cfg
#END
