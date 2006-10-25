#!/bin/csh
cd $LS_SUBCWD
eval `scramv1 runtime -csh`
cmsRun l1CaloExample.cfg
#END
