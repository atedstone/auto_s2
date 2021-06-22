#!/bin/bash
export PATH=$PATH:/apps/bin
export MODULEPATH=/apps/modules:$MODULEPATH
/opt/sge/bin/linux-x64/qsub /home/tedstona/scripts/auto_s2/run_auto_s2.sh
