#!/bin/bash
#$ -N auto_s2
#$ -S /bin/sh
#$ -pe smp 2
#$ -q new.q
#$ -o /flash/tedstona/auto_s2/s2.out
#$ -e /flash/tedstona/auto_s2/s2.err
#$ -m eas
#$ -M andrew.tedstone@unifr.ch

source /home/tedstona/scripts/auto_s2/setup_sentinel2.sh

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/tedstona/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/tedstona/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/tedstona/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/tedstona/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate auto_s2

python /home/tedstona/scripts/auto_s2/auto_sentinel2.py

