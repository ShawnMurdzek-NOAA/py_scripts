#!/bin/sh

# Change directory to where this script is located
fullpath=$( readlink -f "${BASH_SOURCE[0]}" )
path=$( dirname "${fullpath}" )
cd $path

# Add ability to load modules and use slurm (Hera and Jet)
#source /etc/profile.d/modules.sh
#source /etc/profile.d/slurm.sh

# Add ability to load modules and use slurm (Orion)
source /etc/profile

module load contrib/0.1
module load rocoto
rocotorun -w perfect_ob_workflow_WINTER.xml -d perfect_ob_workflow_WINTER.db
