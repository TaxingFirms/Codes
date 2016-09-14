#!/bin/sh

# Use Intel Compiler
source /usr/global/intel/bin/compilervars.sh intel64

#$ -pe threaded 60

# Exports environment variables
#$ -V

# Execute from current directory
cd /home/dwills/Codes

# MyProgram
. /etc/profile.d/modules.sh
module load julia

julia -p60 runWelfareMaximizeHawk.jl >outJulia.txt
