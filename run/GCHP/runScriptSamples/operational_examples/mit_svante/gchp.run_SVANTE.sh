#!/bin/bash

#SBATCH --ntasks-per-node=48
#SBATCH -N 6
#SBATCH -c 1
#SBATCH -t 1-12:00
#SBATCH -p hdr
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH -J GCHP

# Define log name to include simulation start date
start_str=$(sed 's/ /_/g' cap_restart)
log=gchp.${start_str:0:13}z.log
echo "Writing output to ${log}"

# Update config files, set restart symlink, load run env, and do sanity checks
source setCommonRunSettings.sh > ${log}
source setRestartLink.sh >> ${log}
source gchp.env >> ${log}
source checkRunSettings.sh >> ${log}

# Run GCHP and evenly distribute tasks across nodes
NX=$( grep NX GCHP.rc | awk '{print $2}' )
NY=$( grep NY GCHP.rc | awk '{print $2}' )
coreCount=$(( ${NX} * ${NY} ))
time mpirun -n ${coreCount} ./gchp >> ${log}

# Rename mid-run checkpoint files, if any. Discard file if time corresponds
# to run start time since duplicate with initial restart file.
chkpnts=$(ls Restarts)
for chkpnt in ${chkpnts}
do
    if [[ "$chkpnt" == *"gcchem_internal_checkpoint."* ]]; then
       chkpnt_time=${chkpnt:27:13}
       if [[ "${chkpnt_time}" = "${start_str:0:13}" ]]; then
          rm ./Restarts/${chkpnt}
       else
          new_chkpnt=./Restarts/GEOSChem.Restart.${chkpnt_time}z.c${N}.nc4
          mv ./Restarts/${chkpnt} ${new_chkpnt}
       fi
    fi
done

# Rename restart file and update restart symlink if new start time ok
new_start_str=$(sed 's/ /_/g' cap_restart)
if [[ "${new_start_str}" = "${start_str}" || "${new_start_str}" = "" ]]; then
   echo "ERROR: cap_restart either did not change or is empty."
   rm -f Restarts/gcchem_internal_checkpoint
   exit 1
else
    N=$(grep "CS_RES=" setCommonRunSettings.sh | cut -c 8- | xargs )    
    mv Restarts/gcchem_internal_checkpoint Restarts/GEOSChem.Restart.${new_start_str:0:13}z.c${N}.nc4
    source setRestartLink.sh
fi

exit 0
