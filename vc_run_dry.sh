#!/bin/csh

which $shell
echo "----------------------------------------"
echo "--------  Dry FMS Main script ----------"
echo "this script was modified by Lei Wang 11/12/2021"
echo "Starting at `date`"
echo "----------------------------------------"

# this is the latest version of the runscript

limit stacksize unlimited
module purge

module load intel openmpi netcdf-fortran/4.4.4
module list
setenv MALLOC_CHECK_ 0

set echo
rm *.txt

#####################################################
############  1/6: Configuration  ###################
#####################################################

# --- you may modify below---
# computer information
setenv account         eapsdept
setenv npes               16
setenv ntasks_per_node    16
setenv nodes              1
#setenv mem_per_node       MB
setenv run_time           1-00:00 #days-hours:minutes
# the number of CPUs
setenv exp_name           exp2_NCEPsymm


######################################################
#prepare the source script
sbatch --parsable --time=$run_time --ntasks-per-node=$ntasks_per_node --account=$account --job-name=$exp_name --ntasks=$npes --nodes=$nodes vc_dry_job.sbatch

######################################################
unset echo
echo "Finish at `date`"
