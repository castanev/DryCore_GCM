#!/bin/csh

#SBATCH -A eapsdept
#SBATCH --nodes=2 --ntasks=15
#SBATCH --time=1:30:00
#SBATCH --output=output_post_error.txt
#SBATCH --error=post_error.txt
#SBATCH --mail-user=castanev@purdue.edu #email this user
#SBATCH --mail-type=END #mail at end of run
#SBATCH --mail-type=BEGIN #mail at beginning of run
set echo
  python3 Post_processing.py
end