#$ -S /bin/sh

WORK='/storage_scratch/users/adriana.pitea/GenoCNA'

R --no-save < $WORK/ > $WORK/$SGE_TASK_ID.out
