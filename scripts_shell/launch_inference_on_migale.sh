#!/bin/bash

## For qsub using SGE on migale with 3 threads reserved
CMD="qsub -S /bin/bash -q long.q -cwd -V -M paul.bastide@agroparistech.fr -m ae -pe thread 3"

rscript="../simulations_study/dummy.R"
for i in {1..10}
do
    SCRIPT="inference_script.sh"
    rscript_i=${rscript%.R}_${i}.R
    echo "#!/bin/bash/" > $SCRIPT
    echo "R --vanilla < ${rscript_i}" >> $SCRIPT

    OUT="inferences_script_${i}.out"
    ERR="inferences_script_${i}.err"
    echo "$CMD $SCRIPT -o $OUT -e $ERR"
done
