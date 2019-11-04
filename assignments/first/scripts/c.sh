mpic++ $1 -o $2
qsub $3.pbs
tail $3.stdout

