#mpic++ $1 -o $2

for nodeNum in 1 2 3 4 5
do
    for pNum in 1 2 3 4 5 6 7 8
    do
      > p$nodeNum-$pNum.stdout
      > p$nodeNum-$pNum.pbs
	echo "#PBS -j oe -o p$nodeNum-$pNum.stdout -l nodes=$nodeNum:ppn=$pNum -q pp
	mpiexec -machinefile \$PBS_NODEFILE /home/s19026416/parall/assignments/first/first 500000000" >>  p$nodeNum-$pNum.pbs
 	qsub p$nodeNum-$pNum.pbs
    done
done


