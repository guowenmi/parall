#mpic++ $1 -o $2

for nodeNum in 1 2 3 4 5
do
    for pNum in 1 2 3 4 5 6 7 8
    do
      > p$nodeNum-$pNum.stdout
      > p$nodeNum-$pNum.pbs
	echo "#PBS -j oe -o p$nodeNum-$pNum.stdout -l nodes=$nodeNum:ppn=$pNum -q pp
	mpiexec -machinefile \$PBS_NODEFILE /home/s19026416/parall/assignments/second/sort 9999999999" >>  p$nodeNum-$pNum.pbs
 	qsub p$nodeNum-$pNum.pbs
    done
done

#qsub p$1-$2.pbs

#sleep 15 &

#echo "this is p$1-$2 --------------------" >> all.log
#cat p$1-$2.stdout >> all.log

#tail $3.stdout

