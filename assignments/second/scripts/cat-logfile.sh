#mpic++ $1 -o $2

for nodeNum in 1 2 3 4
do
    for pNum in 1 2 3 4 5 6 7 8
    do
	echo "this is p$nodeNum-$pNum --------------------" >> all.log
	cat p$nodeNum-$pNum.stdout >> all.log
    done
done

#qsub p$1-$2.pbs

#sleep 15 &

#echo "this is p$1-$2 --------------------" >> all.log
#cat p$1-$2.stdout >> all.log

#tail $3.stdout

