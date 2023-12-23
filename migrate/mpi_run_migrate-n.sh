#!
for r in */
	do
		cd $r
		echo $r
		date
		date > date.txt
			for m in */
			  do
				cd $m
				  echo $m
				  date
				  mpirun -np 110 hostfile ~/hosts ~/migrate-4.4.4/src/migrate-n-mpi parmfile
				  sleep 1
				cd ..
			  done
		cd ..
	done


