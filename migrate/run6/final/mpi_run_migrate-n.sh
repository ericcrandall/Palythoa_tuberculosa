#!
for r in  rep*
	do
		cd $r
		echo $r
		date
		date > date.txt
			for m in */
			  do
				cd $m
				  echo $m > current.model.txt
				  date
				  mpirun -np 110 --hostfile ~/hosts ~/migrate-5.0.4/src/migrate-n-mpi parmfile
				  sleep 10
				cd ..
			  done
		cd ..
	done


