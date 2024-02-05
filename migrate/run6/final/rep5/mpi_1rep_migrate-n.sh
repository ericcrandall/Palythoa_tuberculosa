#!

			for m in   stepping.stone.KO.D/  
			  do
				cd $m
				  echo $m
				  date
				  mpirun -np 110 --hostfile ~/hosts ~/migrate-5.0.4/src/migrate-n-mpi parmfile
				  sleep 10
				cd ..
			  done



