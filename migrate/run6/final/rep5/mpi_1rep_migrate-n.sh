#!

			for m in  stepping.stone.splits.d/  stepping.stone.splitsD.KOsplitd/ stepping.stone.oneway/  stepping.stone.splitsD/
			  do
				cd $m
				  echo $m
				  date
				  mpirun -np 110 --hostfile ~/hosts ~/migrate-5.0.4/src/migrate-n-mpi parmfile
				  sleep 10
				cd ..
			  done



