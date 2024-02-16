#!
for r in rep*/
	do
		cd $r
		echo $r
		date
		date > 2way_run.txt
			for m in  stepping.stone.splitsD.2way.NS.KO.D/ stepping.stone.splitsD.2way.SN.KO.D/
			  do
				echo $m > current_model.txt
				cd $m
				  echo $m
				  date
				  mpirun -np 110 --hostfile ~/hosts ~/migrate-5.0.4/src/migrate-n-mpi parmfile
				  sleep 10
				cd ..
			  done
		cd ..
	done


