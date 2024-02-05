#!
for r in */
	do
		cd $r
		echo $r
		date
		date > 2000_bins_rerun.txt
			for m in  panmixia/ island/ NWHI_MHI/ stepping.stone.breaks/ stepping.stone.breaks.KO/
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


