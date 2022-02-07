#!
for m in */
  do
    cd $m
      echo $m
      mpirun -np 32 ~/migrate-4.4.4/src/migrate-n-mpi parmfile
      sleep 1
    cd ..
  done


