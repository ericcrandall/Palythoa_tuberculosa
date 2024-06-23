#!
for m in stepping.stone.breaks/ island/
          do
          cd $m
          echo $m > ../current.model.txt
          echo date > date.txt
          mpirun -np 110 --hostfile ~/hosts ~/migrate-5.0.6/src/migrate-n-mpi parmfile
          sleep 10
          cd ..
      done
