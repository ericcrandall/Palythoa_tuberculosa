#!
#for s in */
#do
#       cd $s

        for m in */
        do
                cd $m
                echo $m
                screen -S $m
                exec mpirun -np 32 ~/migrate-4.4.4/src/migrate-n-mpi parmfile
                sleep 1
                cd ..
        done
#       cd ..
#done


