#!
for sim in C*/
  do
    echo $sim
    cd $sim
    ~/installers/IBDsim/IBDSim 
    cd ..
  done
