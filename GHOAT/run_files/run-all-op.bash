cd rest
cd a-comp
sbatch SLURMM-run
cd ../
cd l-comp
sbatch SLURMM-run
cd ../
cd t-comp
sbatch SLURMM-run
cd ../
cd c-comp
sbatch SLURMM-run
cd ../
cd r-comp
sbatch SLURMM-run
cd ../
cd ../

cd sdr

cd e-comp
x=0
while [  $x -lt 10 ]; do
cd e0$x
sbatch SLURMM-run
cd ../
let x=x+1
done
if [ $x -ge 10 ]; then
while [  $x -lt 12 ]; do
cd e$x
sbatch SLURMM-run
cd ../
let x=x+1
done
fi
cd ../

cd v-comp
x=0
while [  $x -lt 10 ]; do
cd v0$x
sbatch SLURMM-run
cd ../
let x=x+1
done
if [ $x -ge 10 ]; then
while [  $x -lt 12 ]; do
cd v$x
sbatch SLURMM-run
cd ../
let x=x+1
done
fi
cd ../

cd ../

