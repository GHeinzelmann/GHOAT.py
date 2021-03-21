FILE=md02.nc

x=0
while [  $x -lt 10 ]; do
cd t0$x
if test -f "$FILE"; then
    echo "t0$x $FILE exists."
else
echo "t0$x $FILE does not exist."
qsub PBS-run
fi
cd ../
cd l0$x
if test -f "$FILE"; then
    echo "l0$x $FILE exists."
else
echo "l0$x $FILE does not exist."
qsub PBS-run
fi
cd ../
cd a0$x
if test -f "$FILE"; then
    echo "a0$x $FILE exists."
else
echo "a0$x $FILE does not exist."
qsub PBS-run
fi
cd ../
cd r0$x
if test -f "$FILE"; then
    echo "r0$x $FILE exists."
else
echo "r0$x $FILE does not exist."
qsub PBS-run
fi
cd ../
let x=x+1
done

if [ $x -ge 10 ]; then
while [  $x -lt 16 ]; do
cd t$x
if test -f "$FILE"; then
    echo "t$x $FILE exists."
else
echo "t$x $FILE does not exist."
qsub PBS-run
fi
cd ../
cd l$x
if test -f "$FILE"; then
    echo "l$x $FILE exists."
else
echo "l$x $FILE does not exist."
qsub PBS-run
fi
cd ../
cd a$x
if test -f "$FILE"; then
    echo "a$x $FILE exists."
else
echo "a$x $FILE does not exist."
qsub PBS-run
fi
cd ../
cd r$x
if test -f "$FILE"; then
    echo "r$x $FILE exists."
else
echo "r$x $FILE does not exist."
qsub PBS-run
fi
cd ../
let x=x+1
done
fi

x=0
while [  $x -lt 10 ]; do
cd c0$x
if test -f "$FILE"; then
    echo "c0$x $FILE exists."
else
echo "c0$x $FILE does not exist."
qsub PBS-run
fi
cd ../
let x=x+1
done

if [ $x -ge 10 ]; then
while [  $x -lt 16 ]; do
cd c$x
if test -f "$FILE"; then
    echo "c$x $FILE exists."
else
echo "c$x $FILE does not exist."
qsub PBS-run
fi
cd ../
let x=x+1
done
fi



