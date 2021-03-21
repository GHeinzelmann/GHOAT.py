FILE=md02.nc


x=0
while [  $x -lt 10 ]; do
cd e0$x
if test -f "$FILE"; then
    echo "e0$x $FILE exists."
else
echo "e0$x $FILE does not exist."
qsub PBS-run
fi
cd ../
cd v0$x
if test -f "$FILE"; then
    echo "v0$x $FILE exists."
else
echo "v0$x $FILE does not exist."
qsub PBS-run
fi
cd ../
let x=x+1
done

if [ $x -ge 10 ]; then
while [  $x -lt 12 ]; do
cd e$x
if test -f "$FILE"; then
    echo "e$x $FILE exists."
else
echo "e$x $FILE does not exist."
qsub PBS-run
fi
cd ../
cd v$x
if test -f "$FILE"; then
    echo "v$x $FILE exists."
else
echo "v$x $FILE does not exist."
qsub PBS-run
fi
cd ../
let x=x+1
done
fi
