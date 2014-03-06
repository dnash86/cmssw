for i in {0..48}
do
    value=$((11400 + $i))
    runTheMatrix.py --what upgrade -l $value --dryRun >> List.txt
done
