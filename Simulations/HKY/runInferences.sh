files=inf*.xml
c=1
for f in $files
do
    if [ "$c" -eq "1" ]; then
        jobid=$(sbatch -J runsims --mem=8gb -o ${f/.xml/.out} ~/jar/random Recombination $c $f)
        echo $jobid
    else
        name="$(cut -d' ' -f4 <<<$jobid)"
        jobid=$(sbatch --dependency=afterany:$name -c 1 -J runsims --mem=8gb -o ${f/.xml/.out} ~/jar/random Recombination $c $f)
    fi
    
    let c=$c+1
    echo $c
    
    if [ "$c" -gt "15" ]; then
        c=1
    fi

done
