files=sim*.xml
for f in $files
do
	sbatch -c 1 -J anoph -N 1 -o ${f/.xml/.out} --time=1:00:00 ~/jar/resume Coevo  $f
done
