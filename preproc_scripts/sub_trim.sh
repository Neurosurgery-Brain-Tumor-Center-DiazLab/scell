export indir="/diazlab/BRAIN_data/raw_data/${1}/"
a=()
for f in `ls $indir`
do
    if [ ! -d $indir$f/trimmed ]
    then
	a+=($f)
    fi
done
to_proc=''
rng="0-$((${#a[@]}-1))%50"
for f in ${a[@]};
do
    to_proc=$f:$to_proc
done
export to_proc
nm="trim"
script="/diazlab/BRAIN_data/scripts/trim.pbs"
qsub -e "${nm}-${1}.err" -N "${nm}-${1}" -t $rng -v indir,to_proc ${script}
echo "qsub -e "${nm}-${1}.err" -N "${nm}-${1}" -t $rng -v indir,to_proc ${script}"