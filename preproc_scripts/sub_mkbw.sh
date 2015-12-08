export indir="/diazlab/BRAIN_data/tophat_out/${1}/"
export otdir_pref="/diazlab/BRAIN_data/bigwigs2/"
export otdir="$otdir_pref${1}/"
if [ ! -d "$otdir_pref${1}" ]
then
    mkdir "$otdir_pref${1}"
fi
nm="mkbw_${1}"
script="bigwigs_gen.pbs"
a=()
for f in `ls $indir`
do
    if [ -d $indir$f ] && [ -f $indir$f/accepted_hits.stats ]
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

qsub -e "${nm}.err" -N "${nm}" -t $rng -v indir,otdir,to_proc ${script}
