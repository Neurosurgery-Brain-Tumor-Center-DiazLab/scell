printf "cell"
IFS=$'\n' read -d '' -r -a lines < /diazlab/refs/human_virus/viral_ids.txt
for g in "${lines[@]}"
do
    printf "\t$g"
done
printf "\n" 
for f in Sample*
do
    if [ ! -f "$f/accepted_hits.sort.bam" ]; then
	continue
    fi
    printf $f
    for g in "${lines[@]}"
    do
	printf "\t`samtools view -c $f/accepted_hits.sort.bam $g`"
    done
    printf "\n"
done