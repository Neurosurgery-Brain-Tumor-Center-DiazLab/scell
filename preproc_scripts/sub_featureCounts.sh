export indr="/diazlab/BRAIN_data/tophat_out/${1}/"
export otdr="/diazlab/BRAIN_data/featureCounts_out/"
export fname="${1}_counts.txt"

nm="count_${1}"
script="run_featureCounts.pbs"
qsub -e "${nm}.err" -N "${nm}" -v indr,otdr,fname ${script}
