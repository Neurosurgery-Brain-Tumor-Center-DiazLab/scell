export indir="/diazlab/BRAIN_data/raw_data/${1}/"
export otdir_pref="/diazlab/BRAIN_data/tophat_out/${1}/"
if [ ! -d "/diazlab/BRAIN_data/tophat_out/${1}" ]
then
    mkdir "/diazlab/BRAIN_data/tophat_out/${1}"
fi
nm="align_${1}"
#human
export bowtie_idx_path="/diazlab/refs/Homo_sapiens/NCBI/build37.2/Sequence/Bowtie2Index/"
export transcriptome_idx_path="/diazlab/refs/Homo_sapiens/NCBI/build37.2/Annotation/Genes/genes_ercc"
export genome_name="genome_ercc"
script="align.pbs"

#chimp
#export bowtie_idx_path="/diazlab/refs/Pan_troglodytes/NCBI/build3.1/Sequence/Bowtie2Index/"
#export transcriptome_idx_path="/diazlab/refs/Pan_troglodytes/NCBI/build3.1/Annotation/Genes/genes_ercc"
#export genome_name="genome_ercc"
#script="align.pbs"

#baboon
#export bowtie_idx_path="/diazlab/refs/baboon/"
#export transcriptome_idx_path="/diazlab/refs/baboon/Annotation/papAnu2_genes_ercc"
#export genome_name="papAnu2_ercc"
#script="align.pbs"

#human_virus
#export bowtie_idx_path="/diazlab/refs/human_virus/"
#export genome_name="human_virus"
#script="align_virus.pbs"

a=()
for f in `ls $indir`
do
    if [ ! -d $otdir_pref$f ] || [ ! -f $otdir_pref$f/accepted_hits.stats ]
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

qsub -e "${nm}.err" -N "${nm}" -t $rng -v indir,otdir_pref,bowtie_idx_path,transcriptome_idx_path,genome_name,to_proc ${script}
