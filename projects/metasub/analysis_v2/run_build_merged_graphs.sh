set -e

K=19

mem=500000
threads=8
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
seqdir=${basedir}/metasub/graphs/output_k${K}_cleaned
outdir=${basedir}/metasub/graphs/output_k${K}_cleaned_merged
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_clean/metagraph/build/metagraph

metadata=$(pwd)/../complete_metadata_extended.clean.v2.csv
chunksize=50
files=""
cnt=0
c=1
while IFS=',' read -r -a line || [[ -n "$line" ]]
do
    uuid=${line[1]}
    if [ "$uuid" == "uuid" ]
    then
        continue
    fi
    if [ ! -f ${seqdir}/${uuid}.clean.fasta.gz ]
    then
        echo "sequence file incomplete for $uuid"
        echo "$uuid complete"
        continue
    fi
    files="$files ${seqdir}/${uuid}.clean.fasta.gz"
    cnt=$(($cnt + 1))
    if [ "$cnt" == "$chunksize" ]
    then
        outgraph=${outdir}/merged_graph_k${K}.level1.chunk${c}.dbg
        outseq=${outdir}/merged_graph_k${K}.level1.chunk${c}
        logname=${outdir}/merged_graph_k${K}.level1.chunk${c}.lsf.log
        if [ ! -f ${outgraph} ]
        then
            cmd="/usr/bin/time -v $metagraph build --parallel $threads -k $K --mem-cap-gb 40 -o ${outgraph} $files && /usr/bin/time -v $metagraph assemble --parallel $threads --unitigs -v -o ${outseq} ${outgraph}"
        elif [ ! -f ${outseq}.fasta.gz ]
        then
            cmd="/usr/bin/time -v $metagraph assemble --unitigs -v -o ${outseq} ${outgraph}"
        else 
            cmd=""
        fi
        if [ ! -z "$cmd" ]
        then
            echo "$cmd" | bsub -J ms_merge -o ${logname} -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
        fi
        files=""
        cnt=0
        c=$(($c + 1))
    fi
done < <(cat $metadata)

if [ "$cnt" -gt 0 ]
then
    outgraph=${outdir}/merged_graph_k${K}.level1.chunk${c}.dbg
    outseq=${outdir}/merged_graph_k${K}.level1.chunk${c}
    logname=${outdir}/merged_graph_k${K}.level1.chunk${c}.lsf.log
    if [ ! -f ${outgraph} ]
    then
        cmd="/usr/bin/time -v $metagraph build --parallel $threads -k $K --mem-cap-gb 40 -o ${outgraph} $files && /usr/bin/time -v $metagraph assemble --parallel $threads --unitigs -v -o ${outseq} ${outgraph}"
    elif [ ! -f ${outseq}.fasta.gz ]
    then
        cmd="/usr/bin/time -v $metagraph assemble --unitigs -v -o ${outseq} ${outgraph}"
    else 
        cmd=""
    fi
    if [ ! -z "$cmd" ]
    then
        echo "$cmd" | bsub -J ms_merge -o ${logname} -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
    fi
fi
