#!/bin/bash

docker pull anu9109/capp-seq

# Run bbduk to trim adapters
echo -e "\e[0;36mRunning bbduk trimming on bam \e[0m"
docker run -v /data:/data -v /tmp:/tmp  -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c './opt/software/bbmap/bbduk.sh in=$infastq1 in2=$infastq2 out=$fastq1 out2=$fastq2 \ 
ref=/opt/software/bbmap/resources/truseq.fa.gz stats=$stats bhist=$bhist qhist=$qhist aqhist=$aqhist lhist=$lhist gchist=$gchist gcbins=auto threads=$cpu minlen=25 qtrim=rl trimq=10 ktrim=r k=25 \ 
mink=11 hdist=1 overwrite=true tbo=t tpe=t'

# Run FastQC
echo -e "\e[0;36mRunning FastQC \e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'fastqc -o $fastqc_dir $fastq1 $fastq2 -t $cpu'

# Alignment
echo -e "\e[0;36mCreating BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -v /dev:/dev-i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'bwa mem -M -t $cpu -R $RGR $bwa_index $fastq1 $fastq2 | \
samblaster --splitterFile >(samtools view -S -u /dev/stdin | sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $tmpsortbam /dev/stdin) \
--discordantFile >(samtools view -S -u /dev/stdin | sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $tmpsortbam /dev/stdin) | \
samtools view -S -u /dev/stdin | \
sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $tmpsortbam /dev/stdin'

# Index bam file
echo -e "\e[0;36mIndexing BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'samtools index $tmpsortbam'

# Mark duplicates
echo -e "\e[0;36mMarking duplicates in BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'sambamba markdup -t $cpu --overflow-list-size 1000000 --hash-table-size 1000000 \
$tmpsortbam $tmpdedupbam'

# Index de-duplicated bam file (Note: Add clean up step)
echo -e "\e[0;36mIndexing de-duplicated BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'samtools index $tmpdedupbam'

# Indel realignment with Abra (Note: Add clean up step)
echo -e "\e[0;36mPerforming indel realignment on BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'java -Xmx16G -jar /opt/software/abra.jar --in $tmpdedupbam --out $tmprealnbam \
--ref $bwa_index --targets $regions --threads $cpu --working $tmpabra > $tmpabralog 2>&1'

# Sort and index re-aligned bam file (Note: Add clean up step)
echo -e "\e[0;36mSorting re-aligned BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $outbam $tmprealnbam'
echo -e "\e[0;36mIndexing re-aligned BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'samtools index $outbam'

# Run Qualimap
echo -e "\e[0;36mRunning Qualimap on BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'qualimap bamqc --java-mem-size=10G -gd HUMAN -sd -gff $regions -bam $outbam -outdir $qm_dir --outformat HTML'

# Run VarDict
echo -e "\e[0;36mCalling variants with VarDict\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'VarDict -th $cpu -Q 10 -q 20 -G $bwa_index -f 0.01 -t -N $sample -b $outbam -c 1 -S 2 -E 3 -g 4 $regions | /opt/software/VarDictJava/VarDict/teststrandbias.R  | /opt/software/VarDictJava/VarDict/var2vcf_valid.pl -N $sample -f 0.01 -P 0.9 -m 4.25 > $outvardict'

echo -e "\e[0;36mCreating compressed VCF\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'bgzip -f $outvardict'
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'tabix -f -p vcf $outvardictgz'

# Cleaning temp folder
echo -e "\e[0;36mCleaning /tmp space used by Abra \e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'rm -r $tmpabra; rm /tmp/$sample*'

echo -e "\e[0;36mDone! \e[0m"
exit 0


