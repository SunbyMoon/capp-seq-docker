#!/bin/bash

docker pull anu9109/capp-seq

# Alignment
echo -e "\e[0;36mCreating tumor BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /dev:/dev -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'bwa mem -M -t $cpu -R $tRGR $bwa_index  $tfastq1 $tfastq2 | \
#samblaster --splitterFile >(samtools view -S -u /dev/stdin | sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $ttmpsortbam /dev/stdin) \
#--discordantFile >(samtools view -S -u /dev/stdin | sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $ttmpsortbam /dev/stdin) | \
#samtools view -S -u /dev/stdin | \
#sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $ttmpsortbam /dev/stdin'
echo -e "\e[0;36mCreating normal BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /dev:/dev -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'bwa mem -M -t $cpu -R $gRGR $bwa_index  $gfastq1 $gfastq2 | \
#samblaster --splitterFile >(samtools view -S -u /dev/stdin | sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $gtmpsortbam /dev/stdin) \
#--discordantFile >(samtools view -S -u /dev/stdin | sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $gtmpsortbam /dev/stdin) | \
#samtools view -S -u /dev/stdin | \
#sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $gtmpsortbam /dev/stdin'

# Index bam file
echo -e "\e[0;36mIndexing tumor BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'samtools index $ttmpsortbam'
echo -e "\e[0;36mIndexing normal BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'samtools index $gtmpsortbam'

# Mark duplicates
echo -e "\e[0;36mMarking duplicates in tumor BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'sambamba markdup -t $cpu --overflow-list-size 1000000 --hash-table-size 1000000 $ttmpsortbam $ttmpdedupbam'
echo -e "\e[0;36mMarking duplicates in normal BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'sambamba markdup -t $cpu --overflow-list-size 1000000 --hash-table-size 1000000 $gtmpsortbam $gtmpdedupbam'

# Index de-duplicated bam file (Note: Add clean up step)
echo -e "\e[0;36mIndexing de-duplicated tumor BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'samtools index $ttmpdedupbam'
echo -e "\e[0;36mIndexing de-duplicated normal BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'samtools index $gtmpdedupbam'

# Indel realignment with Abra (Note: Add clean up step)
echo -e "\e[0;36mPerforming indel realignment on tumor BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'java -Xmx16G -jar /opt/software/abra.jar --in $ttmpdedupbam --out $ttmprealnbam \
#--ref $bwa_index --targets $regions --threads $cpu --working $ttmpabra > $ttmpabralog 2>&1'
echo -e "\e[0;36mPerforming indel realignment on normal BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'java -Xmx16G -jar /opt/software/abra.jar --in $gtmpdedupbam --out $gtmprealnbam \
#--ref $bwa_index --targets $regions --threads $cpu --working $gtmpabra > $gtmpabralog 2>&1'

# Sort and index re-aligned bam file (Note: Add clean up step)
echo -e "\e[0;36mSorting re-aligned tumor BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $toutbam $ttmprealnbam'
echo -e "\e[0;36mSorting re-aligned normal BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $goutbam $gtmprealnbam'
echo -e "\e[0;36mIndexing re-aligned tumor BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'samtools index $toutbam'
echo -e "\e[0;36mIndexing re-aligned normal BAM\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'samtools index $goutbam'

# Create pileup 
echo -e "\e[0;36mCreating tumor pileup\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'samtools mpileup -f $bwa_index -l $regions --output $tpileup $toutbam'
echo -e "\e[0;36mCreating normal pileup\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'samtools mpileup -f $bwa_index -l $regions --output $gpileup $goutbam'

# Call variants
echo -e "\e[0;36mCalling somatic variants\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'java -jar /opt/software/varscan/VarScan.v2.4.2.jar somatic $gpileup $tpileup $outsom --min-var-freq 0.01 --output-vcf 1'

# Process variant output
echo -e "\e[0;36mCreating SNP VCF\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'java -jar /opt/software/varscan/VarScan.v2.4.2.jar processSomatic $outvcfsnp'
echo -e "\e[0;36mCreating indel VCF\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'java -jar /opt/software/varscan/VarScan.v2.4.2.jar processSomatic $outvcfindel'
echo -e "\e[0;36mCombining SNP and indel VCF\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'vcfcombine $outvcfsnphc $outvcfindelhc > $outvcf'
echo -e "\e[0;36mCreating compressed VCF\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'bgzip -f $outvcf'
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'tabix -f -p vcf $outvcfgz'


# Run VarDict
echo -e "\e[0;36mCalling variants with VarDict\e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'VarDict -th $cpu -Q 10 -q 20 -G $bwa_index -f 0.01 -t -N $sample -b "$toutbam|$goutbam" -c 1 -S 2 -E 3 -g 4 $regions | /opt/software/VarDictJava/VarDict/testsomatic.R | /opt/software/VarDictJava/VarDict/var2vcf_somatic.pl -N "$tsample|$gsample" -f 0.01 > $outvardict'

# Run contanimation check script
echo -e "\e[0;36mRunning contamination check \e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'Rscript $cont -t $toutbam -g $goutbam -n $sample  -p $contpanel -o $contout'

# Run CNVkit
echo -e "\e[0;36mRunning CNVkit \e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'cnvkit.py batch -p 0 $toutbam --normal $goutbam --targets $regions --fasta $bwa_index --split --output-reference $refcnvkit --output-dir $cnv_dir --scatter'

# Run CNVkit filter
echo -e "\e[0;36mFiltering CNVkit results \e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c 'Rscript $cnvkit_filter $sample $cnvdir'

# Run Qualimap
echo -e "\e[0;36mRunning Qualimap on tumor bam \e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c ' qualimap bamqc --java-mem-size=10G -gd HUMAN -sd -gff $regions -bam $toutbam -outdir $qm_dir --outfile $tqmout'
echo -e "\e[0;36mRunning Qualimap on normal bam \e[0m"
#docker run -v /data:/data -v /tmp:/tmp -v /home:/home -it --env-file /home/anu/capp-seq-analysis/test.list anu9109/capp-seq bash -c ' qualimap bamqc --java-mem-size=10G -gd HUMAN -sd -gff $regions -bam $goutbam -outdir $qm_dir --outfile $gqmout'

# Run FastQC


echo -e "\e[0;36mDone! \e[0m"
exit 0

