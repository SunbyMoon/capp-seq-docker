#!/bin/bash

docker pull anu9109/capp-seq

# Run bbduk to trim adapters
echo -e "\e[0;36mRunning bbduk trimming on tumor bam \e[0m"
docker run -v /data:/data -v /tmp:/tmp  -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c './opt/software/bbmap/bbduk.sh in=$tinfastq1 in2=$tinfastq2 out=$tfastq1 out2=$tfastq2 \ 
ref=/opt/software/bbmap/resources/truseq.fa.gz stats=$tstats bhist=$tbhist qhist=$tqhist aqhist=$taqhist lhist=$tlhist gchist=$tgchist gcbins=auto threads=$cpu minlen=25 qtrim=rl trimq=10 ktrim=r k=25 \ 
mink=11 hdist=1 overwrite=true tbo=t tpe=t'
echo -e "\e[0;36mRunning bbduk trimming on normal bam \e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c './opt/software/bbmap/bbduk.sh in=$ginfastq1 in2=$ginfastq2 out=$gfastq1 out2=$gfastq2 \
ref=/opt/software/bbmap/resources/truseq.fa.gz stats=$gstats bhist=$gbhist qhist=$gqhist aqhist=$gaqhist lhist=$glhist gchist=$ggchist gcbins=auto threads=$cpu minlen=25 qtrim=rl trimq=10 ktrim=r k=25 \
mink=11 hdist=1 overwrite=true tbo=t tpe=t'

# Run FastQC
echo -e "\e[0;36mRunning FastQC \e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'fastqc -o $fastqc_dir $tfastq1 $tfastq2 $gfastq1 $gfastq2 -t $cpu'

# Alignment
echo -e "\e[0;36mCreating tumor BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -v /dev:/dev-i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'bwa mem -M -t $cpu -R $tRGR $bwa_index  $tfastq1 $tfastq2 | \
samblaster --splitterFile >(samtools view -S -u /dev/stdin | sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $ttmpsortbam /dev/stdin) \
--discordantFile >(samtools view -S -u /dev/stdin | sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $ttmpsortbam /dev/stdin) | \
samtools view -S -u /dev/stdin | \
sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $ttmpsortbam /dev/stdin'
echo -e "\e[0;36mCreating normal BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -v /dev:/dev-i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'bwa mem -M -t $cpu -R $gRGR $bwa_index  $gfastq1 $gfastq2 | \
samblaster --splitterFile >(samtools view -S -u /dev/stdin | sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $gtmpsortbam /dev/stdin) \
--discordantFile >(samtools view -S -u /dev/stdin | sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $gtmpsortbam /dev/stdin) | \
samtools view -S -u /dev/stdin | \
sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $gtmpsortbam /dev/stdin'

# Index bam file
echo -e "\e[0;36mIndexing tumor BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'samtools index $ttmpsortbam'
echo -e "\e[0;36mIndexing normal BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'samtools index $gtmpsortbam'

# Mark duplicates
echo -e "\e[0;36mMarking duplicates in tumor BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'sambamba markdup -t $cpu --overflow-list-size 1000000 --hash-table-size 1000000 \
$ttmpsortbam $ttmpdedupbam'
echo -e "\e[0;36mMarking duplicates in normal BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'sambamba markdup -t $cpu --overflow-list-size 1000000 --hash-table-size 1000000 \
$gtmpsortbam $gtmpdedupbam'

# Index de-duplicated bam file (Note: Add clean up step)
echo -e "\e[0;36mIndexing de-duplicated tumor BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'samtools index $ttmpdedupbam'
echo -e "\e[0;36mIndexing de-duplicated normal BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'samtools index $gtmpdedupbam'

# Indel realignment with Abra (Note: Add clean up step)
echo -e "\e[0;36mPerforming indel realignment on tumor BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'java -Xmx16G -jar /opt/software/abra.jar --in $ttmpdedupbam --out $ttmprealnbam \
--ref $bwa_index --targets $regions --threads $cpu --working $ttmpabra > $ttmpabralog 2>&1'
echo -e "\e[0;36mPerforming indel realignment on normal BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'java -Xmx16G -jar /opt/software/abra.jar --in $gtmpdedupbam --out $gtmprealnbam \
--ref $bwa_index --targets $regions --threads $cpu --working $gtmpabra > $gtmpabralog 2>&1'

# Sort and index re-aligned bam file (Note: Add clean up step)
echo -e "\e[0;36mSorting re-aligned tumor BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $toutbam $ttmprealnbam'
echo -e "\e[0;36mSorting re-aligned normal BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'sambamba sort -t $cpu -m 10G --tmpdir $temp_dir -o $goutbam $gtmprealnbam'
echo -e "\e[0;36mIndexing re-aligned tumor BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'samtools index $toutbam'
echo -e "\e[0;36mIndexing re-aligned normal BAM\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'samtools index $goutbam'

# Run Qualimap
echo -e "\e[0;36mRunning Qualimap on tumor bam \e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'qualimap bamqc --java-mem-size=10G -gd HUMAN -sd -gff $regions -bam $toutbam -outdir $tqm_dir --outformat HTML'
echo -e "\e[0;36mRunning Qualimap on normal bam \e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'qualimap bamqc --java-mem-size=10G -gd HUMAN -sd -gff $regions -bam $goutbam -outdir $gqm_dir --outformat HTML'

# Run VarDict
echo -e "\e[0;36mCalling variants with VarDict\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'VarDict -th $cpu -Q 10 -q 20 -G $bwa_index -f 0.01 -t -N $sample -b "$toutbam|$goutbam" -x 2000 -c 1 -S 2 -E 3 -g 4 $regions | /opt/software/VarDictJava/VarDict/testsomatic.R | /opt/software/VarDictJava/VarDict/var2vcf_somatic.pl -N "$tsample|$gsample" -f 0.01 -P 0.9 -m 4.25 > $outvardict'

echo -e "\e[0;36mCreating compressed VCF\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'bgzip -f $outvardict'
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'tabix -f -p vcf $outvardictgz'

# Annotate VarDict VCF with Annovar
echo -e "\e[0;36mAnnotating variants with Annovar\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'perl $anno/table_annovar.pl $outvardictgz $annodb -buildver hg19 -out $annovarout -protocol refGene,cosmic70,clinvar_20160302,icgc21,nci60,exac03,snp142,1000g2015aug_all,ljb26_all -operation g,f,f,f,f,f,f,f,f -nastring . -vcfinput --thread $cpu'

# Annotate VarDict VCF with Oncotator
echo -e "\e[0;36mAnnotating variants with Oncotator\e[0m"
docker run -v /var/run/docker.sock:/var/run/docker.sock  -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'docker run -i broadinstitute/oncotator:1.9.2.0 /root/oncotator_venv/bin/Oncotator -v -i VCF -o VCF --db-dir $oncodb $outvardictgz $oncoout hg19'
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'bgzip -f $oncoout'
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'tabix -f -p vcf $oncooutgz'

# Filter variants
echo -e "\e[0;36mFiltering variants\e[0m"
docker run -v /data:/data -v /tmp:/tmp  -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'Rscript /home/variant_filter.R $annovcf $filtvcf'

# Run contamination check script
echo -e "\e[0;36mRunning contamination check \e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'Rscript /home/cont.R -t $toutbam -g $goutbam -n $sample  -p /home/contPanel.csv -o $contout'

# Run CNVkit
echo -e "\e[0;36mRunning CNVkit \e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'cnvkit.py batch -p 0 $toutbam --normal $goutbam --targets $cnvregions --fasta $bwa_index --split --output-reference $refcnvkit --output-dir $cnv_dir --scatter'

# Run CNVkit filter
echo -e "\e[0;36mFiltering CNVkit results \e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'Rscript /home/cnvkit_filter.R $sample $cnvkitcns $cnvkitout'

# Generate patient report
echo -e "\e[0;36mCreating patient report \e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list  anu9109/capp-seq bash -c 'Rscript --vanilla /home/arep.R -v $filtvcf -r $annovcf -C $cnvkitout -i $sample -o $report_dir -d /home'

# Generate QC report
echo -e "\e[0;36mCreating QC report \e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'Rscript --vanilla /home/createQCreport.R --fastqc $fastqc_dir --qualimap $qm_dir --id $sample --id_tumor $id_tumor --id_blood $id_blood --sample_config $qc_config --baseSpace_token 2434e39ec0774994a5ac22544cdca40f --sample_tracking $contout --outdir $report_dir --scriptdir /home'

# Cleaning temp folder
echo -e "\e[0;36mCleaning /tmp space used by Abra \e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'rm -r $ttmpabra; rm -r $gtmpabra;rm /tmp/$tsample*; rm /tmp/$gsample*'

# Docker cleanup
echo -e "\e[0;36mRemoving unused Docker containers & images \e[0m"
docker rm $(docker ps -aq -f status=exited)
#docker rmi $(docker images -f "dangling=true" -q)

echo -e "\e[0;36mDone! \e[0m"
exit 0


