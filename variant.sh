#!/bin/bash

docker pull anu9109/capp-seq

# Run VarDict
echo -e "\e[0;36mCalling variants with VarDict\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c '/opt/software/VarDictJava/build/install/VarDict/bin/VarDict -th $cpu -Q 10 -q 20 -G $bwa_index -f 0.01 -t -N $sample -b "$toutbam|$goutbam" -x 2000 -c 1 -S 2 -E 3 -g 4 $regions | /opt/software/VarDictJava/VarDict/testsomatic.R | /opt/software/VarDictJava/VarDict/var2vcf_paired.pl -N "$tsample|$gsample" -f 0.01 -P 0.9 -m 4.25 > $outvardict'

echo -e "\e[0;36mCreating compressed VCF\e[0m"
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'bgzip -f $outvardict'
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'tabix -f -p vcf $outvardictgz'

# Annotate VarDict VCF with Annovar
echo -e "\e[0;36mAnnotating variants with Annovar\e[0m"
docker run -v /home:/home -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'perl $anno/table_annovar.pl $outvardictgz $annodb -buildver hg19 -out $annovarout -protocol refGene,cosmic70,clinvar_20160302,icgc21,nci60,exac03,snp142,1000g2015aug_all,ljb26_all -operation g,f,f,f,f,f,f,f,f -nastring . -vcfinput --thread $cpu'

# Annotate VarDict VCF with Oncotator
echo -e "\e[0;36mAnnotating variants with Oncotator\e[0m"
docker run -v /var/run/docker.sock:/var/run/docker.sock  -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'docker run -i broadinstitute/oncotator:1.9.2.0 /root/oncotator_venv/bin/Oncotator -v -i VCF -o VCF --db-dir $oncodb $outvardictgz $oncoout hg19'
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'bgzip -f $oncoout'
docker run -v /data:/data -v /tmp:/tmp -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'tabix -f -p vcf $oncooutgz'

# Filter variants
echo -e "\e[0;36mFiltering variants\e[0m"
docker run -v /data:/data -v /tmp:/tmp  -i --env-file /home/anu/capp-seq-docker/env.list anu9109/capp-seq bash -c 'Rscript /home/variant_filter.R $annovcf $filtvcf /home'

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
