#!/bin/sh
# this script is a shell script to run through the whole mutation calling pipeline
# input will be
# -s as sample name  ;
# -i as raw reads directory  --absolute directory here 
# -o as output_directory   --absolute directory here
# -d design file --a bed file indicating which regions are captured --majorly for QC purpose
# -n sample number(S1,S2,...depend on the name of sample file)
# -g gene list file for filtering   --absolute directory here
# -c SeqCNV control file  --absolute directory here
# -w run WES pipeline
# -f compressed format(miseq(default),nextseq,bz)
# -t frequency cutoff(default:0.005)
shellPath=$(cd "$(dirname "$0")"; pwd)
export PATH=${shellPath}/bin/jdk1.8.0_121/bin/:$PATH
fastpPath=${shellPath}/bin/fastp
bwaPath=${shellPath}/bin/bwa
samtoolsPath=${shellPath}/bin/samtools-1.2/samtools
gatkPath=${shellPath}/bin/gatk4/gatk
geneFile="ig"
compressedFormat="miseq"
cutoff=0.005
while getopts ":s:i:o:d:n:g:c:wf:t:" opt
do
    case $opt in
        s)
            sName=$OPTARG
        ;;
        i)
            iPath=$OPTARG
        ;;
        o)
            oPath=$OPTARG
        ;;
        d)
            designFile=$OPTARG
        ;;
        n)
            sNum=$OPTARG
        ;;
        g)
            geneFile=$OPTARG
        ;;
        c)
            CNVFile=$OPTARG
        ;;
        w)
            ifWES="Y"
        ;;
        f)
            compressedFormat=$OPTARG
        ;;
        t)
            cutoff=$OPTARG
        ;;
        ?)
            echo "unknown argument"
            exit 1
        ;;
    esac
done
if [ ! ${sName} ]; then
    echo "-s invalid"
    exit 1
fi
if [ ! ${iPath} ]; then
    echo "-i invalid"
    exit 1
fi
if [ ! ${oPath} ]; then
    echo "-o invalid"
    exit 1
fi
if [ ! ${designFile} ]; then
    echo "-d invalid"
    exit 1
fi


mkdir -p ${oPath}/${sName}
exec 1>>${oPath}/${sName}/${sName}_log.txt
exec 2>>${oPath}/${sName}/${sName}_log.txt

echo -e "${iPath}/${sName} startTime:::::::\c" ; date

echo "## Pipeline Starts to Run##"

touch ${oPath}/${sName}/${sName}_status.txt
date_start=$(date +%s)

echo "## Uncompress Fastq File ##" > ${oPath}/${sName}/${sName}_status.txt

if [ ${compressedFormat} == "bz" ]; then
    bunzip2 -k ${iPath}/${sName}.read1.bz2
    bunzip2 -k ${iPath}/${sName}.read2.bz2
    mv ${iPath}/${sName}.read1 ${iPath}/${sName}.read1.fq
    mv ${iPath}/${sName}.read2 ${iPath}/${sName}.read2.fq
elif [ ${compressedFormat} == "nextseq" ];then
    gunzip -c ${iPath}/${sName}_${sNum}_L001_R1_001.fastq.gz > ${iPath}/${sName}_${sNum}_L001_R1_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L002_R1_001.fastq.gz > ${iPath}/${sName}_${sNum}_L002_R1_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L003_R1_001.fastq.gz > ${iPath}/${sName}_${sNum}_L003_R1_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L004_R1_001.fastq.gz > ${iPath}/${sName}_${sNum}_L004_R1_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L001_R2_001.fastq.gz > ${iPath}/${sName}_${sNum}_L001_R2_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L002_R2_001.fastq.gz > ${iPath}/${sName}_${sNum}_L002_R2_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L003_R2_001.fastq.gz > ${iPath}/${sName}_${sNum}_L003_R2_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L004_R2_001.fastq.gz > ${iPath}/${sName}_${sNum}_L004_R2_001.fastq

    touch ${iPath}/${sName}.read1.fq
    cat ${iPath}/${sName}_${sNum}_L001_R1_001.fastq >> ${iPath}/${sName}.read1.fq
    rm ${iPath}/${sName}_${sNum}_L001_R1_001.fastq
    cat ${iPath}/${sName}_${sNum}_L002_R1_001.fastq >> ${iPath}/${sName}.read1.fq
    rm ${iPath}/${sName}_${sNum}_L002_R1_001.fastq
    cat ${iPath}/${sName}_${sNum}_L003_R1_001.fastq >> ${iPath}/${sName}.read1.fq
    rm ${iPath}/${sName}_${sNum}_L003_R1_001.fastq
    cat ${iPath}/${sName}_${sNum}_L004_R1_001.fastq >> ${iPath}/${sName}.read1.fq
    rm ${iPath}/${sName}_${sNum}_L004_R1_001.fastq

    touch ${iPath}/${sName}.read2.fq
    cat ${iPath}/${sName}_${sNum}_L001_R2_001.fastq >> ${iPath}/${sName}.read2.fq
    rm ${iPath}/${sName}_${sNum}_L001_R2_001.fastq
    cat ${iPath}/${sName}_${sNum}_L002_R2_001.fastq >> ${iPath}/${sName}.read2.fq
    rm ${iPath}/${sName}_${sNum}_L002_R2_001.fastq
    cat ${iPath}/${sName}_${sNum}_L003_R2_001.fastq >> ${iPath}/${sName}.read2.fq
    rm ${iPath}/${sName}_${sNum}_L003_R2_001.fastq
    cat ${iPath}/${sName}_${sNum}_L004_R2_001.fastq >> ${iPath}/${sName}.read2.fq
    rm ${iPath}/${sName}_${sNum}_L004_R2_001.fastq
else
    gunzip -c ${iPath}/${sName}_${sNum}_L001_R1_001.fastq.gz > ${iPath}/${sName}_${sNum}_L001_R1_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L001_R2_001.fastq.gz > ${iPath}/${sName}_${sNum}_L001_R2_001.fastq
    mv ${iPath}/${sName}_${sNum}_L001_R1_001.fastq ${iPath}/${sName}.read1.fq
    mv ${iPath}/${sName}_${sNum}_L001_R2_001.fastq ${iPath}/${sName}.read2.fq
fi

cd ${shellPath}
echo -e '${shellPath} is :::::::::\c';echo "${shellPath}"
echo "## Run fastp ##" > ${oPath}/${sName}/${sName}_status.txt
${fastpPath} -i ${iPath}/${sName}.read1.fq -I ${iPath}/${sName}.read2.fq -o ${oPath}/${sName}/${sName}_1_sequence.sheared.txt -O ${oPath}/${sName}/${sName}_2_sequence.sheared.txt -h ${oPath}/${sName}/${sName}.fastp.html -j ${oPath}/${sName}/${sName}.fastp.json

echo "## Run BWA ##" > ${oPath}/${sName}/${sName}_status.txt
${bwaPath} mem -R '@RG\tID:foo\tSM:'${sName}'\tLB:bar\tPL:illumina\tPU:run_std' -t 4 ${shellPath}/reference/hg19/hg19.fasta  ${oPath}/${sName}/${sName}_1_sequence.sheared.txt  ${oPath}/${sName}/${sName}_2_sequence.sheared.txt|${samtoolsPath} view -bS -o ${oPath}/${sName}/${sName}.bam -

mkdir -p ${oPath}/${sName}/BAM
echo "## Run samtools ##"
${samtoolsPath} sort -n  ${oPath}/${sName}/${sName}.bam  ${oPath}/${sName}/BAM/${sName}.sorted
${samtoolsPath} fixmate ${oPath}/${sName}/BAM/${sName}.sorted.bam  ${oPath}/${sName}/BAM/${sName}.matefixed.bam
${samtoolsPath} sort  ${oPath}/${sName}/BAM/${sName}.matefixed.bam  ${oPath}/${sName}/BAM/${sName}.matefixed.sorted
${gatkPath} --java-options "-Xmx4g" MarkDuplicates -I ${oPath}/${sName}/BAM/${sName}.matefixed.sorted.bam -M ${oPath}/${sName}/BAM/${sName}.metrics_file -O ${oPath}/${sName}/BAM/${sName}.rmdup.bam --VALIDATION_STRINGENCY LENIENT --VERBOSITY INFO --COMPRESSION_LEVEL 5 --CREATE_INDEX true --CREATE_MD5_FILE false --ASSUME_SORTED true --REMOVE_DUPLICATES true
${samtoolsPath} index ${oPath}/${sName}/BAM/${sName}.rmdup.bam

mkdir -p ${oPath}/${sName}/RECALIBRATION
echo "## Realignment ##" > ${oPath}/${sName}/${sName}_status.txt
echo "## Run BaseRecalibrator ##"
${gatkPath} --java-options "-Xmx4g" BaseRecalibrator -R ${shellPath}/reference/hg19/hg19.fasta --known-sites ${shellPath}/reference/dbsnp/dbsnp_132.hg19.vcf -I ${oPath}/${sName}/BAM/${sName}.rmdup.bam -O ${oPath}/${sName}/RECALIBRATION/${sName}.recal_data.cvs
echo "## Run ApplyBQSR ##"
${gatkPath} --java-options "-Xmx4g" ApplyBQSR -R ${shellPath}/reference/hg19/hg19.fasta -I ${oPath}/${sName}/BAM/${sName}.rmdup.bam -bqsr ${oPath}/${sName}/RECALIBRATION/${sName}.recal_data.cvs -O ${oPath}/${sName}/RECALIBRATION/${sName}.recal.bam
${samtoolsPath} index ${oPath}/${sName}/RECALIBRATION/${sName}.recal.bam

mkdir -p ${oPath}/${sName}/REALIGNMENT
#Indel realignment is no longer necessary for variant discovery if you plan to use a variant caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2.
${samtoolsPath} sort  ${oPath}/${sName}/RECALIBRATION/${sName}.recal.bam  ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted
${samtoolsPath} index ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam

#generate chrX reads ratio for gender QC
${samtoolsPath} view -o ${oPath}/${sName}/REALIGNMENT/${sName}.X.sam ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam chrX
${samtoolsPath} view -o ${oPath}/${sName}/REALIGNMENT/${sName}.sam ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam
wc -l ${oPath}/${sName}/REALIGNMENT/${sName}.X.sam > ${oPath}/${sName}/REALIGNMENT/${sName}.X.ratio
wc -l ${oPath}/${sName}/REALIGNMENT/${sName}.sam >> ${oPath}/${sName}/REALIGNMENT/${sName}.X.ratio
rm ${oPath}/${sName}/REALIGNMENT/${sName}.X.sam
rm ${oPath}/${sName}/REALIGNMENT/${sName}.sam

echo "## Run GATK HaplotypeCaller ##" > ${oPath}/${sName}/${sName}_status.txt
${gatkPath} --java-options "-Xmx8g" HaplotypeCaller -R ${shellPath}/reference/hg19/hg19.fasta -I ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.g.vcf -ERC GVCF -bamout ${oPath}/${sName}/REALIGNMENT/${sName}.HaplotypeCaller.bam
${gatkPath} --java-options "-Xmx4g" GenotypeGVCFs -R ${shellPath}/reference/hg19/hg19.fasta --variant ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.g.vcf -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.vcf
python ${shellPath}/bin/GATKScript/gatk_cis_trans.py -i ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.vcf -o ${oPath}/${sName}/${sName}.GATK.phase.txt
${gatkPath} SelectVariants -R ${shellPath}/reference/hg19/hg19.fasta -V ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.vcf -select-type SNP -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.vcf
${gatkPath} VariantFiltration -R ${shellPath}/reference/hg19/hg19.fasta -V ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "snp_filter" -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.mark.vcf
${gatkPath} SelectVariants -R ${shellPath}/reference/hg19/hg19.fasta -V ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.vcf -select-type INDEL -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.vcf
${gatkPath} VariantFiltration -R ${shellPath}/reference/hg19/hg19.fasta -V ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "indel_filter" -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.mark.vcf
${gatkPath} MergeVcfs -I ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.mark.vcf -I ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.mark.vcf -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.mark.vcf
rm ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.*
rm ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.*

mkdir -p ${oPath}/${sName}/SNP_CALL
python ${shellPath}/bin/GATKScript/format_GATK_SNP.py -i ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.mark.vcf -o ${oPath}/${sName}/SNP_CALL/${sName}.snp.vcf
python ${shellPath}/bin/annotate_filter/pipeline_filter_and_annotate.py -i ${oPath}/${sName}/SNP_CALL/${sName}.snp.vcf -o ${oPath}/${sName}/SNP_CALL/${sName}.snp.filt+anno.vcf -d ${sName} -t SNP -c ${cutoff}

mkdir -p ${oPath}/${sName}/INDEL_CALL
python ${shellPath}/bin/GATKScript/format_GATK_INDEL.py -i ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.mark.vcf -o ${oPath}/${sName}/INDEL_CALL/${sName}.indel.vcf
python ${shellPath}/bin/annotate_filter/pipeline_filter_and_annotate.py -i ${oPath}/${sName}/INDEL_CALL/${sName}.indel.vcf -o ${oPath}/${sName}/INDEL_CALL/${sName}.indel.filt+anno.vcf -d ${sName} -t INDEL -c ${cutoff}

if [ ${ifWES} ]; then
    python ${shellPath}/bin/annotate_filter/convert_to_readable_form.py -w -s ${oPath}/${sName}/SNP_CALL/${sName}.snp.filt+anno.vcf -i ${oPath}/${sName}/INDEL_CALL/${sName}.indel.filt+anno.vcf -o ${oPath}/${sName}/${sName}.format.analysis
    
    #sort the results by chr & pos
    sed '1d' ${oPath}/${sName}/${sName}.format.analysis > ${oPath}/${sName}/${sName}.format.header.removed.analysis
    head -n 1 ${oPath}/${sName}/${sName}.format.analysis > ${oPath}/${sName}/${sName}.format.sorted.analysis
    sort -k1,1 -k2n,2 ${oPath}/${sName}/${sName}.format.header.removed.analysis > ${oPath}/${sName}/${sName}.format.header.removed.sorted.analysis
    cat ${oPath}/${sName}/${sName}.format.header.removed.sorted.analysis >> ${oPath}/${sName}/${sName}.format.sorted.analysis
    rm ${oPath}/${sName}/${sName}.format.header.removed.analysis
    rm ${oPath}/${sName}/${sName}.format.header.removed.sorted.analysis
else
    #force call
    rm -rf ${oPath}/${sName}/FORCE
    mkdir ${oPath}/${sName}/FORCE

    #mpileup
    ${samtoolsPath} mpileup -Q 0 -f ${shellPath}/reference/hg19/hg19.fasta ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam > ${oPath}/${sName}/FORCE/${sName}.mpileup.txt

    #set the region you want to force call(here we set "chrX:38144608-38146560")
    python ${shellPath}/bin/force_call/force_call.py ${oPath}/${sName}/FORCE/${sName}.mpileup.txt chrX 38144608 38146560 > ${oPath}/${sName}/FORCE/${sName}.force_call

    #divide it into $1.force_call.snp.vcf and $1.force_call.indel.vcf
    python ${shellPath}/bin/force_call/format.py ${oPath}/${sName}/FORCE/${sName}.force_call

    #add annotations
    python ${shellPath}/bin/annotate_filter/pipeline_filter_and_annotate.py -i ${oPath}/${sName}/FORCE/${sName}.force_call.snp.vcf  -o ${oPath}/${sName}/FORCE/${sName}.force_call.snp.anno.vcf -d ${sName} -t SNP -f 1.1 -c 1 -s 1 -p
    python ${shellPath}/bin/annotate_filter/pipeline_filter_and_annotate.py -i ${oPath}/${sName}/FORCE/${sName}.force_call.indel.vcf  -o ${oPath}/${sName}/FORCE/${sName}.force_call.indel.anno.vcf -d ${sName} -t INDEL -f 1.1 -c 1 -s 1 -p

    #combine two parts as final results
    python ${shellPath}/bin/force_call/attach.py ${oPath}/${sName}/FORCE/${sName}.force_call.snp.anno.vcf ${oPath}/${sName}/SNP_CALL/${sName}.snp.filt+anno.vcf
    python ${shellPath}/bin/force_call/attach.py ${oPath}/${sName}/FORCE/${sName}.force_call.indel.anno.vcf ${oPath}/${sName}/INDEL_CALL/${sName}.indel.filt+anno.vcf
    python ${shellPath}/bin/annotate_filter/convert_to_readable_form.py -s ${oPath}/${sName}/SNP_CALL/${sName}.snp.filt+anno.vcf -i ${oPath}/${sName}/INDEL_CALL/${sName}.indel.filt+anno.vcf -o ${oPath}/${sName}/${sName}.format.sorted.analysis
fi

echo "## Run QC ##" > ${oPath}/${sName}/${sName}_status.txt
sh ${shellPath}/qc.sh  -s ${sName} -d ${designFile}  -i ${iPath} -o ${oPath}

echo "## Filter and Annotate ##" > ${oPath}/${sName}/${sName}_status.txt
#use gnomad to filter variants(max allele frequency > 0.005)
python ${shellPath}/bin/otherScript/gnomad_filter.py ${oPath}/${sName}/${sName}.format.sorted.analysis > ${oPath}/${sName}/${sName}.origin.analysis

#modify hgmd annotation
python ${shellPath}/bin/otherScript/hgmd_patch.py ${oPath}/${sName}/${sName}.origin.analysis > ${oPath}/${sName}/${sName}.hgmd.analysis

#modify Annovar wrong annotation and add DXDB annotation
python ${shellPath}/bin/otherScript/anno_modify.py ${oPath}/${sName}/${sName}.hgmd.analysis > ${oPath}/${sName}/${sName}.anno.modify.analysis

#add ClinVar annotation
python ${shellPath}/bin/otherScript/clinvar_anno.py ${oPath}/${sName}/${sName}.anno.modify.analysis > ${oPath}/${sName}/${sName}.anno.clinvar.analysis

if [ ${ifWES} ]; then

    python ${shellPath}/bin/otherScript/variant_filter_wes.py ${geneFile} ${shellPath}/reference/variant_filtered.txt ${oPath}/${sName}/${sName}.anno.clinvar.analysis > ${oPath}/${sName}/${sName}.filter.analysis 

    #internal filter(Baylor)
    python ${shellPath}/bin/otherScript/wes_internal_filter.py ${oPath}/${sName}/${sName}.filter.analysis > ${oPath}/${sName}/${sName}.inter.filter.analysis

    #internal annotation(Clinbytes)
    python ${shellPath}/bin/otherScript/wes_internal_annotate.py ${oPath}/${sName}/${sName}.inter.filter.analysis > ${oPath}/${sName}/${sName}.inter.filter.anno.analysis

    #rf percentage annotation
    python ${shellPath}/bin/otherScript/rf_anno.py ${oPath}/${sName}/${sName}.inter.filter.anno.analysis > ${oPath}/${sName}/${sName}.inter.filter.anno.rf.analysis
   
    #ExAC LOF annotation
    python ${shellPath}/bin/otherScript/ExAC_annotate.py ${oPath}/${sName}/${sName}.inter.filter.anno.rf.analysis > ${oPath}/${sName}/${sName}.inter.filter.anno.exac.analysis

    #GATK phase annotation
    python ${shellPath}/bin/GATKScript/GATK_phase_annotate.py ${oPath}/${sName}/${sName}.inter.filter.anno.exac.analysis  ${oPath}/${sName}/${sName}.GATK.phase.txt > ${oPath}/${sName}/${sName}.inter.filter.anno.exac.phase.analysis

    cp ${oPath}/${sName}/${sName}.inter.filter.anno.exac.phase.analysis ${oPath}/${sName}/${sName}.unclassified.analysis

    #classify analysis file into WES format
    python ${shellPath}/bin/otherScript/divide_WES.py ${oPath}/${sName}/${sName}.unclassified.analysis > ${oPath}/${sName}/${sName}.classified.analysis

    #change column order
    python ${shellPath}/bin/otherScript/column_change_wes.py ${oPath}/${sName}/${sName}.classified.analysis > ${oPath}/${sName}/${sName}.analysis

else
    
    python ${shellPath}/bin/otherScript/variant_filter.py ${geneFile} ${shellPath}/reference/variant_filtered.txt ${oPath}/${sName}/${sName}.anno.clinvar.analysis > ${oPath}/${sName}/${sName}.filter.analysis
    
    #change column order
    python ${shellPath}/bin/otherScript/column_change_cap.py ${oPath}/${sName}/${sName}.filter.analysis > ${oPath}/${sName}/${sName}.analysis

fi

python ${shellPath}/bin/otherScript/simplify.py ${oPath}/${sName}/${sName}.filter.analysis > ${oPath}/${sName}/${sName}.simplify.analysis

#SeqCNV
if [ ${CNVFile} ] ; then

    mkdir ${oPath}/${sName}/CNV/
    sh ${shellPath}/bin/SeqCNV/cnv.sh ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam ${CNVFile} ${oPath}/${sName}/CNV/ ${designFile}
    cp ${oPath}/${sName}/CNV/report/CNV_report.txt ${oPath}/${sName}/${sName}_CNV_report.txt

fi

#rm -rf ${oPath}/${sName}/CNV
rm -rf ${oPath}/${sName}/FORCE
rm ${oPath}/${sName}/*sheared.* 
rm ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.g.vcf*
rm ${oPath}/${sName}/${sName}.format*analysis
rm ${oPath}/${sName}/${sName}.hgmd.analysis
rm ${oPath}/${sName}/${sName}*filter*analysis
rm ${oPath}/${sName}/${sName}.anno.modify.analysis
rm ${oPath}/${sName}/${sName}.anno.clinvar.analysis
rm ${oPath}/${sName}/${sName}.bam
rm -rf ${oPath}/${sName}/BAM
rm -rf ${oPath}/${sName}/RECALIBRATION

rm ${iPath}/${sName}.read1.fq
rm ${iPath}/${sName}.read2.fq

echo -e "${iPath}/${sName} endTime:::::::\c" ; date

date_end=$(date +%s)

echo "## Completed(without batch CNV) ##" > ${oPath}/${sName}/${sName}_status.txt
echo "total time: $((date_end-date_start)) s"

