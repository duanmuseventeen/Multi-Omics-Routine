# 质控----------------------------------------------------------------------------------------------------------------------
# https://blog.csdn.net/sinat_32872729/article/details/94440265
# https://wenku.csdn.net/answer/00a814c941204123bd92917378236616
# https://www.jianshu.com/p/f83626fd1fa1

# 查看文件
tree -h
gunzip 1.fastq/*
  fastqc -t 24 1.fastq/* -o 2.fastqc/
  conda activate python3.11
multiqc 2.fastqc/
  
  # 根据质控结果，去接头---------------------------------------------------------------------------------------------------
for i in IR1 IR2 IR3 IR4 IR5 SHAM1 SHAM2 SHAM3 SHAM4 SHAM5 ;
do
fastp -g\
-i 1.fastq/${i}_1.fq.gz -I 1.fastq/${i}_2.fq.gz \
-o 3.fastp/${i}_1.fq.gz -O 3.fastp/${i}_2.fq.gz
done


# 序列比对------------------------------------------------------------------------------------------------------------------
# https://www.jianshu.com/p/681e02e7f9af
# https://daehwankimlab.github.io/hisat2/howto/
# https://www.zhihu.com/question/24259519/answer/2158984133

# build index---------------------------------------------------------------------------------------------------------------
If you use --snp, --ss, and/or --exon, hisat2-build will need about 200GB RAM for the human genome size as index building involves a graph construction. 
Otherwise, you will be able to build an index on your desktop with 8GB RAM.
# gencode human---------------------------------------------------------------------------------------------------------
conda activate python3.11
gunzip gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip GRCh38.p14.genome.fa.gz

hisat2_extract_splice_sites.py gencode.v47.chr_patch_hapl_scaff.annotation.gtf > genome.ss
hisat2_extract_exons.py gencode.v47.chr_patch_hapl_scaff.annotation.gtf > genome.exon

hisat2-build -p 24 --exon genome.exon --ss genome.ss GRCh38.p14.genome.fa GRCh38.p14.genome_tran
# hisat2-build -p 16 --large-index genome.fa genome.fa

# gencode mouse---------------------------------------------------------------------------------------------------------
gunzip gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip GRCm39.genome.fa.gz

hisat2_extract_splice_sites.py gencode.vM36.chr_patch_hapl_scaff.annotation.gtf > genome.ss
hisat2_extract_exons.py gencode.vM36.chr_patch_hapl_scaff.annotation.gtf > genome.exon

hisat2-build -p 24 --exon genome.exon --ss genome.ss GRCm39.genome.fa GRCm39.genome_tran

# esembl human-----------------------------------------------------------------------------------------------------------
# Build HGFM index with SNPs and transcripts
gunzip *
hisat2_extract_splice_sites.py Mus_musculus.GRCm39.113.gtf > genome.ss
hisat2_extract_exons.py Mus_musculus.GRCm39.113.gtf > genome.exon

hisat2-build -p 24 --snp genome.snp --haplotype genome.haplotype 
--exon genome.exon --ss genome.ss Mus_musculus.GRCm39.dna.toplevel.fa Mus_musculus.GRCm39.genome_snp_tran



for i in IR1 IR2 IR3 IR4 IR5 SHAM1 SHAM2 SHAM3 SHAM4 SHAM5 ;
do
nohup hisat2 -p 16 -t -x /home/FYM/genomeRef/hisat2/mouse/grcm38_snp_tran/genome_snp_tran -1 3.fastp/${i}_1.fq.gz -2 3.fastp/${i}_2.fq.gz -S 4.aligned/${i}.sam
done

nohup hisat2 -p 2 -t -x /home/FYM/genomeRef/hisat2/mouse/grcm38_snp_tran/genome_snp_tran -1 3.fastp/IR1_1.fq.gz -2 3.fastp/IR1_2.fq.gz -S 4.aligned/IR_1.sam &
nohup hisat2 -p 2 -t -x /home/FYM/genomeRef/hisat2/mouse/grcm38_snp_tran/genome_snp_tran -1 3.fastp/IR2_1.fq.gz -2 3.fastp/IR2_2.fq.gz -S 4.aligned/IR_2.sam &
nohup hisat2 -p 2 -t -x /home/FYM/genomeRef/hisat2/mouse/grcm38_snp_tran/genome_snp_tran -1 3.fastp/IR3_1.fq.gz -2 3.fastp/IR3_2.fq.gz -S 4.aligned/IR_3.sam &
nohup hisat2 -p 2 -t -x /home/FYM/genomeRef/hisat2/mouse/grcm38_snp_tran/genome_snp_tran -1 3.fastp/IR4_1.fq.gz -2 3.fastp/IR4_2.fq.gz -S 4.aligned/IR_4.sam &
nohup hisat2 -p 2 -t -x /home/FYM/genomeRef/hisat2/mouse/grcm38_snp_tran/genome_snp_tran -1 3.fastp/IR5_1.fq.gz -2 3.fastp/IR5_2.fq.gz -S 4.aligned/IR_5.sam &

nohup hisat2 -p 2 -t -x /home/FYM/genomeRef/hisat2/mouse/grcm38_snp_tran/genome_snp_tran -1 3.fastp/SHAM1_1.fq.gz -2 3.fastp/SHAM1_2.fq.gz -S 4.aligned/SHAM_1.sam &
nohup hisat2 -p 2 -t -x /home/FYM/genomeRef/hisat2/mouse/grcm38_snp_tran/genome_snp_tran -1 3.fastp/SHAM2_1.fq.gz -2 3.fastp/SHAM2_2.fq.gz -S 4.aligned/SHAM_2.sam &
nohup hisat2 -p 2 -t -x /home/FYM/genomeRef/hisat2/mouse/grcm38_snp_tran/genome_snp_tran -1 3.fastp/SHAM3_1.fq.gz -2 3.fastp/SHAM3_2.fq.gz -S 4.aligned/SHAM_3.sam &
nohup hisat2 -p 2 -t -x /home/FYM/genomeRef/hisat2/mouse/grcm38_snp_tran/genome_snp_tran -1 3.fastp/SHAM4_1.fq.gz -2 3.fastp/SHAM4_2.fq.gz -S 4.aligned/SHAM_4.sam &
nohup hisat2 -p 2 -t -x /home/FYM/genomeRef/hisat2/mouse/grcm38_snp_tran/genome_snp_tran -1 3.fastp/SHAM5_1.fq.gz -2 3.fastp/SHAM5_2.fq.gz -S 4.aligned/SHAM_5.sam &
  
# samtools
for i in IR_1 IR_2 IR_3 IR_4 IR_5 SHAM_1 SHAM_2 SHAM_3 SHAM_4 SHAM_5 ;
do
samtools view -S --threads 12 4.aligned/${i}.sam -b > 4.aligned/${i}.bam
samtools sort --threads 12 4.aligned/${i}.bam -o 4.aligned/${i}_sorted.bam # 如果后续用htseq处理，要加-n参数 https://www.jianshu.com/p/108b198adfc6
samtools index 4.aligned/${i}_sorted.bam
done

# reads计数-------------------------------------------------------------------------------------------------------------------
# https://www.jianshu.com/p/e34ab865f055
# https://www.plob.org/article/11504.html
# http://www.biocloudservice.com/wordpress/?p=58919
# 在基因水平上，常用的软件为HTSeq-count，featureCounts，BEDTools, Qualimap, Rsubread, GenomicRanges等。
# HTSeq 无特殊情况，选择union比对模式
1. HTSeq是对有参考基因组的转录组测序数据进行表达量分析的，其输入文件必须有SAM和GTF文件。
2. 一般情况下HTSeq得到的Counts结果会用于下一步不同样品间的基因表达量差异分析，而不是一个样品内部基因的表达量比较。因此，HTSeq设置了-a参数的默认值10，来忽略掉比对到多个位置的reads信息，其结果有利于后续的差异分析。
3. 输入的GTF文件中不能包含可变剪接信息，否则HTSeq会认为每个可变剪接都是单独的基因，导致能比对到多个可变剪接转录本上的reads的计算结果是ambiguous，从而不能计算到基因的count中。即使设置-i参数的值为transcript_id，其结果一样是不准确的，只是得到transcripts的表达量。

# htseq seperately==============================================================
# https://www.jianshu.com/p/682916b324a6
# http://www.chenlianfu.com/?p=2438
for i in IR_1 IR_2 IR_3 IR_4 IR_5 SHAM_1 SHAM_2 SHAM_3 SHAM_4 SHAM_5 ;
do
htseq-count -f bam -r pos  -i gene_name -s no -m union 5.bam/${i}_sorted.bam /home/FYM/genomeRef/hisat2/mouse/Mus_musculus.GRCm38.84.gtf > 6.counts/${i}_count.txt
done
# htseq multi-thread============================================================
htseq-count -n 10 -f bam -r pos  -i gene_name -s no -m union \
5.bam/IR_1_sorted.bam \
5.bam/IR_2_sorted.bam \
5.bam/IR_3_sorted.bam \
5.bam/IR_4_sorted.bam \
5.bam/IR_5_sorted.bam \
5.bam/SHAM_1_sorted.bam \
5.bam/SHAM_2_sorted.bam \
5.bam/SHAM_3_sorted.bam \
5.bam/SHAM_4_sorted.bam \
5.bam/SHAM_5_sorted.bam \
/home/FYM/genomeRef/hisat2/mouse/Mus_musculus.GRCm38.84.gtf \
--with-header > 7.counts/counts.tsv
# feature-counts================================================================
# https://subread.sourceforge.net/featureCounts.html
-a <string>         Name of an annotation file. GTF/GFF format by default. See
-F option for more format information. Inbuilt annotations
(SAF format) is available in 'annotation' directory of the
package. Gzipped file is also accepted.

-o <string>         Name of output file including read counts. A separate file
including summary statistics of counting results is also
included in the output ('<string>.summary'). Both files
are in tab delimited format.

input_file1 [input_file2] ...   A list of SAM or BAM format files. They can be
either name or location sorted. If no files provided,
<stdin> input is expected. Location-sorted paired-end reads
are automatically sorted by read names.

-p                  If specified, libraries are assumed to contain paired-end
reads. For any library that contains paired-end reads, the
'countReadPairs' parameter controls if read pairs or reads
should be counted.

-B                  Only count read pairs that have both ends aligned.

-C                  Do not count read pairs that have their two ends mapping
to different chromosomes or mapping to same chromosome
but on different strands.

# Long reads
-L                  Count long reads such as Nanopore and PacBio reads. Long
read counting can only run in one thread and only reads
(not read-pairs) can be counted. There is no limitation on
the number of 'M' operations allowed in a CIGAR string in
long read counting.

-t <string>         Specify feature type(s) in a GTF annotation. If multiple
types are provided, they should be separated by ',' with
no space in between. 'exon' by default. Rows in the
annotation with a matched feature will be extracted and
used for read mapping.

# -t gene or -t exon (default) are both ok
featureCounts -a /home/FYM/genomeRef/hisat2/mouse/Mus_musculus.GRCm38.84.gtf -o 8.featurecounts/counts.txt -g 'gene_name' -T 10 -p -B -C 5.bam/IR_1_sorted.bam 5.bam/IR_2_sorted.bam 5.bam/IR_3_sorted.bam 5.bam/IR_4_sorted.bam 5.bam/IR_5_sorted.bam 5.bam/SHAM_1_sorted.bam 5.bam/SHAM_2_sorted.bam 5.bam/SHAM_3_sorted.bam 5.bam/SHAM_4_sorted.bam 5.bam/SHAM_5_sorted.bam

multiqc result.txt.summary
# to understand the param -t -f
featureCounts -a /home/FYM/genomeRef/hisat2/mouse/Mus_musculus.GRCm38.84.gtf -o 8.featurecounts/counts_f.txt -g 'gene_name' -T 10 -p -B -C -f 5.bam/IR_1_sorted.bam
featureCounts -a /home/FYM/genomeRef/hisat2/mouse/Mus_musculus.GRCm38.84.gtf -o 8.featurecounts/counts_texon.txt -g 'gene_name' -T 10 -p -B -C -t exon 5.bam/IR_1_sorted.bam 
featureCounts -a /home/FYM/genomeRef/hisat2/mouse/Mus_musculus.GRCm38.84.gtf -o 8.featurecounts/counts_ttranscript.txt -g 'gene_name' -T 10 -p -B -C -t transcript 5.bam/IR_1_sorted.bam 
featureCounts -a /home/FYM/genomeRef/hisat2/mouse/Mus_musculus.GRCm38.84.gtf -o 8.featurecounts/counts.txt -g 'gene_name' -T 10 -p -B -C -f -t exon 5.bam/IR_1_sorted.bam
featureCounts -a /home/FYM/genomeRef/hisat2/mouse/Mus_musculus.GRCm38.84.gtf -o 8.featurecounts/counts_tgene.txt -g 'gene_name' -T 10 -p -B -C -t gene 5.bam/IR_1_sorted.bam 















