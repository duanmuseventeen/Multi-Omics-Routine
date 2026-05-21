#!/bin/bash

INPUT_DIR="./2.fastq"
OUTPUT_DIR="./2.fastq_gsm"
# mkdir -p $OUTPUT_DIR

# 读取对应关系表（跳过注释行）
grep -v '^#' samples.txt | while read SAMPLE_ID SRR_LIST; do
    echo "正在处理样本: $SAMPLE_ID ..."
    
    R1_FILES=$(echo $SRR_LIST | sed "s|,|_1.fastq.gz $INPUT_DIR/|g" | sed "s|^|$INPUT_DIR/|")_1.fastq.gz
    cat $R1_FILES > ${OUTPUT_DIR}/${SAMPLE_ID}_1.fastq.gz
    
    R2_FILES=$(echo $SRR_LIST | sed "s|,|_2.fastq.gz $INPUT_DIR/|g" | sed "s|^|$INPUT_DIR/|")_2.fastq.gz
    cat $R2_FILES > ${OUTPUT_DIR}/${SAMPLE_ID}_2.fastq.gz
    
    echo "样本 $SAMPLE_ID 合并完成！"
done

echo "所有任务已结束。"
