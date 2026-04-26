#!/bin/bash
set -e

echo "Script is starting..."

MAX_JOBS=6
SRA_DIR="1.sra"
OUTPUT_DIR="2.fastq/"

FIFO_FILE="job_control.fifo"
rm -f   "$FIFO_FILE"
mkfifo  "$FIFO_FILE"
exec 6<>"$FIFO_FILE"
for ((i=0; i<MAX_JOBS; i++)); do 
    echo >&6
done

echo "Starting SRA file parallel dump with $MAX_JOBS cores."
START_TIME=$(date +%s)

for sra_file in "$SRA_DIR"/*.sra; do read -u6
    {
        echo "--> Processing: $(basename "$sra_file")"
        
        /home/fanyuanming/soft/anaconda3/envs/seq_py3.7/bin/fastq-dump --gzip --split-3 "$sra_file" --outdir "$OUTPUT_DIR"
        
        echo >&6
    } & done

wait

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "All files processed."
echo "Total runtime: $((DURATION / 60)) minutes and $((DURATION % 60)) seconds."

exec 6>&-
rm -f "$FIFO_FILE"
