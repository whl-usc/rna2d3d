#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=32G

# This part for prigap1/prigapm/prihomo reads.
extract_random_reads() {
	local input_bam="$1"
	local num_reads="$2"
	local read_type="$3"
	local output_bam="$4"
	{
		samtools view -@ 8 -H "$input_bam"
		case "$read_type" in
			all)
				samtools view -@ 8 -F 4 "$input_bam" | shuf | head -n "$num_reads"
				;;
			rRNA)
				samtools view -@ 8 -F 4 "$input_bam" | awk '$3 == "hs45S"' | shuf | head -n "$num_reads"
				;;
			no_rRNA)
				samtools view -@ 8 -F 4 "$input_bam" | awk '$3 != "hs45S"' | shuf | head -n "$num_reads"
				;;
			*)
				echo "Unknown read_type: $read_type" >&2
				exit 1
				;;
			esac
	} | samtools view -@ 8 -bS -o "$output_bam" -
}

# === Main ===
input_bam="HS-HEK293T_AMT-0.5_T4-24h_exo-0h_prigap1.bam"


max_reads=$(samtools view -@ 8 -c "$input_bam")
thresholds=(300000000 200000000 100000000 50000000)
output_files=("prigap1_300M_all.bam" "prigap1_200M_all.bam" "prigap1_100M_all.bam" "prigap1_50M_all.bam")

for i in "${!thresholds[@]}"; do
    if [ "$max_reads" -ge "${thresholds[$i]}" ]; then
        extract_random_reads "$input_bam" "${thresholds[$i]}" all "${output_files[$i]}"
    fi
done

extract_random_reads $input_bam 10000000 all prigap1_10M_all.bam
extract_random_reads $input_bam 1000000 all prigap1_1M_all.bam
extract_random_reads $input_bam 100000 all prigap1_100K_all.bam
extract_random_reads $input_bam 10000 all prigap1_10K_all.bam
extract_random_reads $input_bam 1000 all prigap1_1K_all.bam

max_reads=$(samtools view -@ 8 -F 4 "$input_bam" | awk '$3 == "hs45S"' | wc -l)
thresholds=(300000000 200000000 100000000 50000000)
output_files=("prigap1_300M_rRNA.bam" "prigap1_200M_rRNA.bam" "prigap1_100M_rRNA.bam" "prigap1_50M_rRNA.bam")

for i in "${!thresholds[@]}"; do
    if [ "$max_reads" -ge "${thresholds[$i]}" ]; then
        extract_random_reads "$input_bam" "${thresholds[$i]}" all "${output_files[$i]}"
    fi
done

extract_random_reads $input_bam 10000000 rRNA prigap1_10M_rRNA.bam
extract_random_reads $input_bam 1000000 rRNA prigap1_1M_rRNA.bam
extract_random_reads $input_bam 100000 rRNA prigap1_100K_rRNA.bam
extract_random_reads $input_bam 10000 rRNA prigap1_10K_rRNA.bam
extract_random_reads $input_bam 1000 rRNA prigap1_1K_rRNA.bam

max_reads=$(samtools view -@ 8 -F 4 "$input_bam" | awk '$3 != "hs45S"' | wc -l)

max_reads=$(samtools view -@ 8 -c "$input_bam")
thresholds=(100000000 10000000)
output_files=("prigap1_100M_rRNA.bam" "prigap1_10M_rRNA.bam")

for i in "${!thresholds[@]}"; do
    if [ "$max_reads" -ge "${thresholds[$i]}" ]; then
        extract_random_reads "$input_bam" "${thresholds[$i]}" all "${output_files[$i]}"
    fi
done

extract_random_reads $input_bam 1000000 no_rRNA prigap1_1M_norRNA.bam
extract_random_reads $input_bam 100000 no_rRNA prigap1_100K_norRNA.bam
extract_random_reads $input_bam 10000 no_rRNA prigap1_10K_norRNA.bam
extract_random_reads $input_bam 1000 no_rRNA prigap1_1K_norRNA.bam

# Check the statistics...
print_bam_stats() {
    local file="$1"
    
    echo "$file:" 
    echo "Size bam (bytes): "
    samtools view "$file" | wc -c

    stats=$(samtools stats "$file")

    echo -n "Average length: "
    echo "$stats" | awk -F '\t' '/^SN\taverage length:/ { print $3 }'
    
    echo -n "Mapped reads: "
    echo "$stats" | awk -F '\t' '/^SN\treads mapped:/ { print $3 }'

    echo ""
}

# Process prigap1_*_all.bam files
for file in prigap1_*_all.bam; do 
    print_bam_stats "$file"
done

# Process prigap1_*_rRNA.bam files
for file in prigap1_*_rRNA.bam; do 
    print_bam_stats "$file"
done

# Process prigap1_*_norRNA.bam files
for file in prigap1_*_norRNA.bam; do 
    print_bam_stats "$file"
done

# Compress all the bam files into bedz format
for file in prigap1_*.bam; do
	/project/zhipengl_72/wilsonhl/whl/pipeline/bin/rna2d3d/scripts/bam2bedz -r gap1 $file
done 

for file in prigap1_*_all.bed.gz; do
	echo "$file:"
	echo "Size bed.gz:" && less $file | wc -l
done

for file in prigap1_*_rRNA.bed.gz; do
	echo "$file:" 
	echo "Size bed.gz:" && less $file | wc -l
done

for file in prigap1_*_norRNA.bed.gz; do
	echo "$file:"
	echo "Size bed.gz:" && less $file | wc -l
done

# This only works for pritrans reads.
extract_random_read_pairs() {
	local input_bam="$1"
	local num_pairs="$2"
	local output_bam="$3"

	tmp_bam_sorted=$(mktemp --suffix=.bam)
	samtools sort -n -@ 8 -o "$tmp_bam_sorted" "$input_bam"
	samtools view -@ 8 -H "$input_bam" > header.sam
	samtools view -@ 8 "$tmp_bam_sorted" \
	| paste - - \
	| shuf \
	| head -n "$num_pairs" \
	| tr '\t' '\n' >> body.sam
	cat header.sam body.sam | samtools view -@ 8 -bS -o "$output_bam" -
	rm header.sam body.sam "$tmp_bam_sorted"
}


for line in HS-HEK293T_AMT-0.5_T4-24h_exo-0h_prigap1_filtered.bam; { samtools view -H input.bam; samtools view -F 4 input.bam; head -n 1 } | samtools view -bS -o 1.bam -

