#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G

################   sub-function (Don't change any parameters)   ##########################################
function DG_extract() {
	## extract the DGs based on number
	awk '$20 == "DG:Z:18S,5.8S,'$DG'" {print}' *.cliques.t_o0.2.sam > 'DG'$DG'_18s_58s.sam'
	awk '$0~/^@/ || $1~/-701/' 'DG'$DG'_18s_58s.sam' > DG$DG'_c18s_58s.sam'
	awk '$0~/^@/ || $1~/-702/' 'DG'$DG'_18s_58s.sam' > DG$DG'_e18s_58s.sam'
   	cat header_only.sam DG$DG'_c18s_58s.sam' > DG$DG'_ctrl_18s_58s.sam' 
   	cat header_only.sam DG$DG'_e18s_58s.sam' > DG$DG'_exo_18s_58s.sam' 
  	
	python3 endpoints.py DG$DG'_ctrl_18s_58s.sam'
	python3 endpoints.py DG$DG'_exo_18s_58s.sam'       

	awk '$20 == "DG:Z:28SS,5.8S'$DG'" {print}' *.cliques.t_o0.2.sam > 'DG'$DG'_58s_28s.sam'
	awk '$0~/^@/ || $1~/-701/' 'DG'$DG'_58s_28s.sam' > DG$DG'_c58s_28s.sam'
	awk '$0~/^@/ || $1~/-702/' 'DG'$DG'_58s_28s.sam' > DG$DG'_e58s_28s.sam'
	cat header_only.sam DG$DG'_c58s_28s.sam' > DG$DG'_ctrl_28s_58s.sam' 
	cat header_only.sam DG$DG'_e58s_28s.sam' > DG$DG'_exo_28s_58s.sam' 
  	
	python3 endpoints.py DG$DG'_ctrl_28s_58s.sam'
	python3 endpoints.py DG$DG'_exo_28s_58s.sam'
	
	awk '$20 == "DG:Z:18S,28SS,'$DG'" {print}' *.cliques.t_o0.2.sam > 'DG'$DG'_18s_28s.sam'
	awk '$0~/^@/ || $1~/-701/' 'DG'$DG'_18s_28s.sam' > DG$DG'_c18s_28s.sam'
	awk '$0~/^@/ || $1~/-702/' 'DG'$DG'_18s_28s.sam' > DG$DG'_e18s_28s.sam'
	cat header_only.sam DG$DG'_c18s_28s.sam' > DG$DG'_ctrl_18s_28s.sam'
	cat header_only.sam DG$DG'_e18s_28s.sam' > DG$DG'_exo_18s_28s.sam'

	python3 endpoints.py DG$DG'_ctrl_18s_28s.sam'
	python3 endpoints.py DG$DG'_exo_18s_28s.sam'

	python3 endpoint_check.py $DG

	mkdir DG_$DG
	mv DG$DG* DG_$DG/
}

#################################
# rRNA ranges:			#
# 3654 - 5523, 18S		#
# 6600 - 6757, 5.8S		#
# 7924 - 12994, "28SS"		#
#################################
for n in {0..700};
do 
    DG=$n
    DG_extract
done

#################################

## generated sorted.bam for IGV ##

#    samtools view -bS -o 'DG'$DG'_ctrl.bam' 'DG'$DG'_ctrl.sam'
#    samtools sort -o 'DG'$DG'_ctrl_sorted.bam' 'DG'$DG'_ctrl.bam'
#    samtools index 'DG'$DG'_ctrl_sorted.bam'

#    samtools view -bS -o 'DG'$DG'_exo.bam' 'DG'$DG'_exo.sam'
#    samtools sort -o 'DG'$DG'_exo_sorted.bam' 'DG'$DG'_exo.bam'
#    samtools index 'DG'$DG'_exo_sorted.bam'
