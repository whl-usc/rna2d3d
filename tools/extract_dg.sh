#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G

################   sub-function (Don't change any parameters)   ##########################################
function DG_extract() {
	
    ## extract the DGs based on number
	awk '$20 == "DG:Z:18S,18S,'$DG'" {print}' *.cliques.t_o0.2.sam > 'DG'$DG'_18s.sam'
	awk '$0~/^@/ || $1~/-control/' 'DG'$DG'_18s.sam' > DG$DG'_c18s.sam'
	awk '$0~/^@/ || $1~/-exo/' 'DG'$DG'_18s.sam' > DG$DG'_e18s.sam'
	cat header_only.sam DG$DG'_c18s.sam' > DG$DG'_ctrl_18s.sam' 
   	cat header_only.sam DG$DG'_e18s.sam' > DG$DG'_exo_18s.sam' 
  	
	python3 endpoints.py DG$DG'_ctrl_18s.sam'
	python3 endpoints.py DG$DG'_exo_18s.sam'       

	awk '$20 == "DG:Z:5.8S,5.8S,'$DG'" {print}' *.cliques.t_o0.2.sam > 'DG'$DG'_5.8s.sam'
	awk '$0~/^@/ || $1~/-control/' 'DG'$DG'_5.8s.sam' > DG$DG'_c5.8s.sam'
	awk '$0~/^@/ || $1~/-exo/' 'DG'$DG'_5.8s.sam' > DG$DG'_e5.8s.sam'
	cat header_only.sam DG$DG'_c5.8s.sam' > DG$DG'_ctrl_5.8s.sam' 
	cat header_only.sam DG$DG'_e5.8s.sam' > DG$DG'_exo_5.8s.sam' 
  	
	python3 endpoints.py DG$DG'_ctrl_5.8s.sam'
	python3 endpoints.py DG$DG'_exo_5.8s.sam'

	awk '$20 == "DG:Z:28SS,28SS,'$DG'" {print}' *.cliques.t_o0.2.sam > 'DG'$DG'_28s.sam'
	awk '$0~/^@/ || $1~/-control/' 'DG'$DG'_28s.sam' > DG$DG'_c28s.sam'
	awk '$0~/^@/ || $1~/-exo/' 'DG'$DG'_28s.sam' > DG$DG'_e28s.sam'
	cat header_only.sam DG$DG'_c28s.sam' > DG$DG'_ctrl_28s.sam' 
	cat header_only.sam DG$DG'_e28s.sam' > DG$DG'_exo_28s.sam' 
  	
	python3 endpoints.py DG$DG'_ctrl_28s.sam'
	python3 endpoints.py DG$DG'_exo_28s.sam'
	
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
for n in {0..100};
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
