#!/bin/bash
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH -N 1

# this code is to analysis the genetic variants effect on DNA methylation. calculate the average methylation and depth in genetic variants regions
# prepare paramter file before run this code: each row contains "file" and "binID", sperate by "\t"
# run this code like: sbatch --array=1-n%m this_code.sh this_code.paramter


read file binID < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

eval $( spack load --sh bedtools2)

pair=$(basename ${file} |awk -F '_' '{print $1"_"$2}')
binSize="100"
binNum="101"

sam1=$(echo "${pair}" |awk -F '_' '{print $1}')
sam2=$(echo "${pair}" |awk -F '_' '{print $2}')

# m1C1="/scratch/twlab/shihua/pro01-PanGenome/2.1_iPSC_WGBS/call_methlation_by_methyGraph/06_iPSC_methyGraphLiftoverToHg38_DeduplicatedFiles/${sam1}-Clo1_chr1-22xy_methyGraLifToHg38_Deduplicated.methylC.bed"
# m1C2=$(echo "${m1C1}" | sed "s/${sam1}-Clo1/${sam1}-Clo2/g")
# m2C1=$(echo "${m1C1}" | sed "s/${sam1}-Clo1/${sam2}-Clo1/g")
# m2C2=$(echo "${m1C1}" | sed "s/${sam1}-Clo1/${sam2}-Clo2/g")

m1C1="/scratch/twlab/shihua/pro01-PanGenome/2.1_iPSC_WGBS/call_methlation_by_methyGraph/06_iPSC_methyGraphLiftoverToHg38_DeduplicatedFiles/PAN001-Clo1_chr1-22xy_methyGraLifToHg38_Deduplicated.methylC.bed"
m1C2="/scratch/twlab/shihua/pro01-PanGenome/2.1_iPSC_WGBS/call_methlation_by_methyGraph/06_iPSC_methyGraphLiftoverToHg38_DeduplicatedFiles/PAN001-Clo2_chr1-22xy_methyGraLifToHg38_Deduplicated.methylC.bed"
m2C1="/scratch/twlab/shihua/pro01-PanGenome/2.1_iPSC_WGBS/call_methlation_by_methyGraph/06_iPSC_methyGraphLiftoverToHg38_DeduplicatedFiles/PAN027-Clo1_chr1-22xy_methyGraLifToHg38_Deduplicated.methylC.bed"
m2C2="/scratch/twlab/shihua/pro01-PanGenome/2.1_iPSC_WGBS/call_methlation_by_methyGraph/06_iPSC_methyGraphLiftoverToHg38_DeduplicatedFiles/PAN027-Clo2_chr1-22xy_methyGraLifToHg38_Deduplicated.methylC.bed"

pd=$(realpath ./)
inter_dir="${pd}/inter-files"
results_dir="${pd}/results"
mkdir ${inter_dir}
mkdir ${results_dir}

name=$(basename $file |sed "s/.bed/_${binNum}bins_bin${binSize}bp_binID-${binID}/g")

echo -e "${file}\n${name}"

echo "#=========================================================start step 1: generate 100 bins ============================================="
out1="${inter_dir}/step1_${name}.bed"
shuf -n 100000 ${file} |awk -v num="${binNum}" -v binSize="${binSize}" -v binID="${binID}" '{start=($2+($3-$2)/2)-binSize/2+binSize*(binID - (num+1)/2);end=($2+($3-$2)/2)+binSize/2+binSize*(binID- (num+1)/2);binPos=binSize*(binID - (num+1)/2);printf("%s\t%.0f\t%.0f\t%s\t%s\n",$1,start,end,$4,binPos)}' | sortBed -i > ${out1} 
echo -e "step 1 done, output: ${out1}"
date

echo "#========================================================= step 2: call methylation ============================================="
out2_1="${inter_dir}/step2_${name}_${sam1}_Clo1_Cal-meth.txt"
out2_2="${inter_dir}/step2_${name}_${sam1}_Clo2_Cal-meth.txt"
out2_3="${inter_dir}/step2_${name}_${sam2}_Clo1_Cal-meth.txt"
out2_4="${inter_dir}/step2_${name}_${sam2}_Clo2_Cal-meth.txt"
out2_5="${inter_dir}/step2_${name}_meth-dept-diff"
out2="${results_dir}/step2_${name}_meth-diff_summary.txt"

process_file() {
	input_file=$1
	output_file=$2
	awk '$NF>=5 && $NF <= 150' "${input_file}" | bedtools map -a "${out1}" -b - -c 5,7 -o mean |awk '$NF!="."' |awk '{printf("%s\t%s\t%s\t%s\t%s\t%.4f\t%.4f\n", $1,$2,$3,$4,$5,$(NF-1),$NF)}' > "${output_file}"
}

process_file "${m1C1}" "${out2_1}"
process_file "${m1C2}" "${out2_2}"
process_file "${m2C1}" "${out2_3}"
process_file "${m2C2}" "${out2_4}"


#################################################################### Merge by R, R code ###################################################################
Rscript -e ' 
# R code start here
## merge files by position

rm(list = ls())
set.seed(0)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to extract parameter value
extract_parameter <- function(parameter_name) {
  index <- match(parameter_name, args)
  if (!is.na(index) && index < length(args)) {
    return(args[index + 1])
  } else {
    return("")  # Default value if not provided
  }
}

# Extract parameter values
name <- extract_parameter("--name")
file1 <- extract_parameter("--file1")
file2 <- extract_parameter("--file2")
file3 <- extract_parameter("--file3")
file4 <- extract_parameter("--file4")
out1 <- extract_parameter("--out1")
out2 <- extract_parameter("--out2")

cat("# R parameters:",name,file1,file2,file3,file4,out1,out2,sep="\n")

sam1C1 <- read.table(file1, sep = "\t", header = F, stringsAsFactors = FALSE)
sam1C2 <- read.table(file2, sep = "\t", header = F, stringsAsFactors = FALSE)
sam2C1 <- read.table(file3, sep = "\t", header = F, stringsAsFactors = FALSE)
sam2C2 <- read.table(file4, sep = "\t", header = F, stringsAsFactors = FALSE)

cat("# head sam1C1 \n")
head(sam1C1)
cat("# head sam1C2 \n")
head(sam1C2)
cat("# head sam2C1 \n")
head(sam2C1)
cat("# head sam2C2 \n")
head(sam2C1)


m1 <- merge(sam1C1, sam1C2, by = c("V1", "V2","V3","V4","V5"))
m1 <- merge(m1, sam2C1, by = c("V1", "V2","V3","V4","V5"))
m1 <- merge(m1, sam2C2, by = c("V1", "V2","V3","V4","V5"))
colnames(m1) <- c("chr","start","end","gv","bin","sam1C1_Meth","sam1C1_Dept","sam1C2_Meth","sam1C2_Dept","sam2C1_Meth","sam2C1_Dept","sam2C2_Meth","sam2C2_Dept")

cat("# head merged \n")
head(m1)
write.table(m1, out2, row.names = F, col.names = T, sep = "\t",quote =F)
cat ("# R merge finished! \n")
m1 <- data.frame(m1)

gv <- m1[1,4]
bin <- m1[1,5]

m1$meth_diff <- round((m1$sam1C1_Meth + m1$sam1C2_Meth - m1$sam2C1_Meth - m1$sam2C2_Meth)/2,4)
m1$dept_diff <- round((m1$sam1C1_Dept + m1$sam1C2_Dept - m1$sam2C1_Dept - m1$sam2C2_Dept)/2,4)

cat("# head merged with meth-dif and dept-dif \n")
head(m1)

me_df <- round(mean(m1$meth_diff),4)
dp_df <- round(mean(m1$dept_diff),4)

final <- paste(name,gv,bin,me_df,dp_df,sep="\t")
write.table(final, out1, row.names = F, col.names = F, sep = "\t",quote =F)
cat("# final diff mean results \n",final,"\n R finished")
# R code stop here
' --name ${name}  --file1 ${out2_1} --file2 ${out2_2} --file3 ${out2_3} --file4 ${out2_4} --out2 ${out2_5} --out1 ${out2_5}_tem

######################################################################################################################################################################

echo -e "step 2 done, output:\n${out2_1}\n${out2_2}\n${out2_3}\n${out2_4}\n${out2_5}\n${out1}"
date
echo "#========================================================= step 3: calculate methylation difference and mean ============================================="
out3="${results_dir}/step3_${name}_means.txt"

cat ${out2_5}_tem >> ${out3}_tem
sort -k 2 -n ${out3}_tem  > ${out3}

rm ${out2_1} ${out2_2} ${out2_3} ${out2_4} ${out2_5}_tem ${out3}_tem ${out1}_tem 

echo -e "step 3 done, output:\n${out3}"
echo -e "all done!"

