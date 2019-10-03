
#Star alignment in sudo's computer 
# can use  "source activate aligners" to activate STAR

mkdir $output_path
genomedir=~/Apps/STAR
dataprefix=Cecilia_RNA #prefix for the data sets
rawfile_path=~ # file path for the rawfile 
output_path=~/star_daokun_ouput  #general file path for the STAR agliment output
refbed_path=~/Rseqc_ref #Reference file path for RseQC
output_path1=${output_path}"/File_Based_QC"        # file path for input OC
output_path2=${output_path}"/Sample_Based_QC"   # file path for sample QC
cd $rawfile_path
for dir in ${dataprefix}*
do 
	cd $dir
	mkdir $output_path/$dir
	for file in *
	do 
		cd $file 
		i=1
		for f1 in *_R1.fastq.gz
		do
			f2=${f1%%R1.fastq.gz}"R2.fastq.gz"
			f3=${f1%%_P0001_R1.fastq.gz}
			f4=${f1%%_L*}
			~/Apps/anaconda/anaconda3/bin/STAR --genomeDir $genomedir --runThreadN 8 --readFilesCommand zcat --readFilesIn $f1 $f2 --outFileNamePrefix ${f3}_ --outSAMtype BAM SortedByCoordinate           			
			if [ $i == 1 ]
		    then 
		        mkdir $output_path/$dir/$f4
            fi		
		    i=$i+1				
			mkdir $output_path/$dir/$f4/$f3
			for f in *out*
			do
				mv $f $output_path/$dir/$f4/$f3
			done
		done
		cd ..		
	done 
	cd ..
done


#after transfer to my computer do the QC in using hpcc

export PATH=/shared/app/anaconda2/bin/:$PATH
export LD_LIBRARY_PATH=/shared/app/anaconda2/lib:$LD_LIBRARY_PATH
export PATH=/shared/app/SAMTOOLS:$PATH
export PATH=/shared/app/FastQC:$PATH
#define the file path. 



# mkdir $rawfile_path if needed
#do fastqc for each input using for loop 
cd $rawfile_path
for dir in ${dataprefix}*
do
	cd $dir
	mkdir rawdata_fastqcoutput
	for file in *-*
	do 
		cd $file 
		for f in *.fastq.gz
		do
			fastqc ${f} --outdir=$rawfile_path/$dir/rawdata_fastqcoutput --thread=8 
		done
		cd ..
	done
	cd rawdata_fastqcoutput 
	rm *.zip
	cd ../..
done

---------------------------------------------------------------
#copy the 9 fastqc output into 9 different folders
cd $output_path
mkdir fastqc
cd fastqc
mkdir duplication_levels
mkdir kmer_profiles
mkdir per_base_gc_content
mkdir per_base_n_content
mkdir per_base_quality
mkdir per_base_sequence_content
mkdir per_sequence_gc_content
mkdir per_sequence_quality
mkdir sequence_length_distribution
cd $rawfile_path
for dir in ${dataprefix}*
do 
	cd $dir 
	cd rawdata_fastqcoutput_repeat
		for file in *_L*
		do
			cd $file
			cd Images
			cp duplication_levels.png $output_path/fastqc/duplication_levels/${file%%P00*}"duplication_levels.png"
			cp kmer_profiles.png $output_path/fastqc/kmer_profiles/${file%%P00*}"kmer_profiles.png" # wrong
			cp per_base_gc_content.png $output_path/fastqc/per_base_gc_content/${file%%P00*}"per_base_gc_content.png"
			cp per_base_n_content.png $output_path/fastqc/per_base_n_content/${file%%P00*}"per_base_n_content.png"
			cp per_base_quality.png $output_path/fastqc/per_base_quality/${file%%P00*}"per_base_quality.png"
			cp per_base_sequence_content.png $output_path/fastqc/per_base_sequence_content/${file%%P00*}"per_base_sequence_content.png"
			cp per_sequence_gc_content.png $output_path/fastqc/per_sequence_gc_content/${file%%P00*}"per_sequence_gc_content.png"
			cp per_sequence_quality.png $output_path/fastqc/per_sequence_quality/${file%%P00*}"per_sequence_quality.png"
			cp sequence_length_distribution.png $output_path/fastqc/sequence_length_distribution/${file%%P00*}"sequence_length_distribution.png"
			cd ../..
		done
	cd ../..	 		   
done


----------------------------------------------------------------



cd $output_path
for dir in ${dataprefix}*
do 
	cd $dir
	for file in *_S*
	do
		cd $file
		for f in *_L*
		do 
			cd $f
			
			#Do insert_size and duplication in picard in a for loop
			# I is input O is output prefix
			for bam in *.bam
			do 
				java -jar /shared/app/picard-tools-1/picard-tools-1.96/CollectInsertSizeMetrics.jar I=$bam O=${f}_insert_size.xls H=${f}_insert_size_histogram.pdf M=0.5

				java -jar /shared/app/picard-tools-1/picard-tools-1.96/EstimateLibraryComplexity.jar I=$bam O=${f}_complex_metrics.xls
			done
			
			#Do QC using Rseqc for each input using for loop: 1, read dupliation 2, junction saturation 3, read distribution 
			# 4, juntion annotaion 5, ribosome contamination   
			# -i: input files   
			# -o: output prefix
			read_duplication.py -i *.bam -o ${f}
			junction_saturation.py -i *.bam -r ${refbed_path}/mm10_RefSeq.bed -o ${f}
			read_distribution.py  -i *.bam -r ${refbed_path}/mm10_RefSeq.bed >${f}"_read_distribution".txt
			junction_annotation.py -i *.bam -o ${f} -r  ${refbed_path}/mm10_RefSeq.bed >>${f}"_junction".txt 2>&1 # 2>&1 means put error output with screen output
			infer_experiment.py -r ${refbed_path}/mm10_RefSeq.bed -i *.bam >${f}"_strandness".txt
			split_bam.py -i *.bam  -r ${refbed_path}/*rRNA.bed -o $f > ${f}"_ribosome".txt
			rm *.in.bam
			rm *.ex.bam
			rm *.junk.bam
			samtools index *.bam ${f}"_Aligned.sortedByCoord.out.bam.bai"
			cd ..
		done
		cd ..
	done
	cd ..
done




# generate the QC output and the files will be read into R and generate the final QC output

mkdir $output_path1
for dir in ${dataprefix}*
do
	cd $dir
	for file in *_S*
	do
		cd $file
		for f in *L*
		do 
			cd $f
			#1 for ribosome contamination
			echo ${f}>>$output_path1/ribosome_names.csv
			sed '2!d' *_ribosome.txt  | sed "s/${f}.in.bam (Reads consumed by input gene list)://g" >> $output_path1/ribosome_rRNA.csv
			sed '3!d' *_ribosome.txt  | sed "s/${f}.ex.bam (Reads not consumed by input gene list)://g" >>$output_path1/ribosome_gene.csv

			#2 strandness
			echo ${f}>>$output_path1/strandness_1.csv
			sed '5!d' *_strandness.txt  | sed "s/Fraction of reads explained by \"1++,1--,2+-,2-+\": //g" >> $output_path1/strandness_2.csv
			sed '6!d' *_strandness.txt  | sed "s/Fraction of reads explained by \"1+-,1-+,2++,2--\": //g" >> $output_path1/strandness_3.csv
			sed '4!d' *_strandness.txt  | sed "s/Fraction of reads failed to determine: //g" >> $output_path1/strandness_4.csv 

			#3 read distribution
			echo ${f}>>$output_path1/read_distribution_names.csv
			awk 'FNR == 8  { print $3}' *read_distribution.txt>>$output_path1/3_UTR.csv
			awk 'FNR == 7  { print $3}' *read_distribution.txt>>$output_path1/5_UTR.csv
			awk 'FNR == 6  { print $3}' *read_distribution.txt>>$output_path1/CDS_.csv
			awk 'FNR == 9  { print $3}' *read_distribution.txt>>$output_path1/Intron_.csv
			awk 'FNR == 12 { print $3}' *read_distribution.txt>>$output_path1/TSS_10kb.csv
			awk 'FNR == 15 { print $3}' *read_distribution.txt>>$output_path1/TES_10kb.csv
			awk 'FNR == 1  { print $3}' *read_distribution.txt>>$output_path1/Total_Reads.csv
			awk 'FNR == 2  { print $3}' *read_distribution.txt>>$output_path1/Total_tags.csv
			awk 'FNR == 3  { print $4}' *read_distribution.txt>>$output_path1/Total_assigned_tags.csv

			#4 mapping statistics for every input and mapping statistics #for every sample(combine replicates)
			echo ${f}>>$output_path1/mapping_names.csv
			echo ${file}>>$output_path1/Sample_mapping_names.csv
			awk 'FNR == 6  { print $6}' *Log.final.out>>$output_path1/number_of_input_reads.csv
			awk 'FNR == 9  { print $6}' *Log.final.out>>$output_path1/number_of_unique_reads.csv
			awk 'FNR == 24  { print $9}' *Log.final.out>>$output_path1/number_of_multiple_reads.csv

			#5 junction
			echo ${f}>>$output_path1/junction_names.csv
			awk 'FNR == 5  { print $4}' *_junction.txt>>$output_path1/Total_Splicing_Events.csv
			awk 'FNR == 6  { print $4}' *_junction.txt>>$output_path1/Known_Splicing_Events.csv
			awk 'FNR == 7  { print $5}' *_junction.txt>>$output_path1/Partial_Novel_Splicing_Events.csv
			awk 'FNR == 8  { print $4}' *_junction.txt>>$output_path1/Novel_Splicing_Events.csv
			awk 'FNR == 10  { print $4}' *_junction.txt>>$output_path1/Total_Splicing_Junctions.csv
			awk 'FNR == 11  { print $4}' *_junction.txt>>$output_path1/Known_Splicing_Junctions.csv
			awk 'FNR == 12  { print $5}' *_junction.txt>>$output_path1/Partial_Novel_Splicing_Junctions.csv
			awk 'FNR == 13  { print $4}' *_junction.txt>>$output_path1/Novel_Splicing_Junctions.csv

			cd ..
		done
		cd ..
	done
	cd ..
done


#6 insert size
i=1
for dir in ${dataprefix}*
do
	cd $dir
	for file in *_S*
	do
		cd $file
		for f in *L*
		do 
			cd $f
			if [ $i == 1 ]
			then 
				echo "Sample_Name">>$output_path1/insert_size_names.txt #generate the header
				awk 'FNR == 7  {$19=$20=$21=""; print $0}' *insert_size.xls>>$output_path1/insert_size.txt  #generate the header
			fi
			i=$i+1
			echo ${f}>>$output_path1/insert_size_names.txt
			sed '8!d' *insert_size.xls>>$output_path1/insert_size.txt
			cd ..
		done
		cd ..
	done
	cd ..
done

#7 duplication rate
i=1
for dir in ${dataprefix}*
do
	cd $dir
	for file in *_S*
	do
		cd $file
		for f in *L*
		do 
			cd $f
			if [ $i == 1 ]
			then
				echo "Sample_Name">>$output_path1/duplication_names.txt
				sed '7!d' *complex_metrics.xls>>$output_path1/duplication_.txt
			fi
			i=$i+1
			echo ${f}>>$output_path1/duplication_names.txt
			sed '8!d' *complex_metrics.xls>>$output_path1/duplication_.txt
			cd ..
		done
		cd ..
	done
	cd ..
done

# for the Rscript,need to set the working dir as $output_path1
Rscript $output_path/file_mapping_stats.R
rm $output_path1/*_*.csv
rm $output_path1/*_*.txt


#--------------------------------------------
# generate the txt file for samtools merge and use -b in samtools
# to merge the sam files from same sample
for dir in ${dataprefix}*
do 
	cd $dir
	for file in *_S*
	do
		cd $file
		for f in *_L*
		do
			cd $f
			for bam in *.bam
			do
				echo $output_path/$dir/$file/$f/$bam >> $output_path/${file}.txt
			done
			cd ..
		done
		cd ..
	done
	cd ..
done

#using samtools to merge the input files from same sample
export PATH=/shared/app/samtools-1.0/bin:$PATH

for dir in ${dataprefix}*
do 
mkdir ${dir}"_merge"
cd $dir
for file in *_S*
do
cd $file
mkdir $output_path/${dir}"_merge"/$file
samtools merge $output_path/${dir}"_merge"/$file/${file%%_S*}.bam -b $output_path/${file}.txt
cd ..
done
cd ..
done

rm $output_path/*_*.txt


#do QC using Rseqc for merged bam files 
#Junction saturation, read distribution,junction annotaion
#strandness,ribosome contamination,duplication,insert size
#and samtool index

export PATH=/shared/app/SAMTOOLS:$PATH

for dir in ${dataprefix}*merge
do 
	cd $dir
	for f in *
	do 
		cd $f
		for bam in *.bam
		do 
			java -jar /shared/app/picard-tools-1/picard-tools-1.96/CollectInsertSizeMetrics.jar I=$bam O=${f%%_S*}_insert_size.xls H=${f%%_S*}_insert_size_histogram.pdf M=0.5

			java -jar /shared/app/picard-tools-1/picard-tools-1.96/EstimateLibraryComplexity.jar I=$bam O=${f%%_S*}_complex_metrics.xls
		done
		junction_saturation.py -i *.bam -r ${refbed_path}/mm10_RefSeq.bed -o ${f%%_S*}
		read_distribution.py  -i *.bam -r ${refbed_path}/mm10_RefSeq.bed >${f%%_S*}"_read_distribution".txt
		junction_annotation.py -i *.bam -o ${f%%_S*} -r  ${refbed_path}/mm10_RefSeq.bed >>${f%%_S*}"_junction".txt 2>&1 
		infer_experiment.py -r ${refbed_path}/mm10_RefSeq.bed -i *.bam >${f%%_S*}"_strandness".txt
		split_bam.py -i *.bam  -r ${refbed_path}/*rRNA.bed -o $f > ${f%%_S*}"_ribosome".txt
		rm *.in.bam
		rm *.ex.bam
		rm *.junk.bam
		samtools index *.bam ${f%%_S*}".bam.bai"
		cd ..
	done
	cd ..
done



# do ribosome contamination, strandness, duplication rate, 
#insert size,read distribution and junction for each sample


mkdir $output_path2
for dir in ${dataprefix}*merge
do
	cd $dir
	for file in *_S*
	do
		cd $file

		#1 for ribosome contamination
		echo ${file%%_S*}>>$output_path2/ribosome_names.csv
		sed '2!d' *_ribosome.txt  | sed "s/${file}.in.bam (Reads consumed by input gene list)://g" >> $output_path2/ribosome_rRNA.csv
		sed '3!d' *_ribosome.txt  | sed "s/${file}.ex.bam (Reads not consumed by input gene list)://g" >>$output_path2/ribosome_gene.csv

		#2 strandness
		echo ${file%%_S*}>>$output_path2/strandness_1.csv
		sed '5!d' *_strandness.txt  | sed "s/Fraction of reads explained by \"1++,1--,2+-,2-+\": //g" >> $output_path2/strandness_2.csv
		sed '6!d' *_strandness.txt  | sed "s/Fraction of reads explained by \"1+-,1-+,2++,2--\": //g" >> $output_path2/strandness_3.csv
		sed '4!d' *_strandness.txt  | sed "s/Fraction of reads failed to determine: //g" >> $output_path2/strandness_4.csv 

		#3 read distribution
		echo ${file%%_S*}>>$output_path2/read_distribution_names.csv
		awk 'FNR == 8  { print $3}' *read_distribution.txt>>$output_path2/3_UTR.csv
		awk 'FNR == 7  { print $3}' *read_distribution.txt>>$output_path2/5_UTR.csv
		awk 'FNR == 6  { print $3}' *read_distribution.txt>>$output_path2/CDS_.csv
		awk 'FNR == 9  { print $3}' *read_distribution.txt>>$output_path2/Intron_.csv
		awk 'FNR == 12 { print $3}' *read_distribution.txt>>$output_path2/TSS_10kb.csv
		awk 'FNR == 15 { print $3}' *read_distribution.txt>>$output_path2/TES_10kb.csv
		awk 'FNR == 1  { print $3}' *read_distribution.txt>>$output_path2/Total_Reads.csv
		awk 'FNR == 2  { print $3}' *read_distribution.txt>>$output_path2/Total_tags.csv
		awk 'FNR == 3  { print $4}' *read_distribution.txt>>$output_path2/Total_assigned_tags.csv


		#4 junction
		echo ${file%%_S*}>>$output_path2/junction_names.csv
		awk 'FNR == 5  { print $4}' *_junction.txt>>$output_path2/Total_Splicing_Events.csv
		awk 'FNR == 6  { print $4}' *_junction.txt>>$output_path2/Known_Splicing_Events.csv
		awk 'FNR == 7  { print $5}' *_junction.txt>>$output_path2/Partial_Novel_Splicing_Events.csv
		awk 'FNR == 8  { print $4}' *_junction.txt>>$output_path2/Novel_Splicing_Events.csv
		awk 'FNR == 10  { print $4}' *_junction.txt>>$output_path2/Total_Splicing_Junctions.csv
		awk 'FNR == 11  { print $4}' *_junction.txt>>$output_path2/Known_Splicing_Junctions.csv
		awk 'FNR == 12  { print $5}' *_junction.txt>>$output_path2/Partial_Novel_Splicing_Junctions.csv
		awk 'FNR == 13  { print $4}' *_junction.txt>>$output_path2/Novel_Splicing_Junctions.csv
           # awk 'FNR == 13,14  { print $4,$5}' can print 13,14 row and 4,5 column
		cd ..
	done
	cd ..
done


#5 insert size
i=1
for dir in ${dataprefix}*merge
do
	cd $dir
	for file in *_S*
	do
		cd $file
		if [ $i == 1 ]
		then 
			echo "Sample_Name">>$output_path2/insert_size_names.txt
			awk 'FNR == 7  {$19=$20=$21=""; print $0}' *insert_size.xls>>$output_path2/insert_size.txt
		fi
		i=$i+1
		echo ${file%%_S*}>>$output_path2/insert_size_names.txt
		sed '8!d' *insert_size.xls>>$output_path2/insert_size.txt
		cd ..
	done
	cd ..
done

#6 duplication rate
i=1
for dir in ${dataprefix}*merge
do
	cd $dir
	for file in *_S*
	do
		cd $file
		if [ $i == 1 ]
		then 
			echo "Sample_Name">>$output_path2/duplication_names.txt
			sed '7!d' *complex_metrics.xls>>$output_path2/duplication_.txt
		fi
		i=$i+1
		echo ${file%%_S*}>>$output_path2/duplication_names.txt
		sed '8!d' *complex_metrics.xls>>$output_path2/duplication_.txt
		cd ..
	done
	cd ..
done

#7 junction saturation
for dir in ${dataprefix}*merge
do
	cd $dir
	for file in *_S*
	do
		cd $file
		echo ${file%%_S*} ${file%%_S*}"_all_junctions">>$output_path2/junction_saturation_names.txt
		echo ${file%%_S*} ${file%%_S*}"_known_junctions">>$output_path2/junction_saturation_names.txt
		echo ${file%%_S*} ${file%%_S*}"_novel_junctions">>$output_path2/junction_saturation_names.txt
		sed '3!d' *junctionSaturation_plot.r | sed "s/y=c(//g" | sed "s/)//g" | sed "s/,/ /g">>$output_path2/junction_saturation.txt
		sed '4!d' *junctionSaturation_plot.r | sed "s/z=c(//g" | sed "s/)//g" | sed "s/,/ /g">>$output_path2/junction_saturation.txt
		sed '5!d' *junctionSaturation_plot.r | sed "s/w=c(//g" | sed "s/)//g" | sed "s/,/ /g">>$output_path2/junction_saturation.txt
		cd ..
	done
	cd ..
done


# for the Rscript,need to set the working dir as $output_path2
Rscript $output_path/sample_mapping_stats.R
rm $output_path2/*_*.csv
rm $output_path2/*_*.txt

mv $output_path1/mapping*rm-replicate*.csv $output_path2

#copy junction saturation to one folder
mkdir $output_path/junction_saturation
for dir in ${dataprefix}*merge
do
cd $dir
for file in *_S*
do
cd $file
cp *junctionSaturation*pdf $output_path/junction_saturation
cd ..
done
cd ..
done


# geneBody_coverage for each sample(samples in one set will be put in one geneBody_coverage output)
for dir in ${dataprefix}*merge
do
	cd $dir
	echo "$ cat bam_path.txt">>bam_path.txt
	for file in *_S*
	do
		cd $file
		echo $output_path/$dir/$file/*.bam>>$output_path/$dir/bam_path.txt 
		cd ..
	done
	geneBody_coverage.py -r $refbed_path/mm10.House*.bed -i bam_path.txt  -o $output_path2/${dir%%_merge}
	cd ..
done



#use HTseq-count to get the raw counts
gtf_path=~/star_daokun_ouput

#htseq-counts with reverse for set1 and none strandness options for set2,3 # finally use -r pos do continue
i=1
for dir in ${dataprefix}*merge
do
cd $dir
if [ $i == 1 ]
then 
for file in *_S*
do
cd $file 
/shared/app/anaconda2/bin/htseq-count -f bam -r pos --stranded=reverse -t exon -i gene_id --additional-attr=gene_name *.bam  $gtf_path/gencode.v*.gtf > ${file%%_S*}"_counts.txt"
cd ..
done
else 
for file in *_S*
do
cd $file 
/shared/app/anaconda2/bin/htseq-count -f bam -r pos --stranded=no -t exon -i gene_id --additional-attr=gene_name *.bam  $gtf_path/gencode.v*.gtf > ${file%%_S*}"_counts.txt"
cd ..
done		
fi
i=$i+1	
cd ..
done


#generate the excel file for all the raw counts
i=1
for dir in ${dataprefix}*merge
do
	cd $dir
	for file in *_S*
	do
		cd $file
		if [ $i == 1 ]
		then 
			awk '{ print $1}' *counts.txt>>$output_path/raw_counts_gene_id.txt
			awk '{ print $2}' *counts.txt>>$output_path/raw_counts_gene_names.txt
		fi
		i=$i+1
		awk '{ print $3}' *counts.txt>>$output_path/${file%%_S*}_raw_counts.txt
		awk '{ print $2}' *counts.txt>>$output_path/${file%%_S*}_raw_counts_summary.txt
		cd ..
	done
	cd ..
done

# it is very important to do file names in this way since the alphabetic order of file
i=1
for f in *_raw_counts.txt
do 
if [ $i == 1 ]
then 
echo "Gene_ID">>$output_path/raw_counts_names.txt
echo "Gene_Name">>$output_path/raw_counts_names.txt
fi
i=$i+1
echo ${f%%_raw_counts.txt}>>raw_counts_names.txt
done

Rscript $output_path/raw_counts.R
rm *raw_counts*.txt
# end of the script here 




#--------------------------------------------------------------------------------------------------------------------------
#code from thierry
/shared/app/anaconda2/bin/htseq-count -f bam -r name --stranded=reverse -t exon -i gene_id --additional-attr=gene_name   *.bam  ~/dokun/gencode.v*.gtf 

#try thierry's new gtf file
for dir in *Set1*merge
do
cd $dir 
for file in *_S*
do
cd $file 
/shared/app/anaconda2/bin/htseq-count -f bam -r pos --stranded=reverse -t exon -i gene_id --additional-attr=gene_name *.bam  $gtf_path/Mus*E2f3ab.gtf > ${file%%_S*}"_counts.txt"
cd ..
done
cd ..
done


/shared/app/anaconda2/bin/htseq-count -f bam -r pos --stranded=reverse -t exon -i gene_id --additional-attr=gene_name *.bam  $gtf_path/gencode.*E2f3ab.gtf > E3fabcounts.txt

#thiery's way to get E2F3ab what we use to do 4.17.18
		samtools view *.bam "chr13:29985300-29986063"|wc >>$output_path/E2F3ab.txt 
		samtools view *.bam "chr13:29984128-29984391"|wc >>$output_path/E2F3ab.txt


# get counts for the set1 RNAseq data without merge files
gtf_path=~/star_daokun_ouput

for file in 8-*_S*
do
cd $file 
for f in *L*
do 
cd $f
/shared/app/anaconda2/bin/htseq-count -f bam -r pos --stranded=reverse -t exon -i gene_id --additional-attr=gene_name *.bam  $gtf_path/gencode.v*tion.gtf > ${f}"_counts.txt"
cd ..
done
cd ..
done			

for file in 5_*_S*
do
cd $file 
for f in *L*
do 
cd $f
/shared/app/anaconda2/bin/htseq-count -f bam -r pos --stranded=reverse -t exon -i gene_id --additional-attr=gene_name *.bam  $gtf_path/gencode.v*tion.gtf > ${f}"_counts.txt"
cd ..
done
cd ..
done

i=1
for dir in *_S*
do
	cd $dir
	for file in *L*
	do
		cd $file
		if [ $i == 1 ]
		then 			
			awk '{ print $1}' *counts.txt>>$output_path/raw_counts_gene_id.txt			
			awk '{ print $2}' *counts.txt>>$output_path/raw_counts_gene_names.txt
		fi
		i=$i+1
		awk '{ print $3}' *counts.txt>>$output_path/${file}_raw_counts.txt
		awk '{ print $2}' *counts.txt>>$output_path/${file}_raw_counts_summary.txt
		cd ..
	done
	cd ..
done
# it is very important to do in this way since the alphabetic order of file
i=1
for f in *_raw_counts.txt
do 
if [ $i == 1 ]
then 
echo "Gene_ID">>$output_path/raw_counts_names.txt
echo "Gene_Name">>$output_path/raw_counts_names.txt
fi
i=$i+1
echo ${f%%_raw_counts.txt}>>raw_counts_names.txt
done

Rscript $output_path/raw_counts.R
rm *raw_counts*.txt


#------------------------------------------------------------------
#get counts for set 2 strandness is no
for dir in Cecilia*Set*2017
do
cd $dir
for file in *_S*
do
cd $file 
for f in *L*
do 
cd $f
/shared/app/anaconda2/bin/htseq-count -f bam -r pos --stranded=no -t exon -i gene_id --additional-attr=gene_name *.bam  $gtf_path/gencode.v*tion.gtf > ${f}"_counts.txt"
cd ..
done
cd ..
done
cd ..
done

#get E2F3a and 3b relative proportion
samtools view *.bam "chr13:29985300-29986063"|wc -l >>$output_path/E2F3ab.txt 
samtools view *.bam "chr13:29984128-29984391"|wc -l >>$output_path/E2F3ab.txt

i=1
for f in Cecilia*2017
do 
cd $f
for dir in *_S*
do
	cd $dir
	for file in *L*
	do
		cd $file
		if [ $i == 1 ]
		then 
			awk '{ print $1}' *counts.txt>>$output_path/raw_counts_gene_id.txt
			awk '{ print $2}' *counts.txt>>$output_path/raw_counts_gene_names.txt
		fi
		i=$i+1
		awk '{ print $3}' *counts.txt>>$output_path/${file}_raw_counts.txt
		awk '{ print $2}' *counts.txt>>$output_path/${file}_raw_counts_summary.txt
		cd ..
	done
	cd ..
done
cd ..
done

# it is very important to do in this way since the alphabetic order of file
i=1
for f in *_raw_counts.txt
do 
if [ $i == 1 ]
then 
echo "Gene_ID">>$output_path/raw_counts_names.txt
echo "Gene_Name">>$output_path/raw_counts_names.txt
fi
i=$i+1
echo ${f%%_raw_counts.txt}>>raw_counts_names.txt
done
Rscript $output_path/raw_counts.R
rm *raw_counts*.txt


# get E2F3a and 3b counts seperately
for dir in Cecilia*2017
do 
cd $dir
for file in *_S*
do
cd $file
for f in *L*
do
cd $f
echo ${f%%_Aligned*}>>$output_path/file_names.txt
samtools view *.bam "chr13:29985300-29986063"|wc -l >>$output_path/E2F3ab.txt 
samtools view *.bam "chr13:29984128-29984391"|wc -l >>$output_path/E2F3ab.txt
cd ..
done
cd ..
done
cd ..
done


#get E2F3 counts from all gene counts
grep E2f3 raw-counts.csv | sed '1!d' | sed "s/\"1721\",\"ENSMUSG00000016477.18\",\"E2f3\",//g" | sed "s/,/ /g" > E2F3_counts.txt


