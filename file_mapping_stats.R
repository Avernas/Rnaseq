setwd ("/home/das260/star_daokun_ouput/File_Based_QC")

#ribosome contamination
R1=read.csv("ribosome_names.csv",header=F)                                     
R2=read.csv("ribosome_rRNA.csv",header=F)   
R3=read.csv("ribosome_gene.csv",header=F)
R4=R2/(R2+R3)
R5=cbind(R1,R2,R3,R4)
colnames(R5)=c("Sample_Name","Reads_consumed_by_input_gene_list","Reads_not_consumed_by_input_gene_list","Ribosomal_percentage")
write.csv(R5, file="ribosome-contamination.csv")

#strandness
S1=read.csv("strandness_1.csv",header=F)                                     
S2=read.csv("strandness_2.csv",header=F)   
S3=read.csv("strandness_3.csv",header=F)
S4=read.csv("strandness_4.csv",header=F)
S5=cbind(S1,S2,S3,S4)
colnames(S5)=c("Sample_Name","Fraction_of_reads_explained_by_1++,1--,2+-,2-+","Fraction_of_reads_explained_by_1+-,1-+,2++,2--","Fraction_of_reads_failed_to_determine")
write.csv(S5, file="strandness.csv")

#Read Distribution
Reads_Name=read.csv("read_distribution_names.csv",header=F)
UTR3=read.csv("3_UTR.csv",header=F)                                     
UTR5=read.csv("5_UTR.csv",header=F)   
CDS=read.csv("CDS_.csv",header=F)
Introns=read.csv("Intron_.csv",header=F)
TSS=read.csv("TSS_10kb.csv",header=F)
TES=read.csv("TES_10kb.csv",header=F)
Total_Reads=read.csv("Total_Reads.csv",header=F)
Total_tags=read.csv("Total_tags.csv",header=F)
Total_Assigned_Tags=read.csv("Total_assigned_tags.csv",header=F)
Unassigned_Tags=Total_tags-UTR3-UTR5-CDS-Introns-TSS-TES
exon=(UTR3+UTR5+CDS)/Total_tags
intron=Introns/Total_tags
intergenic=1-exon-intron
RE1=cbind(Reads_Name,UTR3,UTR5,CDS,Introns,TES,TSS,Unassigned_Tags,Total_Reads,Total_Assigned_Tags,Total_tags,Reads_Name,exon,intron,intergenic)
colnames(RE1)=c("Sample_Name","3'UTR","5'UTR","CDS","Introns","TES_10kb","TSS_10kb","Unassigned_Tags","Total_Reads","Total_Assigned_Tags","Total_Tags","Sample_Name","exon","intron","intergenic")
write.csv(RE1, file="read-distribution.csv")

#insert size 
insert_size_name=read.table("insert_size_names.txt",header=T)
insert_size=read.table("insert_size.txt",header=T)
I1=cbind(insert_size_name,insert_size)
write.csv(I1, file="insert-size.csv")

#duplication
duplication_name=read.table("duplication_names.txt",header=T)
duplication=read.table("duplication_.txt",header=T)
attach(duplication)
D1=cbind(duplication_name,READ_PAIRS_EXAMINED,READ_PAIR_DUPLICATES,READ_PAIR_OPTICAL_DUPLICATES,PERCENT_DUPLICATION,ESTIMATED_LIBRARY_SIZE)
write.csv(D1, file="duplication-rate.csv")
detach(duplication)

#junction
Junction_Name=read.csv("junction_names.csv",header=F)
Total_Splicing_Events=read.csv("Total_Splicing_Events.csv",header=F)
Known_Splicing_Events=read.csv("Known_Splicing_Events.csv",header=F)
Partial_Novel_Splicing_Events=read.csv("Partial_Novel_Splicing_Events.csv",header=F)
Novel_Splicing_Events=read.csv("Novel_Splicing_Events.csv",header=F)
Total_Splicing_Junctions=read.csv("Total_Splicing_Junctions.csv",header=F)
Known_Splicing_Junctions=read.csv("Known_Splicing_Junctions.csv",header=F)
Partial_Novel_Splicing_Junctions=read.csv("Partial_Novel_Splicing_Junctions.csv",header=F)
Novel_Splicing_Junctions=read.csv("Novel_Splicing_Junctions.csv",header=F)
J1=cbind(Junction_Name,Total_Splicing_Events,Known_Splicing_Events,Partial_Novel_Splicing_Events,Novel_Splicing_Events,Junction_Name,Total_Splicing_Junctions,Known_Splicing_Junctions,Partial_Novel_Splicing_Junctions,Novel_Splicing_Junctions)
colnames(J1)=c("Junction_Name","Total_Splicing_Events","Known_Splicing_Events","Partial_Novel_Splicing_Events","Novel_Splicing_Events","Junction_Name","Total_Splicing_Junctions","Known_Splicing_Junctions","Partial_Novel_Splicing_Junctions","Novel_Splicing_Junctions")
write.csv(J1, file="junction.csv")

#mapping statistics for every input
Sample_Mapping_Name=read.csv("Sample_mapping_names.csv",header=F)
mapping_names=read.csv("mapping_names.csv",header=F)                                     
input_reads=read.csv("number_of_input_reads.csv",header=F)   
unique_reads=read.csv("number_of_unique_reads.csv",header=F)
multiple_reads=read.csv("number_of_multiple_reads.csv",header=F)
unmapped_reads=input_reads-unique_reads-multiple_reads
unique_ratio=unique_reads/input_reads
multiple_ratio=multiple_reads/input_reads
unmapped_ratio=unmapped_reads/input_reads
M1=cbind(Sample_Mapping_Name,mapping_names,input_reads,unique_reads,multiple_reads,unmapped_reads,unique_ratio,multiple_ratio,unmapped_ratio)
colnames(M1)=c("Sample_Mapping_Names","Sample_Name","Total_Input_Reads","Unique_Reads","Multiple_Reads","Unmapped_Reads","Unique_Mapping_Ratio","Multiple_Mapping_Ratio","Unmapped_Ratio")
write.csv(M1, file="mapping-statistics.csv")

#mapping statistics for every sample(combine replicates)

M1=read.csv("mapping-statistics.csv",header=T) 
Sample_index=array()
Sample_Unique_Name=c()
Sample_Total_Read=array()
Sample_Unique_Read=array()
Sample_Multiple_Read=array()
Sample_Unmapped_Read=array()
i=1
attach(M1)
for (index in M1$X)
{
  
  if ((Sample_Mapping_Names[X==index] != Sample_Mapping_Names[X==index+1]) && index < max(M1$X) )# create flags for each sample by checking the names for each input
  {
    Sample_index[index]=i
    Sample_Unique_Name[i]=toString(Sample_Mapping_Names[X==index])
    i=i+1
  }
  else if (index == max(M1$X)) # the last input need to have a sample unique name so add it here
  {
    Sample_index[index]=i
    Sample_Unique_Name[i]=toString(Sample_Mapping_Names[X==index])
  }
  else 
  {
    Sample_index[index]=i
  }
}
detach(M1)
M1=cbind(M1,Sample_index)
attach(M1)
for (ii in 1:length(Sample_Unique_Name))
{
  Sample_Total_Read[ii]=sum(Total_Input_Reads[M1$Sample_index==ii])
  Sample_Unique_Read[ii]=sum(Unique_Reads[M1$Sample_index==ii])
  Sample_Multiple_Read[ii]=sum(Multiple_Reads[M1$Sample_index==ii])
  Sample_Unmapped_Read[ii]=sum(Unmapped_Reads[M1$Sample_index==ii])
}
Sample_Unique_Mapping_Ratio=Sample_Unique_Read/Sample_Total_Read
Sample_Multiple_Mapping_Ratio=Sample_Multiple_Read/Sample_Total_Read
Sample_Unmapped_Mapping_Ratio=Sample_Unmapped_Read/Sample_Total_Read
M2=cbind(Sample_Unique_Name,Sample_Total_Read,Sample_Unique_Read,Sample_Multiple_Read,Sample_Unmapped_Read,Sample_Unique_Mapping_Ratio,Sample_Multiple_Mapping_Ratio,Sample_Unmapped_Mapping_Ratio)
write.csv(M2,file="mapping-statistics(rm-replicate).csv")
detach(M1)




