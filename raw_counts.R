setwd ("/home/das260/star_daokun_ouput") # set as the output_path

#raw counts files
RC1=read.table("raw_counts_gene_id.txt",header = F)
RC2=read.table("raw_counts_gene_names.txt",header = F)
RC3=read.table("raw_counts_names.txt",header = F)
filelist = list.files(pattern = "*raw_counts.txt")
filelist1 = list.files(pattern = "*raw_counts_summary.txt")
RC5=cbind(RC1,RC2)
RC5=RC5[1:(dim(RC1)[1]-5),]
#or we can use data=lapply(filelist, FUN=read.table, header=TRUE) datafr = do.call("cbind", data)
sample_file_name=c()
for (i in filelist)
{
  tem.data=read.table(i,header=F)
  RC5=cbind(RC5,tem.data)
}
for (j in 1:dim(RC3)[1])
{
  sample_file_name[j]=toString(RC3[j,1])
}
colnames(RC5)=sample_file_name
write.csv(RC5, file="raw-counts.csv")

data=lapply(filelist1, FUN=read.table, header=F) 
datafr = do.call("cbind", data)
RC6=datafr[(dim(datafr)[1]-4):dim(datafr)[1],]
colnames(RC6)=sample_file_name[3:length(sample_file_name)]
rownames(RC6)=c("No_Feature","Ambiguous","Too_Low_Quality","Not_Aligned","Alignment_Not_Unique")
write.csv(RC6, file="raw-counts-summary.csv")
