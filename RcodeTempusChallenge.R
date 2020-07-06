library(stringr)
library(sjmisc)
library(httr)

vcf_header <- c()
vcf_main_info <- c()
vcf_names <- c()
name_List <- c()
nameList_edited <- c()
depth <- c()
chrom <- c()
pos <- c()
ref <- c()
alt <- c()
num_var_reads <- c()
percent_reads <- c()
exac_input <- c()
exac_query <- c()
query_info_extracted <- c()
exac_allele_freq <- c()
exac_ensembl_info <- c()
ensembl_transposed <- c()


input_file_path <- "/Users/carleebettler/Downloads/Challenge_data.vcf"
#CHANGE THIS if running this code on a different computer 

vcf_header <-readLines(input_file_path)
vcf_main_info <-read.table(input_file_path , stringsAsFactors = FALSE)
#Reading in the vcf header information and main body information 


vcf_header  <- vcf_header[-(grep("#CHROM",vcf_header)+1):-(length(vcf_header))]
vcf_names<-unlist(strsplit(vcf_header[length(vcf_header)],"\t"))
names(vcf_main_info)<-vcf_names
vcf_main_info$INFO <- as.character(vcf_main_info$INFO)
#Parsing vcf file and extracting file column names 

for(i in 1:nrow(vcf_main_info)){
  #This loop extracts specific data from the vcf file needed for the annotation info
  chrom[i] <- vcf_main_info$`#CHROM`[[i]]
  pos[i] <- vcf_main_info$POS[[i]]
  ref[i] <- vcf_main_info$REF[[i]]
  
  alt[i] <- strsplit((vcf_main_info$ALT), ",")[[i]]
  #For simplicity's sake, used the first alt nucleotide sequence listed 
  
  variantType <- strsplit((vcf_main_info$INFO), ";")[[i]]
  
  name_List[i] <- str_sub(variantType[41], 6)
  #Pulls the type(s) of variation(s)
  
  depth[i] <- str_sub(variantType[8], 4)
  #DP used for depth of sequence coverage at the site of variation.
  
  num_var_reads[i] <- str_sub(variantType[3], 4)
  #AC used for the number of reads supporting the variant. Not 100% sure if this is right but I wasn't sure what else to use.
  
  percent_reads[i] <- as.numeric(str_sub(variantType[4], 4))*100
  #AF used for percentage of reads supporting the variant. Not 100% sure if this is right but I wasn't sure what else to use.
  
  exac_input[i] <- paste0(chrom[i],"-",pos[i],"-", ref[i],"-",alt[i])
  #Uses this exac input format: A JSON array of variant strings(CHROMOSOME-POSITION-REFERENCE-VARIANT).
}

exac_query <- httr::POST(url="http://exac.hms.harvard.edu/rest/bulk/variant", body=jsonlite::toJSON(as.character(exac_input, encode = "json")))

query_info_extracted <- content(exac_query)

for (i in seq_along(query_info_extracted)) {   
  if(is.null(query_info_extracted[[i]]$variant$allele_freq)){
    exac_allele_freq[i] <- "No data available"
  }else{
    exac_allele_freq[i] <- query_info_extracted[[i]]$variant$allele_freq
  }
}
#Loop pulls exac data for the allele frequency of variant from Broad Institute ExAC Project API

for (i in seq_along(query_info_extracted)) {   
  if(is.null(query_info_extracted[[i]]$variant$genes)){
    exac_ensembl_info[i] <- "No data available"
  }else{
    exac_ensembl_info[i] <- query_info_extracted[[i]]$variant$genes[i]
  }
}

ensembl_transposed <- t(t(exac_ensembl_info))

#Loop pulls exac data for the gene Ensembl ID. I thought this would be a useful metric to get more information. Also, the formatting for this variable was off, so I had to transpose this data. 

for (i in 1:length(name_List)){
  #This for loop annotates the type of variation, choosing which is the most deleterious 
  myStr <- name_List[i]
  str <- c("complex","ins")
  
  if(all(sapply(c("complex"), grepl, myStr))){
    nameList_edited[i] <- "Complex Mutation"
    #If the type of variation contains complex or multiple options and complex, complex has the potential to be the most deleleterious (not always, but potentially) so that becomes the default 
  }
  
  else if(all(sapply(c("ins", "del"), grepl, myStr))){
    nameList_edited[i] <- "Indel"
    #If the type of variation contains and insertion and deletion, regardless of what else is there, these are the most deleterious variations so indel becomes the default
  }
  else if(all(sapply(c("ins"), grepl, myStr))){
    nameList_edited[i] <- "Insertion"
    #If the type of variation contains an insertion, or an insertion and a snp or mnp, the insertion is the most deleterious variation so insertion becomes the default
  }
  else if(all(sapply(c("del"), grepl, myStr))){
    nameList_edited[i] <- "Deletion"
    #If the type of variation contains a deletion, or a deletion and a snp or mnp, the deletion is the most deleterious variation so deletion becomes the default
  }
  else if(all(sapply(c("mnp"), grepl, myStr))){
    nameList_edited[i] <- "Substitution: Multi-Nucleotide Polymorphisms"
    #If the type of variation contains a mnp, or an mnp and a snp, the mnp is the most deleterious variation so mnp becomes the default
  }
  else{
    nameList_edited[i] <- "Substitution: Single-Nucleotide Polymorphism"
    #Only snp(s) are left, everything else has been assigned 
  }
}

output_file_path <- "/Users/carleebettler/annotated_VCF_file.csv"
#CHANGE THIS if running this code on a different computer   


df <- data.frame(Chromosome = chrom, Position = pos, Type_of_Variation = nameList_edited, Depth_of_Sequence_Coverage = depth, Number_of_Reads_Supporting_Variant = gsub(",.*$", "", num_var_reads), Percentage_of_Reads_Supporting_Variant = gsub(",.*$", "",percent_reads), Allele_Frequency_from_ExAC = exac_allele_freq, Ensembl_Gene_ID=as.character(ensembl_transposed))
#Relevant annotation is converted to a data frame so it can be written to a csv file, a file type I chose becauses it is generally easy to work with and parse if further analysis is desired. 


write.csv(df,output_file_path, row.names = FALSE)