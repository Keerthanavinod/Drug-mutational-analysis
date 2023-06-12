setwd("~/Downloads")
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)

#load gdsc data
gdsc<- fread("Cell_listFri Dec 16 05_30_03 2022.csv")

#filter version 2 data
gdsc<- filter(gdsc, gdsc$Datasets=="GDSC2")

#filter AML data
aml_samples<- filter(gdsc, gdsc$`Tissue sub-type`=="acute_myeloid_leukaemia")

ic50_cell_line<- fread("LAML_IC_Fri Dec 16 11_53_49 2022.csv")


#cosmic data
complete_mutation_data<- fread("CosmicCLP_MutantExport (1).tsv/CosmicCLP_MutantExport.tsv")

#filter for aml
filtered__aml<-filter(complete_mutation_data,complete_mutation_data$`Histology subtype 1`=="acute_myeloid_leukaemia")

#rearrage
final_aml_mutation_data<- filtered__aml[,-c(9:12,14,15)]

#selecting aml cell lines
cell_line<- aml_samples$`Cell line Name`

#filtering out only those cell lines
final_aml_mutation_data<- filter(final_aml_mutation_data, final_aml_mutation_data$`Sample name`%in% cell_line)

#removing unwanted mutations
data<- final_aml_mutation_data[final_aml_mutation_data$`Mutation Description`!="Unknown" & final_aml_mutation_data$`Mutation Description`!="Substitution - coding silent" & final_aml_mutation_data$`Mutation Description`!="Nonstop extension",]

#editing gene names
x<-data.frame(data$`Gene name`)
x<-str_split_fixed(x$data..Gene.name.,"_",2)
c_names<-x[1:20631]
data$`Gene name`<-c_names

#creating binary matrix
matrix<-data %>%
  pivot_wider(`Gene name`,
              names_from = `Sample name`,
              values_from =`Sample name`,
              values_fn = list(`Sample name` = length),
              values_fill = list(`Sample name`= 0))

final_data<-data%>%
  distinct(`Gene name`, `Sample name`, .keep_all = TRUE)

final_matrix<-final_data %>%
  pivot_wider(`Gene name`,
              names_from = `Sample name`,
              values_from =`Sample name`,
              values_fn = list(`Sample name` = length),
              values_fill = list(`Sample name`= 0))

#creating a new column total mutations and suming up the mutations row-wise
final_matrix$Total_mutations<-rowSums(final_matrix[,-1])

#arranging in decreasing order
final_matrix<-final_matrix[order(final_matrix$Total_mutations, decreasing = TRUE),]

#selecting a mutation frequency cutoff of min 5
cut_off_5<- filter(final_matrix, final_matrix$Total_mutations>=5)

#removing unwanted genes
remove_unwanted_genes<- cut_off_5[-c(1,2,7,15),]
binary_matrix<- filter(final_matrix, final_matrix$Total_mutations>=5)
binary_matrix<- binary_matrix[-c(1,2,7,15,29),]
gene_name<- binary_matrix$`Gene name`

#if mutated 1 else 0
binary_matrix_final<- as.data.frame(ifelse(binary_matrix[,]==1,'Mut','Wt'))
binary_matrix_final$`Gene name`<-gene_name

binary_matrix_final<- binary_matrix_final[,-25]
binary_matrix_final<- gather(binary_matrix_final, 'Cell Line Name', 'Type',2:24 )
merged<- merge(binary_matrix_final,ic50_cell_line,by="Cell Line Name")
merged<- merged[,c(1:6,10)]
write.csv(merged,"/home/keerthana2/AML/GDSC/boxplot.csv",row.names = F)

##loop for generating boxplots for all gene-drug combinations
tp53<- filter(merged, merged$`Gene name`=="TP53")
unique(tp53$`Drug Name`)
drug_names<- as.character(unique(tp53$`Drug Name`))
length(drug_names)
for(i in 1:length(drug_names)){
  ggplot(tp53, aes(x=`Gene name`, y=IC50, fill=Type)) +
    geom_boxplot()+
    xlab("TP53")+
    ylab("IC50 OF 'i'")+
    ggtitle("Drug sensitivity of cell lines")
  
}




