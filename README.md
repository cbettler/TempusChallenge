# Tempus Challenge

For this challenge, I was tasked with creating a prototype variant annotation tool for a VCF file. All code for this challenge was written in R. The input file was supplied by Tempus. To run this code, simply download it and open it in R, download the VCF file, and then change the file path name. The resulting output is a .csv file. 

As a reference, here are the exact instructions: 

For this challenge, you are asked to prototype a variant annotation tool. We will provide you with a VCF file, and you will create a small software program to output a table annotating each variant in the file. Each variant must be annotated with the following pieces of information:
1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, annotate with the most deleterious possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API
(API documentation is available here: http://exac.hms.harvard.edu/)
6. Additional optional information from ExAC that you feel might be relevant.




