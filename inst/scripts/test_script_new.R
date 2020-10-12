#################### Test script ####################
library(RIPAT)

input_dir = "D:/SMWU/200927/TEST"
input_file = 'A5_15856M_quickmap.txt'
input_file2 = 'I59_15943M_quickmap.txt'
input_file_id = 'A5'
output_dir = input_dir
output_file_name = 'A5_TEST'
output_file_name2 = 'I59_TEST'
random_output_file_name = 'A5_TEST_RANDOM'
random_output_file_name2 = 'I59_TEST_RANDOM'
norandom_output_file_name = "A5_TEST_NORANDOM"

##### 00. Make data files - Do it once at the first time of running RIPAT.
makeData(organism = 'GRCh37')

##### 01. Make R object
### Experimental data
# Single file
blat_obj = makeExpSet(inFile = paste0(input_dir, '/', input_file),
                      mapTool = 'blast', 
                      vectorPos = 'front',
                      outPath = output_dir,
                      outFileName = output_file_name)
blat_obj2 = makeExpSet(inFile = paste0(input_dir, '/', input_file2),
                       mapTool = 'blast', 
                       vectorPos = 'front',
                       outPath = output_dir,
                       outFileName = output_file_name)

# Several files
blat_obj2 = makeInputObj2(inDir = input_dir,
                          id = 'A5',
                          mapTool = 'blat',
                          vectorPos = 'front',
                          outPath = output_dir,
                          outFileName = paste0(output_file_name, '_multiple'))
### Random data
ran_obj = makeRanSet(organism = 'GRCh37',
                     randomSize = 5000,
                     outPath = output_dir,
                     outFileName = 'ran5000')

##### 02. Integration frequency of each chromosome
# Single input
int_on_chr = integPatternByChr(hits = blat_obj, # it can be a list of exp_data
                               ran_hits = ran_obj, # it can be a list of ran_data
                               excelOut = TRUE,
                               isExpList = FALSE, # TRUE if 'hits' is a list
                               isRanList = FALSE, # TRUE if 'ran_hits' is a list
                               outPath = output_dir,
                               outFileName = output_file_name)

##### 03. Gene
blat_gene_random = annoByGene(hits = blat_obj,
                              ran_hits = ran_obj,
                              mapTool = 'blast',
                              organism = 'GRCh37', 
                              interval = 1000,
                              range = c(-5000, 5000),
                              outPath = output_dir,
                              outFileName = random_output_file_name)

blat_gene_random2 = annoByGene(hits = blat_obj2,
                              ran_hits = ran_obj,
                              mapTool = 'blast',
                              organism = 'GRCh37', 
                              interval = 1000,
                              range = c(-5000, 5000),
                              outPath = output_dir,
                              outFileName = random_output_file_name2)


blat_gene_norandom = annoByGene(hits = blat_obj,
                                mapTool = 'blast',
                                organism = 'GRCh37',
                                interval = 1000,
                                range = c(-5000, 5000),
                                outPath = output_dir,
                                outFileName = norandom_output_file_name)

##### 04. CpG
blat_cpg_random = annoByCpG(hits = blat_obj,
                            ran_hits = ran_obj,
                            mapTool = 'blast',
                            organism = 'GRCh37',
                            interval = 1000,
                            range = c(-5000, 5000),
                            outPath = output_dir,
                            outFileName = random_output_file_name)
blat_cpg_random2 = annoByCpG(hits = blat_obj2,
                             ran_hits = ran_obj,
                             mapTool = 'blast',
                             organism = 'GRCh37', 
                             interval = 1000,
                             range = c(-5000, 5000),
                             outPath = output_dir,
                             outFileName = random_output_file_name2)

blat_cpg_norandom = annoByCpG(hits = blat_obj,
                              organism = 'GRCh37',
                              interval = 1000,
                              range = c(-5000, 5000),
                              outPath = output_dir,
                              outFileName = norandom_output_file_name)

##### 05. Pathogenic variant
blat_var_random = annoByVar(hits = blat_obj,
                            ran_hits = ran_obj,
                            organism = 'GRCh37',
                            interval = 1000,
                            range = c(-5000, 5000),
                            outPath = output_dir,
                            outFileName = random_output_file_name)

blat_var_norandom = annoByVar(hits = blat_obj,
                              organism = 'GRCh37',
                              interval = 1000,
                              range = c(-5000, 5000),
                              outPath = output_dir,
                              outFileName = norandom_output_file_name)

##### 06. Repeat
blat_repeat_random = annoByRepeat(hits = blat_obj,
                                  ran_hits = ran_obj,
                                  organism = 'GRCh37',
                                  interval = 1000,
                                  range = c(-5000, 5000),
                                  outPath = output_dir,
                                  outFileName = random_output_file_name)

blat_repeat_norandom = annoByRepeat(hits = blat_obj,
                                    organism = 'GRCh37',
                                    interval = 1000,
                                    range = c(-5000, 5000),
                                    outPath = output_dir,
                                    outFileName = norandom_output_file_name)
##### 07. Cancer gene
blast_cgene_random = annoByCGC(hits = blat_obj,
                               organism = 'GRCh37',
                               interval = 1000,
                               range = c(-5000, 5000),
                               outPath = output_dir,
                               outFileName = random_output_file_name)
blast_cgene_norandom = annoByCGC(hits = blat_obj,
                                 organism = 'GRCh37',
                                 interval = 1000,
                                 range = c(-5000, 5000),
                                 outPath = output_dir,
                                 outFileName = norandom_output_file_name)

##### 08. Karyogram
# Single file
karyo_gene = drawingKaryo(hits = blat_obj,
                          organism = "GRCh37",
                          feature = blat_gene_norandom$Gene_data,  
                          includeUndecided = FALSE,
                          outPath = output_dir,
                          outFileName = paste0(output_file_name, '_single'))

# Several files
karyo_gene2 = drawingKaryo(hits = blat_obj,
                           organism = "GRCh37",
                           feature = blat_gene_norandom$Gene_data,  
                           includeUndecided = FALSE,
                           outPath = output_dir,
                           outFileName = paste0(output_file_name, '_multiple'))

##### 09. Make result files
results = makeDocument(res = blat_gene_random,
                       resType = 'gene',
                       interval = 5000,
                       range = c(-20000, 20000),
                       includeUndecided = FALSE,
                       outPath = output_dir,
                       outFileName = paste0(output_file_name, '_multiple'))


