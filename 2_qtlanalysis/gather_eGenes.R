setwd('/rsrch5/scratch/epi/abhattacharya3/GTEx_compare')
load('/rsrch5/scratch/epi/abhattacharya3/GTEx_compare/ensembl_gene_info.RData')

setwd('/rsrch5/scratch/epi/bhattacharya_lab')
tissues = list.files()

all_egenes = all_egenes_common = data.frame(Tissue = character(),
                        Version = character(),
                        Method = character(),
                        ensembl_gene_id = character(),
                        Total_eGenes = numeric(),
                        Total = numeric())
for (t in tissues){
  
  print(t)
  setwd(file.path('/rsrch5/scratch/epi/bhattacharya_lab/',
                  t))
  
  # List all files in current directory and filter for cisQTL.txt files
  all_files <- list.files(pattern = "\\.txt$")
  
  # Filter for files ending in cisQTL.txt (not cisQTL_nominal.txt)
  cisqtl_files <- all_files[grepl("cisQTL\\.txt$", all_files) & !grepl("nominal", all_files)]
  
  
  # Extract version and method information
  library(stringr)
  
  file_info <- data.frame(
    filename = cisqtl_files,
    version = str_extract(cisqtl_files, "(?<=v)\\d+"),
    method = str_extract(cisqtl_files, "(star|salmon)"),
    stringsAsFactors = FALSE
  )
  
  for (f in 1:nrow(file_info)){
    
    dat = data.table::fread(file_info$filename[f])
    if (!grepl('common',file_info$filename[f])){
      
      all_egenes = rbind(all_egenes,
                         data.frame(Tissue = t,
                                    Version = file_info$version[f],
                                    Method = file_info$method[f],
                                    ensembl_gene_id = dat$V1[dat$V20 < .05],
                                    Total_eGenes = length(dat$V1[dat$V20 < .05]),
                                    Total = nrow(dat)))
      
    } else {
      all_egenes_common = rbind(all_egenes_common,
                                data.frame(Tissue = t,
                                           Version = file_info$version[f],
                                           Method = file_info$method[f],
                                           ensembl_gene_id = dat$V1[dat$V20 < .05],
                                           Total_eGenes = length(dat$V1[dat$V20 < .05]),
                                           Total = nrow(dat)))
    }
    
  }
  
  
  
}

all_egene_num = all_egenes[,c('Tissue','Version','Method',
                              'Total_eGenes','Total')]
all_egene_num = all_egene_num[!duplicated(all_egene_num),]
all_egene_common_num = all_egenes_common[,c('Tissue','Version','Method',
                                            'Total_eGenes','Total')]
all_egene_common_num = all_egene_common_num[!duplicated(all_egene_common_num),]

gene_info_this = merge(all_egenes[,c('Tissue','ensembl_gene_id')],
                       gene_info,by='ensembl_gene_id')
gene_info_this = gene_info_this[!duplicated(gene_info_this),]

saveRDS(gene_info_this,
        '/rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis/eGene_colocalization.RDS')
saveRDS(list(total = all_egenes,
             common = all_egenes_common),
        '/rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis/eGene_intersection_forplots.RDS')