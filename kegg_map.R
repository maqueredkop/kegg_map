
###############################################################################
# This program is used to download the KEGG database by the 'KEGGREST' package.
# The main goal is to get the pathway infomation of every gene.
# The input is a kegg hsa id list file of all human gene, eg 'hsa:222029'
###############################################################################

print(Sys.time())

library(KEGGREST)
library(stringr)

# get the whole human gene hsa_id
hsa_id_list <- read.table('kegg_hsa_id.txt')

##whole_human_gene <- c('hsa:100287010', 'hsa:100506548', 'hsa:100288846', 'hsa:222029')

success_file = 'kegg_map.txt'
fail_file = 'kegg_map_fail.txt'

cat('hsa_id', '|', 'gene_symbol', '|', 'HGNC_id', '|', 'uniprot_id', '|', 'defination', '|', 'pathway', file = success_file, append = TRUE)
cat('\n', file = success_file, append = TRUE)

for (i in 1:length(hsa_id_list$V1)){
    # remove the attributes of the i 
    hsa_id <- as.character(hsa_id_list$V1[i])
    # keggGet...    
    ##query <- keggGet(hsa_id)
    query <- tryCatch(keggGet(hsa_id), error = function(e) return('error'))
    
    if (query == 'error'){
        cat(hsa_id, fill = TRUE, file = fail_file, append = TRUE)
        next
    }
    
    # next we extract the needed info from the query
    gene_symbol <- query[[1]]$NAME
    if (length(gene_symbol) == 0){
        gene_symbol <- ' '
    }
    
    HGNC_id <- str_subset(query[[1]]$DBLINKS, 'HGNC:')
    if (length(HGNC_id) > 0){
        
        HGNC_id <- str_split(HGNC_id, ': ')[[1]][2]
    }else{
        HGNC_id <- ' '
    }
    
    uniprot_id <- str_subset(query[[1]]$DBLINKS, 'UniProt:')
    
    if (length(uniprot_id) > 0){
      
        uniprot_id <- str_split(uniprot_id, ': ')[[1]][2]
    }else{
        uniprot_id <- ' '
    }
    
    defination <- query[[1]]$DEFINITION
    if (length(defination) == 0){
        defination <- ' '
    }
    
    pathway <- as.character(query[[1]]$PATHWAY)
    if (length(pathway) == 0){
        pathway <- ' '
    }else{
        pathway <- str_c(pathway, collapse = '//')
    }
    
    cat(hsa_id, '|', gene_symbol, '|', HGNC_id, '|', uniprot_id, '|', defination, '|', pathway, file = success_file, append = TRUE)
    cat('\n', file = success_file, append = TRUE)

    Sys.sleep(0.5)    
}

print('Done...')

print(Sys.time())
