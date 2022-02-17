options(stringsAsFactors = FALSE)

installed_pkgs <- installed.packages()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  # The following initializes usage of Bioc devel
if (length(setdiff("fgsea", installed_pkgs)) > 0) {
    BiocManager::install("fgsea")
}
if (length(setdiff("plyr", installed_pkgs)) > 0) {
    install.packages("plyr")
}

if (length(setdiff("dplyr", installed_pkgs)) > 0) {
    install.packages("dplyr")
}
if (length(setdiff("data.table", installed_pkgs)) > 0) {
    install.packages("data.table")
}

if (length(setdiff("boot", installed_pkgs)) > 0) {
    install.packages("boot")
}
if (length(setdiff("snow", installed_pkgs)) > 0) {
    install.packages("snow")
}
if (length(setdiff("parallel", installed_pkgs)) > 0) {
    install.packages("parallel")
}
if (length(setdiff("stringr", installed_pkgs)) > 0) {
    install.packages("stringr")
}
if (length(setdiff("optparse", installed_pkgs)) > 0) {
    install.packages("optparse")
}
if (length(setdiff("collections", installed_pkgs)) > 0) {
    install.packages("collections")
}

suppressMessages(library(fgsea))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(snow))
suppressMessages(library(parallel))
suppressMessages(library(stringr))
suppressMessages(library(optparse))
suppressMessages(library(collections))

option_list = list(
    make_option(c("--folder"), action = "store", default = "ex", type = 'character', help = "Project folder to be used for storing data used and output by this program."),
    make_option(c("--name"), action = "store", default = "ex_", type = 'character', help = "Task name prefix."),
    make_option(c("--stage1"), action = "store", default = TRUE, type = 'logical', help = "Whether perform CSEA analysis."),
    make_option(c("--stage2"), action = "store", default = TRUE, type = 'logical', help = "Whether perform mediation analysis."),
    make_option(c("--ga"), action = "store", default = paste0("./ex/mutation_processed_with_annotation.tsv"), type='character', help="Tabular file of genetic alteration data (i.e. mutations)."),
    make_option(c("--exp"), action = "store", default = "./data/expr_processed.tsv", type = 'character', help="Tabular file of gene expression data."),
    make_option(c("--dep"), action = "store", default = "./data/dependency_processed.tsv", type = 'character', help="Tabular file of gene dependency data [Optional]"),
    make_option(c("--cna"), action = "store", default = "./data/CNV_processed.tsv", type = 'character', help="Tabular file of CNV data [Optional]"),
    make_option(c("--cl"), action = "store", default = FALSE, type = 'logical', help = "Whether the data is from DepMap project. [Optional]"),
    make_option(c("--map"), action = "store", default = "./ex_TCGA/TCGA_sample_mapping.tsv", type = 'character', help = "Tabular file of gene sample type mapper: map user-provided to in-house tissue types [Optional]"),
    make_option(c("--ppi"), action = "store", default = "./data/PPI_filtered.txt", type = 'character', help = "Tabular file of two columns of genes, i.e. two-column tab delimited undirected PPI interactions. [Optional]"),
    make_option(c("--start_PPI"), action = "store", default = -1, type = 'integer', help = "Starting line of --ppi file. Set to -1 if no subsetting is needed [Optional]"),
    make_option(c("--end_PPI"), action = "store", default = -1, type = 'integer', help = "Ending line of --ppi file. Set to -1 if no subsetting is needed [Optional]"),
    make_option(c("--directed"), action = "store", default = "FALSE", type = 'logical', help = "If TRUE, only check the regulation of the gene in the first column upon the gene in the second column. [Optional]"),
    make_option(c("--reversedirected"), action="store", default="FALSE", type = 'logical', help = "If TRUE, only check the regulation of the gene in the second column upon the gene in the first column. [Optional]"),
    make_option(c("--sample"), action = "store", default = -1, type = 'integer', help = "Random sampling how many lines from the provided gene file (indicated with --genes). Set to -1 if no sampling is needed [Optional]"),
    make_option(c("--genes"), action = "store", default = "", type = 'character', help="Genes to be subsetted, please also set --sample to a positive integer. [Optional]"),
    make_option(c("--core"), action = "store", default = 1, type = 'integer', help = "Cpu cores used to run this script concurrently. [Optional]"),
    make_option(c("--seed"), action = "store", default = 1, type = 'integer', help = "Seed used in set.seed(). Set to -1 to use rownumbers of PPI file as seeds. [Optional]"),
    make_option(c("--npermut"), action = "store", default = 100, type = 'integer', help = "Number of permutations for fgsea. [Optional]"),
    make_option(c("--nproc"), action = "store", default = 1, type = 'integer', help = "Number of processes to be used for fgsea. [Optional]"),
    make_option(c("--min_size"), action = "store", default = 2, type = 'integer', help = "Min cell line number used in fgsea. [Optional]"),
    make_option(c("--max_size"), action = "store", default = 500, type = 'integer', help = "Max cell line number used in fgsea. [Optional]"),
    make_option(c("--precise_mediation_pvalue"), action = "store", default = FALSE, type = 'logical', help="Whether more precise p-values are needed (TRUE: 100000 permutations; FALSE: 1000 permutations). [Optional]"),
    make_option(c("--remove_cis_confounding"), action = "store", default = TRUE, type = 'logical', help = "Whether to remove trans records with potential cis regulation. [Optional]"),
    make_option(c("--cna_adjust"), action = "store", default = FALSE, type = 'logical', help = "Whether remove CNA confounding effect. Specify --cna also when using this argument [Optional]"),
    make_option(c("--cna_discretization"), action="store", default = FALSE, type = 'logical', help = "Whether discretize the CNA values. Specify --cna and --cna_adjust also when using this argument [Optional]"),
    make_option(c("--cna_permutation_size"), action = "store", default = 50, type = 'integer', help = "Calculate empirical p-values based on the predefined permutation size. Specify --cna and --cna_adjust also when using this argument [Optional]"),
    make_option(c("--cna_adjust_Pvalue_cutoff"), action = "store", default = 0.1, type = 'double', help = "Remove CNA confounding effect according to the P-value cutoff. Specify --cna and --cna_adjust also when using this argument [Optional]"),
    make_option(c("--cna_amplification_cutoff"), action = "store", default = 1, type = 'double', help = "CN > this cutoff is considered as amplication. Specify --cna and --cna_adjust also when using this argument [Optional]"),
    make_option(c("--cna_deletion_cutoff"), action = "store", default = 0.5, type = 'double', help = "CN < this cutoff is considered as deletion. Specify --cna and --cna_adjust also when using this argument [Optional]"),
    make_option(c("--remove_contradiction"), action = "store", default = TRUE, type = 'logical', help = "Whether remove contradiction between CSEA and mediation analysis. [Optional]"),
    make_option(c("--FDR_cutoff_cis"), action = "store", default = 0.1, type = 'double', help = "The FDR cutoff for defining cis genetic alterations. [Optional]"),
    make_option(c("--FDR_cutoff_trans"), action = "store", default = 0.01, type = 'double', help = "The FDR cutoff for defining trans genetic alterations. [Optional]"),
    make_option(c("--CSEA_cis_cutoff"), action = "store", default = 0.25, type = 'double', help = "The CSEA p-value cutoff for removing insignificant genetic alterations before p-value combination. [Optional]"),
		make_option(c("--mediation_cis_cutoff"), action = "store", default = 0.25, type = 'double', help = "The mediation p-value cutoff for removing insignificant genetic alterations before p-value combination. [Optional]"),
		make_option(c("--CSEA_trans_cutoff"), action = "store", default = 0.05, type = 'double', help = "The CSEA p-value cutoff for removing insignificant genetic alterations before p-value combination. [Optional]"),
		make_option(c("--mediation_trans_cutoff"), action = "store", default = 0.05, type = 'double', help = "The mediation p-value cutoff for removing insignificant genetic alterations before p-value combination. [Optional]"),
    make_option(c("--output_mutation_site_centric_table"), action = "store", default = TRUE, type = 'logical', help="Output mutation information tables [Optional]")
)

arguments = parse_args(OptionParser(option_list = option_list))

if(arguments$start_PPI > 0 & arguments$end_PPI > 0){
  arguments$name <- paste0(arguments$name, arguments$start_PPI, "_", arguments$end_PPI, "_")
}


# Print Parameters
cat("PARAMETERS:\n", file = stdout())
for(arg.idx in 1 : length(arguments)) {
  cat(paste0(names(arguments)[arg.idx], "\t", arguments[[arg.idx]], "\n"), file = stdout())
}
cat("------------------\n", file = stdout())

output_dir <- file.path(getwd(), arguments$folder)
if (!dir.exists(output_dir)){
    dir.create(output_dir)
}
cat("Using ",output_dir," as the working directory\n", file = stdout())


######### Part 1: Load omics data #########

cat("1 of 7 Reading data\n", file = stdout())
# Map disease samples to DepMap cell_lines
mapping <- NULL
if(!arguments$cl){
    if(file.exists(arguments$map)){
        mapping <- fread(arguments$map, sep="\t", header = TRUE)
    }else{
        stop("No mapping file for user-provided data.")
    }
}


gene_subsets <- NULL

if(file.exists(arguments$genes)){
    gene_subsets <- read.table(arguments$genes, sep="\t", header = F)
    
    if(arguments$sample > 0){ # sampling
        set.seed(arguments$seed)
        gene_subsets <- as.character(gene_subsets[sample(nrow(gene_subsets), arguments$sample), 1])
        cat("Read", paste0(arguments$genes, " with ", length(gene_subsets), " rows\n"), file = stdout())
    }
    
}

Input.net <- NULL
mutations_dt <- NULL
dependency_dt <- NULL
CNA_dt <- NULL
exp_tpm_dt <- NULL

if(file.exists(arguments$ppi)){
  Input.net <- read.table(arguments$ppi, header = TRUE, sep = '\t')
  if(arguments$start_PPI > 0 & arguments$end_PPI > 0 & arguments$end_PPI > arguments$start_PPI){
    Input.net <- Input.net[arguments$start_PPI : arguments$end_PPI,]
  }else if(arguments$start_PPI < 0 & arguments$end_PPI < 0){
    
  }else{
    stop("End_PPI should be larger than start_PPI and both should be larger than zero.")
  }
  PPI_genes <- unique(c(Input.net$V1, Input.net$V2))
  if(!is.null(gene_subsets)){
    Input.net <- Input.net[Input.net$V1 %in% gene_subsets | Input.net$V2 %in% gene_subsets,]
  }
  
}else{
  stop(paste0("File path ", arguments$ppi, " is not correctly specified."))
}

cat("Read",paste0(arguments$ppi, " with ", nrow(Input.net), " rows after filtering.\n"), file = stdout())

if(file.exists(arguments$ga)){
    mutations_dt <- fread(arguments$ga, header =TRUE, sep ='\t')
    if(!("cell_line" %in% colnames(mutations_dt)) & !("-" %in% colnames(mutations_dt))){
      stop(paste0("Data ", arguments$ga, " does not contain columns of cell_line or disease."))
    }
    
    setkey(mutations_dt, gene_name)
    if(nrow(Input.net) > 0){
        mutations_dt <- mutations_dt[.(unique(c(Input.net$V1, Input.net$V2)))]
    }
    setkey(mutations_dt, gene_name)
}else{
  stop(paste0("File path ", arguments$ga, " is not correctly specified."))
}

cat("Read",paste0(arguments$ga, " with ", nrow(mutations_dt), " rows after filtering\n"), file = stdout())

if(file.exists(arguments$dep)){
    dependency_dt <- fread(arguments$dep, header =TRUE, sep ='\t')
    setkey(dependency_dt, gene)
    if(nrow(Input.net) > 0){
      dependency_dt <- dependency_dt[.(unique(c(Input.net$V1, Input.net$V2)))]
    }
    setkey(dependency_dt, gene)
}else{
    stop(paste0("File path ", arguments$dep, " is not correctly specified."))
}

cat("Read", paste0(arguments$dep, " with ", nrow(dependency_dt), " rows after filtering.\n"), file = stdout())




if(file.exists(arguments$exp)){
    exp_tpm_dt <- fread(arguments$exp, header =TRUE, sep ='\t')
    if(!("cell_line" %in% colnames(exp_tpm_dt)) & !("-" %in% colnames(exp_tpm_dt))){
        stop(paste0("Data ", arguments$exp, " does not contain columns of cell_line or disease."))
    }
    setkey(exp_tpm_dt, gene)
    if(nrow(Input.net) > 0){
        exp_tpm_dt <- exp_tpm_dt[.(unique(c(Input.net$V1, Input.net$V2)))]
    }
    setkey(exp_tpm_dt, gene)
}else{
    stop(paste0("File path ", arguments$exp, " is not correctly specified."))
}
cat("Read",paste0(arguments$exp, " with ", nrow(exp_tpm_dt), " rows after filtering.\n"), file = stdout())

######### Part 2: CSEA and mediation analysis #########
# 2.1 CSEA

mut_types <- c(unique(mutations_dt$var_class), "Any")
mut_types <- unique(mut_types[mut_types != ""])
mut_types <- sort(mut_types)

npermut <- arguments$npermut # Number of permutations for fgsea.
nproc <- arguments$nproc # Number of processes to be used for fgsea.
min.size <- arguments$min_size
max.size <- arguments$max_size


unique_genes <- unique(c(Input.net[,1], Input.net[,2]))

CSEA <-
  function(index = NA, Gene_A = "", Gene_B = "", npermut = 100, nproc = 10, min.size = 1, max.size = 100, seed = 1, include_cis_analysis = TRUE){
    # Gene_A: mutations of multiple genes
    # Gene_B: expression of Gene_B
    tmp_data_frame <- data.frame("index" = vector("numeric", length = 0), "Gene_A" = vector("character", length = 0), "Gene_B" = vector("character", length = 0), "Mutation_type" = vector("character", length = 0), "Regulation" = vector("character", length = 0), "ES" = vector("numeric",length = 0), "NES" = vector("numeric", length = 0), "P_val" = vector("numeric",length = 0), "Leading_edge" = vector("character", length = 0))
    
    if(is.null(key(exp_tpm_dt))){
      setkey(exp_tpm_dt, gene)
    }
    Gene_A <- unlist(strsplit(Gene_A, ","))
    exp_dat_gene_2 <- exp_tpm_dt[.(Gene_B),]
    exp_dat_gene_2 <- exp_dat_gene_2[!duplicated(exp_dat_gene_2),]
    if(nrow(exp_dat_gene_2) < 1){
      return(tmp_data_frame)
    }

    exp_dat_gene_2 <- exp_dat_gene_2[!is.na(exp_dat_gene_2$cell_line),]
    RNK2 <- exp_dat_gene_2$rna_expression
    names(RNK2) <- exp_dat_gene_2$cell_line
    RNK2 <- sort(RNK2, decreasing = FALSE)
    if(is.null(key(mutations_dt))){
      setkey(mutations_dt, gene_name)
    }
    unique_mut_genes <- NULL
    if(include_cis_analysis){
      unique_mut_genes <- unique(c(Gene_A, Gene_B))
    }else{
      unique_mut_genes <- Gene_A
    }
    mutations_genes <- mutations_dt[.(unique_mut_genes), c("gene_name", "var_class", "cell_line")]
    ctx_mut_genes <- list()
    for(uniq in unique_mut_genes){
      mutations_gene1 <- mutations_genes[which(mutations_genes$gene_name %in% uniq),]
      for(i in 1 : length(mut_types)){
        if(mut_types[i] == "Any"){
            ctx_mut_genes[[paste0(uniq, "::Any")]] <- mutations_gene1$cell_line
        }else{
            tmp <- mutations_gene1[mutations_gene1$var_class %in% mut_types[i],]
            ctx_mut_genes[[paste0(uniq, "::", mut_types[i])]] <- tmp$cell_line
        }

      }
    }

    # Perform the fGSEA analysis along all the cell lines
    if(length(RNK2) > 0){
      set.seed(seed)
      res_pos <- fgseaMultilevel(pathways = ctx_mut_genes, stats = RNK2, minSize = min.size, maxSize = max.size, sampleSize = 500, nPermSimple = npermut, nproc = nproc, scoreType = "pos")
      set.seed(seed)
      res_neg <- fgseaMultilevel(pathways = ctx_mut_genes, stats = RNK2, minSize = min.size, maxSize = max.size, sampleSize = 500, nPermSimple = npermut, nproc = nproc, scoreType = "neg")
      
      gene_mut <- sapply(res_pos$pathway, function(z) strsplit(z, split="::")[[1]][1]) # ARRB2::Any -> ARRB2
      pathway1 <- sapply(res_pos$pathway, function(z) strsplit(z, split="::")[[1]][2]) # ARRB2::Any -> Any

      if(nrow(res_pos) > 1){
        for(j in 1 : nrow(res_pos)){
          p1 <- 1
          if(!is.na(res_pos$pval[j])){
            p1 <- res_pos$pval[j]
          }
          p2 <- 1
          if(!is.na(res_neg$pval[j])){
            p2 <- res_neg$pval[j]
          }
          min_index <- which.min(c(p1, p2))
          reg <- "Neutral"
          if(min_index == 1){
            if(res_pos$pval[j] < 0.01){
              reg <- "Up"
            }
            tmp_data_frame <- rbind(tmp_data_frame, data.frame("index" = index, "Gene_A" = gene_mut[j], "Gene_B" = Gene_B, "Mutation_type" = pathway1[j], "Regulation" = reg, "ES" = res_pos$ES[j] , "NES" = res_pos$NES[j], "P_val" = res_pos$pval[j], "Leading_edge" = paste(res_pos$leadingEdge[[j]], collapse = ',')))
          }else{
            if(res_neg$pval[j] < 0.01){
              reg <- "Down"
            } 
            tmp_data_frame <- rbind(tmp_data_frame, data.frame("index" = index, "Gene_A" = gene_mut[j], "Gene_B" = Gene_B, "Mutation_type" = pathway1[j], "Regulation" = reg, "ES" = res_neg$ES[j], "NES" = res_neg$NES[j], "P_val" = res_neg$pval[j], "Leading_edge" = paste(res_neg$leadingEdge[[j]], collapse = ',')))
          }
        }
      }
    }

    return(tmp_data_frame)
  }


params <- NULL
gene_dict <- dict()
include_cis_analysis <- NULL
row_ppi <- nrow(Input.net)
seed_set <- arguments$seed

if(arguments$stage1){
    cat("Step 2 of 7, Constructing parameter tables for CSEA analysis.\n",file = stdout())
    
    for(index in 1:row_ppi){
        
        if(arguments$seed < 0){
            seed_set = index
        }
        
        if(arguments$directed){
          if(gene_dict$has(Input.net[index, 2])){
              include_cis_analysis <- FALSE
          }else{
            include_cis_analysis <- TRUE
            gene_dict$set(Input.net[index, 2], "");
          }
          params <- rbind(params,data.frame("index" = index, "Gene_A" = Input.net[index, 1], "Gene_B" = Input.net[index, 2], "npermut" = npermut, "nproc" = nproc, "min.size" = min.size, "max.size" = max.size, "seed" = seed_set, "include_cis_analysis" = include_cis_analysis))
          
        }
        if(arguments$reversedirected){
          if(gene_dict$has(Input.net[index, 1])){
              include_cis_analysis <- FALSE
          }else{
            include_cis_analysis <- TRUE
            gene_dict$set(Input.net[index, 1], "");
          }
          params <- rbind(params,data.frame("index" = index, "Gene_A" = Input.net[index, 2], "Gene_B" = Input.net[index, 1], "npermut" = npermut, "nproc" = nproc, "min.size" = min.size, "max.size" = max.size, "seed" = seed_set, "include_cis_analysis" = include_cis_analysis))
        }
        if(!arguments$directed & !arguments$reversedirected){
          if(gene_dict$has(Input.net[index, 2])){
              include_cis_analysis <- FALSE
          }else{
              include_cis_analysis <- TRUE
              gene_dict$set(Input.net[index, 2], "");
          }
          params <- rbind(params,data.frame("index" = index, "Gene_A" = Input.net[index, 1], "Gene_B" = Input.net[index, 2], "npermut" = npermut, "nproc" = nproc, "min.size" = min.size, "max.size" = max.size, "seed" = seed_set, "include_cis_analysis" = include_cis_analysis))
          
          if(gene_dict$has(Input.net[index, 1])){
              include_cis_analysis <- FALSE
          }else{
              include_cis_analysis <- TRUE
              gene_dict$set(Input.net[index, 1], "");
          }
          params <- rbind(params,data.frame("index" = index, "Gene_A" = Input.net[index, 2], "Gene_B" = Input.net[index, 1], "npermut" = npermut, "nproc" = nproc, "min.size" = min.size, "max.size" = max.size, "seed" = seed_set, "include_cis_analysis" = include_cis_analysis))
        }
    }
}




CSEA_results <- NULL

if(arguments$stage1){
    cat("Step 3 of 7, Performing CSEA analysis. It may take a long time.\n", file = stdout())
    
    if(arguments$core < 2){
      for(index in 1 : nrow(params)){
          result <- CSEA(params[index, "index"], params[index, "Gene_A"], params[index, "Gene_B"], params[index, "npermut"], params[index, "nproc"], params[index, "min.size"], params[index, "max.size"], params[index, "seed"], params[index, "include_cis_analysis"])
          if(nrow(result) > 0){
              CSEA_results <- rbind(CSEA_results, result)
          }
      }
    }else{
        if(arguments$core > detectCores()){
          stop("Error: specified CPU cores are larger than available!")
        }
        
        options(cl.cores = arguments$core)
        this.cluster <- makeCluster(getOption("cl.cores", 2))
        clusterCall(cl = this.cluster, fun = function(){
          library(data.table);
          library(fgsea)
        })
        if(formalArgs(clusterExport)[2] %in% "list"){
            clusterExport(cl = this.cluster,  list = c("CSEA", "params", "exp_tpm_dt", "mutations_dt", "mut_types"))
        }else{
            clusterExport(cl = this.cluster,varlist=c("CSEA", "params", "exp_tpm_dt", "mutations_dt", "mut_types"))
        }
    
        CSEA_results1 <- 
        parLapply(cl = this.cluster,
          1:nrow(params),
          function(idx) {
            args1 <- as.list(params[idx,])
            formals(CSEA) <- args1
            rtn <- replicate(1, CSEA(), simplify = FALSE)
            return(rtn)
          }
        )
        stopCluster(this.cluster)
        
        for(j in 1:length(CSEA_results1)){
          CSEA_results <- rbind(CSEA_results, CSEA_results1[[j]][[1]])
        }
    }

}else{
  cat("Step 2 and 3 of 7 skipped, because CSEA analysis is disabled.\n", file = stdout())
}


if(arguments$stage1){
    CSEA_results <- CSEA_results[!duplicated(CSEA_results[, c("Gene_A", "Gene_B", "Mutation_type")]),]
    
    
    CSEA_results$P_val_adjusted <- p.adjust(CSEA_results$P_val, method = "BH")
    CSEA_results$ID <- 1 : nrow(CSEA_results)
    
    CSEA_results <- CSEA_results[, c("ID", "index", "Gene_A", "Gene_B", "Mutation_type", "Regulation", "ES", "NES", "P_val", "Leading_edge", "P_val_adjusted")]
    
    write.table(CSEA_results, file = paste0("./",arguments$folder,"/", arguments$name, "CSEA_results.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    rm(params)
}

# 2.2 Mediation
if(arguments$stage2){
    lm_fit_fast <- function(dat, X, out) {
        tryCatch({
          beta = solve(crossprod(X), t(crossprod(dat, X)))
          return(beta[out, 1])
        },
        error = function(e){
          return(NA)
        })
    }
    
    freedman_lane_sim <- function(adj.design, ie.design, alpha2.design,
                                  gamma2.design,gamma3.design,delta.design,
                                  adj.fitted, adj.reduced.fitted, 
                                  adj.residual, adj.reduced.residual,
                                  ie.fitted, ie.reduced.fitted, ie.residual,
                                  ie.reduced.residual, alpha2.fitted,
                                  alpha2.reduced.fitted, alpha2.residual,
                                  alpha2.reduced.residual,
                                  gamma2.fitted, gamma2.reduced.fitted,
                                  gamma2.residual, gamma2.reduced.residual,
                                  gamma3.fitted, gamma3.reduced.fitted,
                                  gamma3.residual, gamma3.reduced.residual,
                                  delta.fitted, delta.residual,
                                  n.greater.med, n.greater.med.Mstar,
                                  n.greater.alpha1, n.greater.alpha2,
                                  n.greater.gamma1, n.greater.gamma2,
                                  n.greater.gamma3, n.greater.delta, n.perms,
                                  adj.stat, ie.stat, alpha2.stat, gamma2.stat,
                                  gamma3.stat, delta.stat,
                                  gamma1.orig, gamma2.orig, gamma3.orig, 
                                  delta.orig, alpha1.orig, alpha2.orig,
                                  med.orig.lr) {

      # number of simulations with a stat greater than the observed.
      nc <- ncol(adj.reduced.residual)
    
      adj.residual  <- t(adj.residual)
      adj.fitted <- t(adj.fitted)
    
      adj.reduced.residual  <- t(adj.reduced.residual)
      adj.reduced.fitted <- t(adj.reduced.fitted)
    
      ie.residual  <- t(ie.residual)
      ie.fitted <- t(ie.fitted)
    
      ie.reduced.residual  <- t(ie.reduced.residual)
      ie.reduced.fitted <- t(ie.reduced.fitted)
      
      alpha2.residual  <- t(alpha2.residual)
      alpha2.fitted <- t(alpha2.fitted)
    
      alpha2.reduced.residual  <- t(alpha2.reduced.residual)
      alpha2.reduced.fitted <- t(alpha2.reduced.fitted)
      
      gamma2.residual  <- t(gamma2.residual)
      gamma2.fitted <- t(gamma2.fitted)
    
      gamma2.reduced.residual  <- t(gamma2.reduced.residual)
      gamma2.reduced.fitted <- t(gamma2.reduced.fitted)
      
      gamma3.residual  <- t(gamma3.residual)
      gamma3.fitted <- t(gamma3.fitted)
    
      gamma3.reduced.residual  <- t(gamma3.reduced.residual)
      gamma3.reduced.fitted <- t(gamma3.reduced.fitted)
      
      delta.residual  <- t(delta.residual)
      delta.fitted <- t(delta.fitted)
    
      for(i in 1:n.perms){
    
          	ystar.null <- NULL
          	mstar.null <- NULL
          	alpha2.null <- NULL
          	gamma1.null <- NULL
          	gamma1.alt  <- NULL
          	gamma2.null <- NULL
          	gamma3.null <- NULL
          	deltastar.null <- NULL

            set.seed(i)
          	perm <- sample(1:nc)
          
          	for(j in 1:nc){
          	      ystar.null[j] <- adj.reduced.fitted[j] + adj.reduced.residual[perm[j]]
          		    mstar.null[j] <- ie.reduced.fitted[j] + ie.reduced.residual[perm[j]]
          		    alpha2.null[j]  <- alpha2.reduced.fitted[j] + alpha2.reduced.residual[perm[j]]
          		    
          		    gamma2.null[j] <- gamma2.reduced.fitted[j] + gamma2.reduced.residual[perm[j]]
          		    gamma3.null[j] <- gamma3.reduced.fitted[j] + gamma3.reduced.residual[perm[j]]
                  deltastar.null[j] <- delta.residual[perm[j]]

          		}
          
          	mstar <- as.matrix(mstar.null)
          	colnames(mstar) <- "m"
          	adj.design.Mstar <- cbind(ie.design, as.matrix(mstar)) # V1, x1,x2,m
               
          	gamma1.star.null <- lm_fit_fast(as.matrix(ystar.null), adj.design, adj.stat)
          	gamma1.star.null.Mstar 	<- lm_fit_fast(as.matrix(ystar.null), adj.design.Mstar, adj.stat)

          	alpha1.star.null <- lm_fit_fast(as.matrix(mstar.null), ie.design, ie.stat)

          	alpha2.star.null <- lm_fit_fast(as.matrix(alpha2.null), alpha2.design, alpha2.stat)

          	gamma2.star.null <- lm_fit_fast(as.matrix(gamma2.null), gamma2.design, gamma2.stat)

          	gamma3.star.null <- lm_fit_fast(as.matrix(gamma3.null), gamma3.design, gamma3.stat)
          	
          	delta.star.null <- lm_fit_fast(as.matrix(deltastar.null), delta.design, delta.stat)
          
          	medstar.null.lr <- alpha1.star.null * gamma1.star.null
          	medstar.null.lr.Mstar <- alpha1.star.null * gamma1.star.null.Mstar 
          	
          	n.greater.med <- n.greater.med + (abs(medstar.null.lr) >= abs(med.orig.lr))
          	n.greater.med.Mstar	<- n.greater.med.Mstar + (abs(medstar.null.lr.Mstar) >= abs(med.orig.lr))
          	if(length(alpha1.star.null) > 0 & !is.na(alpha1.star.null)){
          	    n.greater.alpha1 <- n.greater.alpha1 + (abs(alpha1.star.null) >= abs(alpha1.orig))
          	}else{
          	    n.greater.alpha1 <- n.greater.alpha1 + 1
          	}
          	if(length(alpha2.star.null) > 0 & !is.na(alpha2.star.null)){
          	    n.greater.alpha2 <- n.greater.alpha2 + (abs(alpha2.star.null) >= abs(alpha2.orig))
          	}else{
          	    n.greater.alpha2 <- n.greater.alpha2 + 1
          	}
          	if(length(gamma1.star.null) > 0 & !is.na(gamma1.star.null)){
          	    n.greater.gamma1 <- n.greater.gamma1 + (abs(gamma1.star.null) >= abs(gamma1.orig))
          	}else{
          	    n.greater.gamma1 <- n.greater.gamma1 + 1
          	}
          	if(length(gamma2.star.null) > 0 & !is.na(gamma2.star.null)){
          	    n.greater.gamma2 <- n.greater.gamma2 + (abs(gamma2.star.null) >= abs(gamma2.orig))
          	}else{
          	    n.greater.gamma2 <- n.greater.gamma2 + 1
          	}
            if(length(gamma3.star.null) > 0 & !is.na(gamma3.star.null)){
          	    n.greater.gamma3 <- n.greater.gamma3 + (abs(gamma3.star.null) >= abs(gamma3.orig))
            }else{
          	    n.greater.gamma3 <- n.greater.gamma3 + 1
          	}
            if(length(delta.star.null) > 0 & !is.na(delta.star.null)){
          	    n.greater.delta <- n.greater.delta + (abs(delta.star.null) >= abs(delta.orig))
            }else{
          	    n.greater.delta <- n.greater.delta + 1
          	}
            
          	#medstar.alt.lr[i] <- medstar.alt.lr.i
          	#medstar.alt.lr.Mstar[i] <- medstar.alt.lr.i.Mstar
       }
    
    
      return(list("n.greater.med" = n.greater.med, "n.greater.med.Mstar" = n.greater.med.Mstar,	"n.greater.alpha1" = n.greater.alpha1, "n.greater.alpha2" = n.greater.alpha2, "n.greater.gamma1" = n.greater.gamma1, "n.greater.gamma2" = n.greater.gamma2, "n.greater.gamma3" = n.greater.gamma3, "n.greater.delta" = n.greater.delta))
    
    }
    
    # Freedman Lane Function of smaller p-value
    
    freedman_lane_sim_refine_mediation_p_value <- function(ie.design,
    			 	    adj.reduced.fitted, adj.reduced.residual,
    			 	    ie.reduced.fitted, ie.reduced.residual,
    				    n.greater.med.Mstar, n.greater.alpha1,
    				    n.perms, adj.stat, ie.stat,
    				    gamma1.orig, med.orig.lr) {
     
      # number of simulations with a stat greater than the observed.
      nc <- ncol(adj.reduced.residual)
    
      adj.reduced.residual  <- t(adj.reduced.residual)
      adj.reduced.fitted <- t(adj.reduced.fitted)
    
      ie.reduced.residual  <- t(ie.reduced.residual)
      ie.reduced.fitted <- t(ie.reduced.fitted)
    
      for(i in 1:n.perms){
    
          	ystar.null <- NULL
          	mstar.null <- NULL
    
          	gamma1.null <- NULL
          	gamma1.alt  <- NULL
          	
            set.seed(i)
          	perm <- sample(1 : nc)
          
          	for(j in 1:nc){
          	      ystar.null[j] <- adj.reduced.fitted[j] + adj.reduced.residual[perm[j]]
          		    mstar.null[j] <- ie.reduced.fitted[j] + ie.reduced.residual[perm[j]]
          	}
          
          	mstar <- as.matrix(mstar.null)
          	colnames(mstar) <- "m"
          	adj.design.Mstar <- cbind(ie.design, as.matrix(mstar)) # V1, x1,x2,m
    
          	gamma1.star.null.Mstar 	<- lm_fit_fast(as.matrix(ystar.null), adj.design.Mstar, adj.stat)
    
          	alpha1.star.null <- lm_fit_fast(as.matrix(mstar.null), ie.design, ie.stat)
    
          	medstar.null.lr.Mstar <- alpha1.star.null * gamma1.star.null.Mstar 
    
          	n.greater.med.Mstar	<- n.greater.med.Mstar + (abs(medstar.null.lr.Mstar) >= abs(med.orig.lr))
          	
       }
    
      return(n.greater.med.Mstar)
    
    }
    
    # Freedman Lane Function of smaller p-value
    
    freedman_lane_sim_refine_delta_p_value <- function(delta.design,
    			 	    delta.residual, n.greater.delta, n.perms, delta.stat, delta.orig) {
    
      # number of simulations with a stat greater than the observed.
      nc <- ncol(delta.residual)
    
      delta.residual  <- t(delta.residual)
    
      for(i in 1:n.perms){
    
          	deltastar.null <- NULL
          	
            set.seed(i)
          	perm <- sample(1 : nc)
          
          	for(j in 1:nc){ 
               deltastar.null[j] <- delta.residual[perm[j]]
          	}
          
          	delta.star.null <- lm_fit_fast(as.matrix(deltastar.null), delta.design, delta.stat)
    
            if(length(delta.star.null) > 0 & !is.na(delta.star.null)){
          	  n.greater.delta <- n.greater.delta + (abs(delta.star.null) >= abs(delta.orig))
            }else{
          	  n.greater.delta <- n.greater.delta + 1
          	}
       }
      return(n.greater.delta)
    }
    
    mediation_analysis <- 
      function(dat = NULL, n.perms = 100, Permut1 = FALSE, precise_pvalue = FALSE){
      
      y <- dat$y
      x1 <- dat$x1
      x2 <- dat$x2
      m <- dat$m

      # model fitting
      crude.lr <- lm(y ~ x2 + x1, data = dat, model = FALSE)
      adj.lr   <- update(crude.lr, formula = . ~ m + .)
      ie.lr    <- lm(m ~ x1 + x2, data = dat, model = FALSE)
      x.lr <- lm(x2 ~ x1, data = dat, model = FALSE)
    
      ############################################################
      # Data needs to be in matrix form for the following functions
      dat <- dat[, c(4, 1, 2, 3)]
      my.matrix <- as.matrix(dat)
      y <- as.matrix(my.matrix[, 1]) # y
      adj.variables_null <- as.matrix(my.matrix[, 2:3]) # x1, x2
      adj.variables  <- as.matrix(my.matrix[, 2:4]) # x1, x2, m
      m <- as.matrix(my.matrix[, 4]) # m
      ie.variables_null  <- as.matrix(my.matrix[, 3]) # x2
      ie.variables   <- as.matrix(my.matrix[, 2:3]) # x1, x2
      x1 <- as.matrix(my.matrix[, 2]) # x1
      x2 <- as.matrix(my.matrix[, 3]) # x2
      gamma2.variables <- as.matrix(my.matrix[, 3:4]) # x2, m
      gamma3.variables <- as.matrix(my.matrix[, c(2,4)]) # x1, m
      
      colnames(x1) <- "x1" # X
      colnames(y) <- "y" # Y
      colnames(m) <- "m" # M
      colnames(ie.variables_null) <- "x2" # C
    
      # Set design matrices
      # gamma 1
      adj.design <- cbind(1, adj.variables) # 1 x1, x2, m
      adj.reduced.design <- cbind(1, adj.variables_null) # 1 x1, x2
      
      # alpha 1
      ie.design <- cbind(1, ie.variables) # 1 x1, x2
      ie.reduced.design <- cbind(1, ie.variables_null) # 1 x2
      
      # alpha 2
      alpha2.design <- cbind(1, ie.variables) # 1 x1, x2
      alpha2.reduced.design <- cbind(1, x1) # 1 x1
      
      # gamma 2
      gamma2.design <- cbind(1, adj.variables) # 1 x1, x2, m
      gamma2.reduced.design <- cbind(1, gamma2.variables) # 1 x2, m
      
      # gamma 3
      gamma3.design <- cbind(1, adj.variables) # 1 x1, x2, m
      gamma3.reduced.design <- cbind(1, gamma3.variables) # 1 x1, m
      # delta
      delta.design <- cbind(1, x1) # 1 x1
    
      # Obtain fitted outcomes and residuals
      # gamma 1
      adj.fit <- lm.fit(adj.design, y)
      adj.fitted <- t(adj.fit$fitted.values)
      adj.residual <- t(adj.fit$residuals)
    
      adj.fit0 <- lm.fit(adj.reduced.design, y)
      adj.reduced.fitted <- t(adj.fit0$fitted.values)
      adj.reduced.residual  <- t(adj.fit0$residuals)
      
      # alpha 1
      ie.fit    <- lm.fit(ie.design, m)
      ie.fitted <- t(ie.fit$fitted.values)
      ie.residual  <- t(ie.fit$residuals)
    
      ie.fit0 <- lm.fit(ie.reduced.design, m)
      ie.reduced.fitted <- t(ie.fit0$fitted.values)
      ie.reduced.residual  <- t(ie.fit0$residuals)
      
      # alpha 2
      alpha2.fit    <- lm.fit(alpha2.design, m)
      alpha2.fitted <- t(alpha2.fit$fitted.values)
      alpha2.residual  <- t(alpha2.fit$residuals)
    
      alpha2.fit0 <- lm.fit(alpha2.reduced.design, m)
      alpha2.reduced.fitted <- t(alpha2.fit0$fitted.values)
      alpha2.reduced.residual  <- t(alpha2.fit0$residuals)
      
      #gamma 2
      gamma2.fit <- lm.fit(gamma2.design, y)
      gamma2.fitted <- t(gamma2.fit$fitted.values)
      gamma2.residual <- t(gamma2.fit$residuals)
      
      gamma2.reduced.fit <- lm.fit(gamma2.reduced.design, y)
      gamma2.reduced.fitted <- t(gamma2.reduced.fit$fitted.values)
      gamma2.reduced.residual <- t(gamma2.reduced.fit$residuals)
      
      #gamma 3
      gamma3.fit <- lm.fit(gamma3.design, y)
      gamma3.fitted <- t(gamma3.fit$fitted.values)
      gamma3.residual <- t(gamma3.fit$residuals)
      
      gamma3.reduced.fit <- lm.fit(gamma3.reduced.design, y)
      gamma3.reduced.fitted <- t(gamma3.reduced.fit$fitted.values)
      gamma3.reduced.residual <- t(gamma3.reduced.fit$residuals)
      
      #delta
      delta.fit <- lm.fit(delta.design, x2)
      delta.fitted <- t(delta.fit$fitted.values)
      delta.residual <- t(delta.fit$residuals)
    
      # Calculate test statistic
      adj.stat <- setdiff(colnames(adj.variables), colnames(adj.variables_null))
      gamma1.orig <- adj.fit$coef[adj.stat] #gamma 1
    
      ie.stat	<- setdiff(colnames(ie.variables), colnames(ie.variables_null))
      #  ie.stat	<- "x1"
      alpha1.orig 	<- ie.fit$coef[ie.stat] # alpha 1
      
      alpha2.stat	<- setdiff(colnames(ie.variables), colnames(x1))
      alpha2.orig 	<- alpha2.fit$coef[alpha2.stat] # alpha 2
      
      gamma2.stat <- setdiff(colnames(adj.variables), colnames(gamma2.variables)) # x2
      gamma2.orig 	<- ie.fit$coef[gamma2.stat] # gamma 2
      
      gamma3.stat <- setdiff(colnames(adj.variables), colnames(gamma3.variables)) # x1
      gamma3.orig 	<- ie.fit$coef[gamma3.stat] # gamma 3
      
      delta.stat <- colnames(x1) # x1
      delta.orig 	<- delta.fit$coef[delta.stat] # delta
    
      med.orig.lr <- alpha1.orig * gamma1.orig
    
      # Permutations 
      n.greater.med	<- 0
      n.greater.med.Mstar <- 0
      n.greater.alpha1	<- 0
      n.greater.alpha2	<- 0
      n.greater.gamma1	<- 0
      n.greater.gamma2	<- 0
      n.greater.gamma3	<- 0
      n.greater.delta	<- 0
      
      p.med.null <- NULL
      p.med.null.Mstar <- NULL
      p.alpha1.null <- NULL
      p.alpha2.null <- NULL
      p.gamma1.null <- NULL
      p.gamma2.null <- NULL
      p.gamma3.null <- NULL
      p.delta.null <- NULL
        
      if(Permut1){
        li <- freedman_lane_sim(adj.design, ie.design, alpha2.design, gamma2.design, gamma3.design, delta.design, adj.fitted, adj.reduced.fitted, adj.residual, adj.reduced.residual, ie.fitted, ie.reduced.fitted, ie.residual, ie.reduced.residual, alpha2.fitted, alpha2.reduced.fitted, alpha2.residual, alpha2.reduced.residual, gamma2.fitted, gamma2.reduced.fitted, gamma2.residual, gamma2.reduced.residual, gamma3.fitted, gamma3.reduced.fitted, gamma3.residual, gamma3.reduced.residual, delta.fitted, delta.residual, n.greater.med, n.greater.med.Mstar, n.greater.alpha1, n.greater.alpha2, n.greater.gamma1, n.greater.gamma2, n.greater.gamma3, n.greater.delta, n.perms, adj.stat, ie.stat, alpha2.stat, gamma2.stat, gamma3.stat, delta.stat, gamma1.orig, gamma2.orig, gamma3.orig, delta.orig, alpha1.orig, alpha2.orig, med.orig.lr)
      
        p.med.null <- li$n.greater.med / n.perms
        p.med.null.Mstar <- li$n.greater.med.Mstar / n.perms
        p.alpha1.null <- li$n.greater.alpha1 / n.perms
        p.alpha2.null <- li$n.greater.alpha2 / n.perms
        p.gamma1.null <- li$n.greater.gamma1 / n.perms
        p.gamma2.null <- li$n.greater.gamma2 / n.perms
        p.gamma3.null <- li$n.greater.gamma3 / n.perms
        p.delta.null <- li$n.greater.delta / n.perms
        if(precise_pvalue){
            if(!is.na(p.delta.null)){
        	    if(p.delta.null < 0.001){
        	    	refined_p_delta <- NA
        	    	refined_N_delta <- freedman_lane_sim_refine_delta_p_value(delta.design, delta.residual, n.greater.delta, n.perms * 10, delta.stat, delta.orig)
        	    	refined_p_delta <- (refined_N_delta) / (n.perms * 10)
        	    	if(refined_p_delta < 0.0001){
        	    			refined_N_delta <- freedman_lane_sim_refine_delta_p_value(delta.design, delta.residual, n.greater.delta, n.perms * 100, delta.stat, delta.orig)
        	    			refined_p_delta <- (refined_N_delta) / (n.perms * 100)
        	    	}
        	    	if(refined_p_delta > 0){
        	    		p.delta.null <- refined_p_delta
        	    	}
        	    	
        	    }
            }
            
            if(!is.na(p.med.null.Mstar)){
        	    if(p.med.null.Mstar < 0.001){
        	    	refined_p_med <- NA
        	    	refined_N_med <- freedman_lane_sim_refine_mediation_p_value(ie.design, adj.reduced.fitted, adj.reduced.residual, ie.reduced.fitted, ie.reduced.residual, n.greater.med.Mstar, n.greater.alpha1, n.perms * 10, adj.stat, ie.stat, gamma1.orig, med.orig.lr)
        				refined_p_med <- (refined_N_med) / (n.perms * 10)
        	    	if(refined_p_med < 0.0001){
        	    			refined_N_med <- freedman_lane_sim_refine_mediation_p_value(ie.design, adj.reduced.fitted, adj.reduced.residual, ie.reduced.fitted, ie.reduced.residual, n.greater.med.Mstar, n.greater.alpha1, n.perms * 100, adj.stat, ie.stat, gamma1.orig, med.orig.lr)
        					  refined_p_med <- (refined_N_med) / (n.perms * 100)
        	    	}
        	    	if(refined_p_med > 0){
        	    		p.med.null.Mstar <- refined_p_med
        	    	}
        	    	
        	    }
            }
        }
        
      }
    
      tmp<-data.frame("gamma_1" = coef(adj.lr)[2], "gamma_2" = coef(adj.lr)[3], "gamma_3" = coef(adj.lr)[4], "alpha_1" = coef(ie.lr)[2], "alpha_2" = coef(ie.lr)[3], "delta" = coef(x.lr)[2], "P_gamma_1" = p.gamma1.null, "P_gamma_2" = p.gamma2.null, "P_gamma_3" = p.gamma3.null, "P_alpha_1" = p.alpha1.null, "P_alpha_2" = p.alpha2.null, "P_delta" = p.delta.null, "mediation" = med.orig.lr, "P_mediation" = p.med.null.Mstar)

      return(tmp)
      }
    
    mediation_analysis_exp_dep <-
      function(Gene_A = "", Gene_B = "", Mut_type1 = "", nperm = 1000, ID = 1, Permut1 = FALSE, is_cellline_dat = TRUE, precise_pvalue = FALSE){

        res <- list("result" = data.frame("ID" = ID, "targetA" = Gene_A, "targetB" = Gene_B, "gamma_1" = NA, "gamma_2" = NA, "gamma_3" = NA, "alpha_1" = NA, "alpha_2" = NA, "delta" = NA, "P_gamma_1" = NA, "P_gamma_2" = NA, "P_gamma_3" = NA, "P_alpha_1" = NA, "P_alpha_2" = NA, "P_delta" = NA, "mediation" = NA, "P_mediation" = NA))
        
        expr_dat_gene_1 <- exp_tpm_dt[.(Gene_A)]
        expr_dat_gene_2 <- exp_tpm_dt[.(Gene_B)]
        if(nrow(expr_dat_gene_1) < 2 | nrow(expr_dat_gene_2) < 2){
            return(res)
        }
        if(!is_cellline_dat){
            expr_dat_gene_1 <- merge(expr_dat_gene_1, mapping, by.x="cell_line", by.y="ID", all.x=T)
            expr_dat_gene_1 <- expr_dat_gene_1[, c("gene", "rna_expression", "nearest_cell_line")]
            colnames(expr_dat_gene_1) <- c("gene", "rna_expression_1", "cell_line")
            expr_dat_gene_1 <- aggregate(rna_expression_1 ~ cell_line, data=expr_dat_gene_1, mean)
            colnames(expr_dat_gene_1) <- c("cell_line", "rna_expression_1")
            expr_dat_gene_1 <- as.data.table(expr_dat_gene_1)
            
            expr_dat_gene_2 <- merge(expr_dat_gene_2, mapping, by.x="cell_line", by.y="ID", all.x=T)
            expr_dat_gene_2 <- expr_dat_gene_2[, c("gene", "rna_expression", "nearest_cell_line")]
            colnames(expr_dat_gene_2) <- c("gene", "rna_expression_2", "cell_line")
            expr_dat_gene_2 <- aggregate(rna_expression_2 ~ cell_line, data=expr_dat_gene_2, mean)
            colnames(expr_dat_gene_2) <- c("cell_line", "rna_expression_2")
            expr_dat_gene_2 <- as.data.table(expr_dat_gene_2)
        }else{
            expr_dat_gene_1 <- expr_dat_gene_1[, c("rna_expression", "cell_line")]
            colnames(expr_dat_gene_1) <- c("rna_expression_1", "cell_line")
            expr_dat_gene_2 <- expr_dat_gene_2[, c("rna_expression", "cell_line")]
            colnames(expr_dat_gene_2) <- c("rna_expression_2", "cell_line")
        }
        
        dep_dat_gene_1 <- dependency_dt[.(Gene_A)]
        colnames(dep_dat_gene_1) <- c("gene_1", "dependency_1", "cell_line", "context")
        
        setkey(dep_dat_gene_1, 'cell_line')
        setkey(expr_dat_gene_1, 'cell_line')
    
        dep_dat_gene_2 <- dependency_dt[.(Gene_B)]
        
        colnames(dep_dat_gene_2) <- c("gene_2", "dependency_2", "cell_line", "context")
        
        setkey(dep_dat_gene_2, 'cell_line')
        setkey(expr_dat_gene_2, 'cell_line')
        
        DTlist <- list(expr_dat_gene_1, dep_dat_gene_1, expr_dat_gene_2, dep_dat_gene_2)
        merged_dat_dep_exp_gene_1_2 <-Reduce(function(X,Y) X[Y], DTlist)
        merged_dat_dep_exp_gene_1_2 <- merged_dat_dep_exp_gene_1_2[!is.na(merged_dat_dep_exp_gene_1_2$cell_line),]
        merged_dat_dep_exp_gene_1_2 <- merged_dat_dep_exp_gene_1_2 %>% filter(!is.na(rna_expression_1) & !is.na(rna_expression_2) & !is.na(dependency_1) & !is.na(dependency_2))

        mutations_gene1 <- mutations_dt[.(Gene_A), c("var_class", "cell_line")] # do not use the provided context
        if(!is_cellline_dat){ # For tumor data
            mutations_gene1 <- merge(mutations_gene1, mapping, by.x="cell_line", by.y="ID", all.x=T)
            mutations_gene1 <- mutations_gene1[!is.na(mutations_gene1$nearest_cell_line),]
            mutations_gene1 <- mutations_gene1[, c("var_class", "nearest_cell_line")]
            colnames(mutations_gene1) <- c("var_class", "cell_line")
        }
        
        if(Mut_type1 == "Any"){
            mut_ctx_1 <- mutations_gene1$cell_line
        }else{
            mut_ctx_1 <- mutations_gene1$cell_line[which(mutations_gene1$var_class %in% Mut_type1)]
        }

        if(nrow(merged_dat_dep_exp_gene_1_2) < 1){
           return(res)
        }
        
        if(length(which(merged_dat_dep_exp_gene_1_2$cell_line %in% mut_ctx_1)) > 30){ # At least 30 observations are required to build lr models, according to 10 observations for one variable
            merged_dat_dep_exp_gene_1_2 <- merged_dat_dep_exp_gene_1_2[which(merged_dat_dep_exp_gene_1_2$cell_line %in% mut_ctx_1),] # keep only contexts matching those for genetic alterations
        }else{
            
            ctx <- sapply(mut_ctx_1, function(z) paste(strsplit(z, split="_")[[1]][-1], collapse="_"))
              
            freq <- table(ctx)
            
            if(length(which(merged_dat_dep_exp_gene_1_2$context %in% names(freq[which(freq==max(freq))]))) > length(which(merged_dat_dep_exp_gene_1_2$cell_line %in% mut_ctx_1))){
                merged_dat_dep_exp_gene_1_2 <- merged_dat_dep_exp_gene_1_2[which(merged_dat_dep_exp_gene_1_2$context %in% names(freq[which(freq==max(freq))])),] # keep only top 25% contexts
            }else{
                merged_dat_dep_exp_gene_1_2 <- merged_dat_dep_exp_gene_1_2[which(merged_dat_dep_exp_gene_1_2$cell_line %in% mut_ctx_1),] # keep only contexts matching those for genetic alterations
            }
        }
       
        if(nrow(merged_dat_dep_exp_gene_1_2)<1){
            return(res)
        }
        
        # data merged_dat_dep_exp_gene_1_2
        # col 1: rna_expression # gene 1's expr
        # col 2: dependency # gene 1's dep
        # col 3: i.rna_expression # gene 2's expr
        # col 4: i.dependency # gene 2's dep
    
        merged_dat_dep_exp_gene_1_2 <- merged_dat_dep_exp_gene_1_2[, c("rna_expression_1", "dependency_1", "rna_expression_2", "dependency_2")]
        colnames(merged_dat_dep_exp_gene_1_2) <- c("x1", "x2", "m", "y")
        if(nrow(merged_dat_dep_exp_gene_1_2) > 1){ # lr requires at least two records
            result_gene1_to_2 <- mediation_analysis(dat = merged_dat_dep_exp_gene_1_2, n.perms=nperm, Permut1 = Permut1, precise_pvalue)
            res <- list("result"=data.frame("ID" = ID, "targetA" = Gene_A , "targetB" = Gene_B, "gamma_1" = result_gene1_to_2[1], "gamma_2" = result_gene1_to_2[2], "gamma_3" = result_gene1_to_2[3], "alpha_1" = result_gene1_to_2[4], "alpha_2" = result_gene1_to_2[5], "delta" = result_gene1_to_2[6], "P_gamma_1" = result_gene1_to_2[7], "P_gamma_2" = result_gene1_to_2[8], "P_gamma_3" = result_gene1_to_2[9], "P_alpha_1" = result_gene1_to_2[10], "P_alpha_2" = result_gene1_to_2[11], "P_delta" = result_gene1_to_2[12], "mediation" = result_gene1_to_2[13], "P_mediation" = result_gene1_to_2[14]))
            return(res)
        }else{
            return(res)
        }
      }
}

final_med_results <- data.frame("ID" = vector("numeric",length=0), "targetA" = vector("character", length=0), "targetB" = vector("character",length=0), "gamma_1" = vector("numeric",length=0), "gamma_2" = vector("numeric",length=0), "gamma_3" = vector("numeric",length=0), "alpha_1" = vector("numeric",length=0), "alpha_2" = vector("numeric",length=0), "delta" = vector("numeric",length=0), "P_gamma_1" = vector("numeric",length=0), "P_gamma_2" = vector("numeric",length=0), "P_gamma_3" = vector("numeric",length=0), "P_alpha_1" = vector("numeric",length=0), "P_alpha_2" = vector("numeric",length=0), "P_delta" = vector("numeric",length=0), "mediation" = vector("numeric",length=0), "P_mediation" = vector("numeric",length=0))

if(!exists("CSEA_results")){
    if(file.exists(paste0("./",arguments$folder,"/", arguments$name, "CSEA_results.tsv"))){
        CSEA_results <- fread(paste0("./",arguments$folder,"/", arguments$name, "CSEA_results.tsv"), header =TRUE, sep ='\t')
    }else{
        stop("Stage 2 relies on results from stage 1, please perform stage 1 CSEA analysis first!")
    }
}

records <- data.frame(CSEA_results[, "Gene_A"], CSEA_results[, "Gene_B"], CSEA_results[, "Mutation_type"], "nperm" = 1000, "ID" = CSEA_results[, "ID"], "Permut1" = TRUE, "cellLine_dat" = arguments$cl, "precise_med_p" = arguments$precise_mediation_pvalue)
colnames(records) <- c("Gene_A", "Gene_B", "Mut_type1", "nperm", "ID", "Permut1",  "cellLine_dat",  "precise_med_p")

if(arguments$stage2){
    cat("Step 4 of 7, Performing mediation analysis. It may take a long time.\n", file = stdout())
    
    if(arguments$core < 2){
        tmp_results <- NULL
        
        for(idx in 1 : nrow(records)){
            med <- mediation_analysis_exp_dep(records[idx, 1], records[idx, 2], records[idx, 3], records[idx, 4], records[idx, 5], records[idx, 6], records[idx, 7], records[idx, 8])
            tmp_results <- rbind(tmp_results, med[["result"]])
            if((idx %% 100) == 0){
              # at most 100 records are stored in tmp_results to keep the data frame not two big
              final_med_results <- rbind(final_med_results, tmp_results)
        	    
              tmp_results <- NULL
            }
        }
        final_med_results <- rbind(final_med_results, tmp_results)
    }else{
        if(arguments$core > detectCores()){
          stop("Error: specified CPU cores are larger than available!")
        }
      
        options(cl.cores = arguments$core)
        this.cluster <- makeCluster(getOption("cl.cores", 2))
        clusterCall(cl = this.cluster, fun = function(){
          library(dplyr);
          library(data.table)
        })
        if(formalArgs(clusterExport)[2] %in% "list"){
            clusterExport(cl = this.cluster, list = c("records", "mapping", "mediation_analysis_exp_dep", "mediation_analysis", "freedman_lane_sim", "lm_fit_fast", "dependency_dt", "freedman_lane_sim_refine_delta_p_value", "freedman_lane_sim_refine_mediation_p_value", "mutations_dt", "exp_tpm_dt"))
        }else{
            clusterExport(cl = this.cluster,varlist=c("records", "mapping", "mediation_analysis_exp_dep", "mediation_analysis", "freedman_lane_sim", "lm_fit_fast", "dependency_dt", "freedman_lane_sim_refine_delta_p_value", "freedman_lane_sim_refine_mediation_p_value", "mutations_dt", "exp_tpm_dt"),envir=environment())
        }
        
        mediation_analysis_results <- 
          parLapply(cl = this.cluster,
          1:nrow(records),
          function(idx) {
            #args1 <- as.list(records[idx,])
            #formals(mediation_analysis_exp_dep) <- args1
            #rtn <- replicate(1, mediation_analysis_exp_dep(records[idx,1],records[idx,2],records[idx,3],records[idx,4],records[idx,5],records[idx,6]), simplify = FALSE)
            rtn <- mediation_analysis_exp_dep(records[idx, 1], records[idx, 2], records[idx, 3], records[idx, 4], records[idx, 5], records[idx, 6], records[idx, 7], records[idx, 8])
            return(rtn)
          }
        )
        
        stopCluster(this.cluster)
        
        for( j in 1:length(mediation_analysis_results)){
          final_med_results <- rbind(final_med_results, mediation_analysis_results[[j]][["result"]])
        }
    }
    
    write.table(final_med_results, file = paste0("./",arguments$folder,"/", arguments$name, "mediation_analysis_results.tsv"), sep = "\t", row.names = FALSE, col.names=TRUE, quote=FALSE)
}else{
    cat("Step 4 of 7 is skipped, because mediation analysis is disabled.\n", file = stdout())
}

if(!exists("final_med_results")){
    final_med_results <- fread(paste0("./",arguments$folder,"/", arguments$name, "mediation_analysis_results.tsv"), header =TRUE, sep ='\t')
    strcol<- c("ID", "gamma_1", "gamma_2", "gamma_3", "alpha_1", "alpha_2", "delta", "P_gamma_1", "P_gamma_2", "P_gamma_3", "P_alpha_1", "P_alpha_2", "P_delta", "mediation", "P_mediation")
    for (col in strcol) set(final_med_results, j = col, value = as.numeric(final_med_results[[col]]))
}

if(!exists("CSEA_results")){
    CSEA_results <- fread(paste0("./",arguments$folder,"/", arguments$name, "CSEA_results.tsv"), header =TRUE, sep ='\t')
}

merged_CSEA_Mediation <- NULL
if(length(which(duplicated(CSEA_results[,c("Gene_A", "Gene_B", "Mutation_type")]))) > 0){
    merged_CSEA_Mediation <- merge(CSEA_results,final_med_results, by='ID', all.x=T)
    merged_CSEA_Mediation <- merged_CSEA_Mediation[which(!duplicated(CSEA_results[,c("Gene_A", "Gene_B", "Mutation_type")])),]
    cat(paste0(length(which(duplicated(CSEA_results[,c("Gene_A", "Gene_B", "Mutation_type")])))," records are duplicated and will be deleted\n"),file=stdout())
		merged_CSEA_Mediation$ID <- 1 : nrow(merged_CSEA_Mediation)
}else{
    merged_CSEA_Mediation <- merge(CSEA_results, final_med_results, by='ID', all.x=T)   
		merged_CSEA_Mediation$ID <- 1 : nrow(merged_CSEA_Mediation)
}

if(is.null(merged_CSEA_Mediation)){
		stop(paste0("No result is available for filtering."))
}

merged_CSEA_Mediation$targetA <- NULL
merged_CSEA_Mediation$targetB <- NULL
rm(dependency_dt)

######### Part 3: Filtering data #########

cat("Step 5 of 7, Filtering trans regulation results that can be explained by cis regulation.\n", file = stdout())
merged_CSEA_Mediation_pvals <- merged_CSEA_Mediation[, c("Gene_A", "Gene_B", "Mutation_type", "P_val_adjusted")]

merged_CSEA_Mediation_pvals_cis <- merged_CSEA_Mediation_pvals[which(merged_CSEA_Mediation_pvals$Gene_A == merged_CSEA_Mediation_pvals$Gene_B), c("Gene_A", "Gene_B", "Mutation_type", "P_val_adjusted")]

removed_trans_line_index <- NULL
cis_pvalue <- function(Gene_B = NULL, Mutation_type = NULL) {
    index_select <- which(merged_CSEA_Mediation_pvals_cis$Gene_B %in% Gene_B & merged_CSEA_Mediation_pvals_cis$Mutation_type %in% Mutation_type)
    cis_lines <- merged_CSEA_Mediation_pvals_cis[index_select,]
    if(nrow(cis_lines) < 1){
      return(NA)
    }else{
      return(cis_lines$P_val_adjusted)
    }
}

removed_trans_line_index <- NULL
if(arguments$remove_cis_confounding){
    final_cis_pvalue_results <- NULL
    if(arguments$core < 2){
      for(idx in 1:nrow(merged_CSEA_Mediation_pvals)){
        tmp <- cis_pvalue(merged_CSEA_Mediation_pvals[idx, c("Gene_B")], merged_CSEA_Mediation_pvals[idx, c("Mutation_type")])
        final_cis_pvalue_results <- c(final_cis_pvalue_results, tmp)
      }
    }else{
      options(cl.cores = arguments$core)
      this.cluster <- makeCluster(getOption("cl.cores", 2))
      clusterCall(cl = this.cluster, fun = function(){
        library(data.table)
      })
      if(formalArgs(clusterExport)[2] %in% "list"){
          clusterExport(cl = this.cluster, list = c("cis_pvalue", "merged_CSEA_Mediation_pvals", "merged_CSEA_Mediation_pvals_cis"))
      }else{
          clusterExport(cl = this.cluster,varlist=c("cis_pvalue", "merged_CSEA_Mediation_pvals", "merged_CSEA_Mediation_pvals_cis"))
      }
    
      cis_pval_results <- 
      parLapply(cl = this.cluster,
        1:nrow(merged_CSEA_Mediation_pvals),
        function(idx) {
          args1 <- as.list(merged_CSEA_Mediation_pvals[idx, c("Gene_B", "Mutation_type")])
          formals(cis_pvalue) <- args1
          rtn <- replicate(1, cis_pvalue(), simplify = FALSE)
          return(rtn)
        }
      )
      stopCluster(this.cluster)
      
      final_cis_pvalue_results <- unlist(cis_pval_results)
    }
    
    removed_trans_line_index <- which(final_cis_pvalue_results < 0.10 & merged_CSEA_Mediation_pvals$Gene_A != merged_CSEA_Mediation_pvals$Gene_B)
    if(length(removed_trans_line_index) > 0){
      cat(paste0(length(removed_trans_line_index), " of ", nrow(merged_CSEA_Mediation), " trans records were removed due to potential cis regulation\n"), file = stdout())
    }
}

CNA_confounder_exp_glm_analysis <-
  function(ID = NA, Gene_A = "", Gene_B = "", Mutation_type="", Permut1= 100, is_cellline_dat = FALSE){
      # Gene_A: gene with mutations
      # Gene_B: gene with CNA and exp
      # Mutation_type: gene 1 mutation type

    	exp_dat_gene_2 <- exp_tpm_dt[.(Gene_B),]
			exp_dat_gene_2 <- exp_dat_gene_2[!is.na(exp_dat_gene_2$cell_line), c("gene", "rna_expression", "cell_line")]

			CNV_dat_gene_2 <- CNA_dt[.(Gene_B),]
			CNV_dat_gene_2 <-CNV_dat_gene_2[!is.na(CNV_dat_gene_2$cell_line), c("log_copy_number", "cell_line")]
			colnames(CNV_dat_gene_2) <- c("CNV", "cell_line")
			merged_dat_CNA_exp_gene_2 <- merge(exp_dat_gene_2, CNV_dat_gene_2, by = "cell_line")


			mutations_dat <- mutations_dt[.(Gene_A), c("gene_name", "var_class", "cell_line")]
			merged_dat_CNA_exp_mut_gene_2 <- merge(merged_dat_CNA_exp_gene_2, mutations_dat, by = "cell_line", all=T)

			if(!is_cellline_dat){ # For tumor data
			  mutations_dat <- merge(merged_dat_CNA_exp_mut_gene_2, mapping, by.x = "cell_line", by.y = "ID", all.x = T)
			  mutations_dat <- mutations_dat[!is.na(mutations_dat$nearest_cell_line),]
			  mutations_dat <- mutations_dat[, c("rna_expression", "CNV", "var_class", "nearest_cell_line")]
			  colnames(mutations_dat) <- c("rna_expression", "CNV" , "var_class", "cell_line")
			}else{
			  mutations_dat <- merged_dat_CNA_exp_mut_gene_2[, c("rna_expression", "CNV", "var_class", "cell_line")]
			}
			mutations_dat$Mutation <- 0
			if(Mutation_type == "Any"){
			    #mutations_dat <- mutations_dat[, c("rna_expression","CNV" , "cell_line")]
			    mutations_dat[!is.na(mutations_dat$var_class),] <- 1
			}else{
			    #mutations_dat <- mutations_dat[mutations_dat$var_class %in% Mutation_type, c("rna_expression", "CNV", "cell_line")]
			    mutations_dat$Mutation[mutations_dat$var_class %in% Mutation_type] <- 1
			}

			colnames(mutations_dat) <- c("rna_expression", "CNV", "var_class", "cell_line", "Mutation")
      
      if(nrow(mutations_dat) < 1){
          tmp_data_frame <- data.frame("Gene_A" = Gene_A, "Gene_B" = Gene_B, "Mutation_type" = Mutation_type, "P_value_CNV" = NA)
          return(tmp_data_frame)
      }
    	
      if(arguments$cna_discretization){
          CN <- length(nrow(mutations_dat))
          CN[mutations_dat$CNV >= arguments$cna_amplification_cutoff] <- 1
          CN[mutations_dat$CNV <= arguments$cna_deletion_cutoff] <- -1
          CN[mutations_dat$CNV > arguments$cna_deletion_cutoff & mutations_dat$CNV < arguments$cna_amplification_cutoff] <- 0
          mutations_dat$CNV <- CN # log2CN to binary CN
      }

      model <- glm(rna_expression ~ CNV + Mutation, data = mutations_dat)
      coef_mut <- coef(summary(model))[, 1]["Mutation"]
      coef_CN <- coef(summary(model))[, 1]["CNV"]
      
      smaller_num <- 0
      
      P_permutation <- NA
      
      if(coef_CN >= 0) {
          model_mut_noCNV <- glm(rna_expression ~ Mutation, data = mutations_dat)
          coef_mut_noCNV <- coef(summary(model_mut_noCNV))[, 1]["Mutation"]
          decrease_rate <- abs((coef_mut_noCNV - coef_mut)/coef_mut) # https://nbisweden.github.io/excelerate-scRNAseq/session-normalization/confounding-factors.pdf
          # higher decrease_rate indicates more contribution of CNV to expression
          for(index in 1 : Permut1){
              set.seed(index)
              perm <- sample(1 : nrow(mutations_dat))
              merged_dat1 <- mutations_dat[, c("rna_expression", "CNV")]
              merged_dat1$Mutation <- mutations_dat$Mutation[perm]
              model <- glm(rna_expression ~ CNV + Mutation, data = merged_dat1)
              coef_mut <- coef(summary(model))[, 1]["Mutation"]
              model_mut_noCNV <- glm(rna_expression ~ Mutation, data = merged_dat1)
              coef_mut_noCNV <- coef(summary(model_mut_noCNV))[, 1]["Mutation"]
              decrease_rate_p <- (coef_mut_noCNV - coef_mut)/coef_mut
              if(abs(decrease_rate_p) <= decrease_rate){
                  smaller_num <- smaller_num + 1
              }
          }
          
          P_permutation <- smaller_num/Permut1

          if(P_permutation > 0.05 & P_permutation < 0.2){ # Further refine p-values if close to cutoff 0.1
              for(index in (1 + Permut1) : (Permut1 * 2)){
                  set.seed(index)
                  perm <- sample(1 : nrow(mutations_dat))
                  merged_dat1 <- mutations_dat[, c("rna_expression", "CNV")]
                  merged_dat1$Mutation <- mutations_dat$Mutation[perm]
                  model <- glm(rna_expression ~ CNV + Mutation, data = merged_dat1)
                  coef_mut <- coef(summary(model))[, 1]["Mutation"]
                  
                  model_mut_noCNV <- glm(rna_expression ~ Mutation, data = merged_dat1)
                  coef_mut_noCNV <- coef(summary(model_mut_noCNV))[, 1]["Mutation"]
                  decrease_rate_p <- (coef_mut_noCNV - coef_mut)/coef_mut
                  if(abs(decrease_rate_p) <= decrease_rate){
                      smaller_num <- smaller_num + 1
                  }
              }
              P_permutation <- smaller_num/(Permut1 * 2)
          }
      }

      # P-value indicates significance of CNA contribution to expression
      tmp_data_frame <- data.frame("ID" = ID, "P_value_CNV" = P_permutation)
   
      return(tmp_data_frame)
}

CNA_confounder_analysis_data_frame <- data.frame("ID" = rep(NA, nrow(merged_CSEA_Mediation)), "P_value_CNV" = rep(NA, nrow(merged_CSEA_Mediation)))
removed_index_CNA_confounder <- NULL
records <- NULL
cna_permutation_size <- arguments$cna_permutation_size

if(arguments$cna_adjust){
    if(file.exists(arguments$cna)){
        CNA_dt <- fread(arguments$cna, header =TRUE, sep ='\t')
        if(!("cell_line" %in% colnames(CNA_dt)) & !("-" %in% colnames(CNA_dt))){
            stop(paste0("Data ", arguments$cna, " does not contain columns of cell_line or disease."))
        }
        setkey(CNA_dt, gene_name)
        if(nrow(Input.net) > 0){
          CNA_dt <- CNA_dt[.(unique(c(Input.net$V1, Input.net$V2)))]
        }
        setkey(CNA_dt, gene_name)
    }else{
        stop(paste0("File path ", arguments$cna, " is not correctly specified."))
    }
    cat("Read", paste0(arguments$cna, " with ", nrow(CNA_dt), " rows after filtering.\n"), file = stdout())

  for(index in 1 : nrow(merged_CSEA_Mediation)){

      records <- rbind(records, data.frame("ID" = merged_CSEA_Mediation[index, "ID"], "Gene_A"  = merged_CSEA_Mediation[index, "Gene_A"], "Gene_B" = merged_CSEA_Mediation[index, "Gene_B"], "Mutation_type" = merged_CSEA_Mediation[index, "Mutation_type"], "Permut1" = cna_permutation_size, "is_cellline_dat" = arguments$cl))
  }

  cat("Step 6 of 7, removing CNA confounding.\n", file = stdout())
  if(arguments$core < 2){
      for(index in 1:nrow(records)){
          CNA_confounder_analysis_data_frame[index,] <- CNA_confounder_exp_glm_analysis("ID" = records[index, "ID"], "Gene_A"  = records[index, "Gene_A"], "Gene_B" = records[index, "Gene_B"], "Mutation_type" = records[index, "Mutation_type"], "Permut1" = records[index, "Permut1"], "is_cellline_dat" = records[index, "is_cellline_dat"])
      }
      removed_index_CNA_confounder <- CNA_confounder_analysis_data_frame$ID[which(CNA_confounder_analysis_data_frame$P_value_CNV < arguments$cna_adjust_Pvalue_cutoff)]
    }else{
      options(cl.cores = arguments$core)
      this.cluster <- makeCluster(getOption("cl.cores", 2))
      clusterCall(cl = this.cluster, fun = function(){
        library(data.table); library(dplyr)
      })
      if(formalArgs(clusterExport)[2] %in% "list"){
          clusterExport(cl = this.cluster, list = c("CNA_confounder_exp_glm_analy_dt", "exp_tpm_dt", "mutations_dt", "CNA_dt", "records", "arguments", "mapping"))
      }else{
          clusterExport(cl = this.cluster, varlist = c("CNA_confounder_exp_glm_analysis", "exp_tpm_dt", "mutations_dt", "CNA_dt", "records", "arguments", "mapping"))
      }
      
      CNA_confounder_pval_results <- 
      parLapply(cl = this.cluster,
          1:nrow(records),
          function(idx) {
              args1 <- as.list(c(records[idx, c("ID", "Gene_A", "Gene_B", "Mutation_type", "Permut1", "is_cellline_dat")]))
              formals(CNA_confounder_exp_glm_analysis) <- args1
              rtn <- replicate(1, CNA_confounder_exp_glm_analysis(), simplify = FALSE)
              return(rtn)
          }
      )
      stopCluster(this.cluster)
      for(j in 1:length(CNA_confounder_pval_results)){
          CNA_confounder_analysis_data_frame[j,] <- CNA_confounder_pval_results[[j]][[1]]
      }
      removed_index_CNA_confounder <- CNA_confounder_analysis_data_frame$ID[which(CNA_confounder_analysis_data_frame$P_value_CNV < arguments$cna_adjust_Pvalue_cutoff)]
      rm(CNA_confounder_pval_results)
  }
  removed_index_CNA_confounder_unique <- NULL
  
  if(length(removed_index_CNA_confounder) > 0){
      removed_index_CNA_confounder_unique <- setdiff(removed_index_CNA_confounder, merged_CSEA_Mediation_pvals$ID[removed_trans_line_index])
  }
  cat(paste0(length(removed_index_CNA_confounder_unique), " of ", nrow(CNA_confounder_analysis_data_frame), " records were removed due to potential CNA confounding\n"), file = stdout())
  rm(CNA_dt)

}
rm(records)

merged_CSEA_Mediation <- merged_CSEA_Mediation[!(merged_CSEA_Mediation$ID %in% unique(c(removed_index_CNA_confounder, merged_CSEA_Mediation_pvals$ID[removed_trans_line_index]))),]

if(arguments$remove_contradiction){
    removed_cis_line <- merged_CSEA_Mediation$ID[which(merged_CSEA_Mediation$Gene_A == merged_CSEA_Mediation$Gene_B & merged_CSEA_Mediation$delta > 0)]
    cat(paste0(length(removed_cis_line), " cis records were removed due to unreasonable mediation analysis results.\n"), file = stdout())
    
    removed_trans_line <- merged_CSEA_Mediation$ID[which(merged_CSEA_Mediation$Gene_A !=merged_CSEA_Mediation$Gene_B & merged_CSEA_Mediation$mediation > 0)]
    cat(paste0(length(removed_trans_line), " trans records were removed due to unreasonable mediation analysis results.\n"),file = stdout())
    
    merged_CSEA_Mediation <- merged_CSEA_Mediation[!(merged_CSEA_Mediation$ID %in% c(removed_cis_line, removed_trans_line)),]
    cat(paste0(nrow(merged_CSEA_Mediation), " of ",nrow(CNA_confounder_analysis_data_frame), " records were kept for p-value combination.\n"), file = stdout())
}
rm(exp_tpm_dt)
rm(CNA_confounder_analysis_data_frame)

write.table(merged_CSEA_Mediation, file = paste0("./",arguments$folder,"/", arguments$name, "merged_CSEA_Mediation.tsv"), sep = "\t", row.names = FALSE, col.names=TRUE, quote=FALSE)

######### Part 4: Combining results to scores #########

cat("Step 7 of 7, Combining the results and generating final results.\n", file = stdout())

liptak <- function(p, w) {
    if (missing(w)) {
        w <- rep(1, length(p))/length(p)
    } else {
        if (length(w) != length(p)){
            stop("Length of p and w must equal!")
        }
    }
    Zi <- -qnorm(p) 
    Z  <- sum(w*Zi)/sqrt(sum(w^2))
    p.val <- pnorm(-Z)
    return(p.val)
}

cis_line <- which(merged_CSEA_Mediation$Gene_A == merged_CSEA_Mediation$Gene_B)
trans_line <- which(merged_CSEA_Mediation$Gene_A != merged_CSEA_Mediation$Gene_B)
pvalues_cis_df <- merged_CSEA_Mediation[cis_line, c("P_val", "P_delta", "ID", "Mutation_type")]
pvalues_trans_df <- merged_CSEA_Mediation[trans_line, c("P_val", "P_mediation", "ID", "Mutation_type")]
colnames(pvalues_cis_df) <- c("p_CSEA", "p_cis", "ID", "Mut_type")
colnames(pvalues_trans_df) <- c("p_CSEA", "p_trans", "ID", "Mut_type")

pvalues_cis_df$p_CSEA[which(is.na(pvalues_cis_df$p_CSEA))] <- 1
pvalues_cis_df$p_CSEA[which(pvalues_cis_df$p_CSEA < 0.00001)] <- 0.00001
pvalues_cis_df$p_cis[which(is.na(pvalues_cis_df$p_cis))] <- 1
pvalues_cis_df$p_cis[which(pvalues_cis_df$p_cis < 0.00001)] <- 0.00001

pvalues_trans_df$p_CSEA[which(is.na(pvalues_trans_df$p_CSEA))] <- 1
pvalues_trans_df$p_CSEA[which(pvalues_trans_df$p_CSEA < 0.00001)] <- 0.00001
pvalues_trans_df$p_trans[which(is.na(pvalues_trans_df$p_trans))] <- 1
pvalues_trans_df$p_trans[which(pvalues_trans_df$p_trans < 0.00001)] <- 0.00001

if(!exists("mut_types")){
    mut_types <- unique(merged_CSEA_Mediation$Mutation_type)
    mut_types <- sort(mut_types)
}

weight_2s  <- c(1/10, 1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2, 1, seq(2, 10, by = 1))
params <- expand.grid(weight_2s, mut_types)
colnames(params) <- c("weight_2", "mut_type")

MLFC_metrics <- function(weight_2, mut_type) {
    
MLFC_dat <- NULL
    cat(paste0("Calculating ", mut_type, " ", weight_2, "\n"), file = stdout())
  
    pvalues_cis_df_tmp <- pvalues_cis_df[pvalues_cis_df$Mut_type %in% mut_type,]
  
    pvalues_cis_df_tmp$combined_P_vals <- apply(pvalues_cis_df_tmp[, c("p_CSEA", "p_cis")], 1, liptak, c(1, weight_2))
    
    pvalues_cis_df_tmp$combined_P_vals_adjusted <- p.adjust(pvalues_cis_df_tmp$combined_P_vals, method = "BH")
   
    pvalues_trans_df_tmp <- pvalues_trans_df[pvalues_trans_df$Mut_type %in% mut_type,]
    
    pvalues_trans_df_tmp$combined_P_vals <- apply(pvalues_trans_df_tmp[, c("p_CSEA", "p_trans")], 1, liptak, c(1, weight_2))
    
    pvalues_trans_df_tmp$combined_P_vals_adjusted <- p.adjust(pvalues_trans_df_tmp$combined_P_vals, method = "BH")
    
    #pvalues_cis_df_tmp <- pvalues_cis_df_tmp[pvalues_cis_df_tmp$combined_P_vals_adjusted >= cutoff_FDR,]
    pvalues_cis_df_tmp <- pvalues_cis_df_tmp[order(pvalues_cis_df_tmp$combined_P_vals),]
    
    running_sum <- 0
    for(i in 1:nrow(pvalues_cis_df_tmp)){
      pval_i <- pvalues_cis_df_tmp$combined_P_vals_adjusted[i]
      quantile_i <- i / nrow(pvalues_cis_df_tmp)
      tmp <- abs(log2(pval_i / quantile_i))
      if(!is.nan(tmp)){
        running_sum <- running_sum + tmp
      }
    }
    if(nrow(pvalues_cis_df_tmp) < 1){
    		MLFC_val_cis <- 0
    }else{
    		MLFC_val_cis <- running_sum / nrow(pvalues_cis_df_tmp)
    }
    MLFC_val_cis <- running_sum / nrow(pvalues_cis_df_tmp)
    
    # pvalues_trans_df_tmp <- pvalues_trans_df_tmp[pvalues_trans_df_tmp$combined_P_vals_adjusted >= cutoff_FDR,]
    pvalues_trans_df_tmp <- pvalues_trans_df_tmp[order(pvalues_trans_df_tmp$combined_P_vals),]
    
    running_sum <- 0
    if(nrow(pvalues_trans_df_tmp) > 0){
	    for(i in 1:nrow(pvalues_trans_df_tmp)){
	      pval_i <- pvalues_trans_df_tmp$combined_P_vals_adjusted[i]
	      quantile_i <- i / nrow(pvalues_trans_df_tmp)
	      tmp <- abs(log2(pval_i / quantile_i))
	      if(!is.nan(tmp)){
	        running_sum <- running_sum + tmp
	      }
	    }
    }
    if(nrow(pvalues_trans_df_tmp) < 1){
    		MLFC_val_trans <- 0
    }else{
    		MLFC_val_trans <- running_sum / nrow(pvalues_trans_df_tmp)
    }
    
    MLFC_dat <- data.frame("Mutation_type" = mut_type, "weight_2" = weight_2, "MLFC_cis" = MLFC_val_cis, "MLFC_trans" = MLFC_val_trans)
    return(MLFC_dat)
}

final_MLFC_results <- NULL

if(arguments$core < 2){
	for(i in 1:nrow(params)){
			
			MLFC_results <- MLFC_metrics(weight_2 = params[i, "weight_2"], mut_type = params[i, "mut_type"])
			final_MLFC_results <- rbind(final_MLFC_results, MLFC_results)
	}
}else{

		options(cl.cores = arguments$core)
		# set up a cluster using N-1 cores
		this.cluster <- makeCluster(getOption("cl.cores", 2))
		clusterCall(cl = this.cluster, fun = function(){

		})
		if(formalArgs(clusterExport)[2] %in% "list"){
		    clusterExport(cl = this.cluster,
		                  list = c("liptak", "pvalues_cis_df", 
		                           "pvalues_trans_df", "params", "MLFC_metrics"))
		}else{
		    clusterExport(cl = this.cluster, 
		                  varlist = c("liptak", "pvalues_cis_df", 
		                            "pvalues_trans_df", "params", "MLFC_metrics"))
		}

		MLFC_results <- 
		  parLapply(cl = this.cluster,
		  1:nrow(params),
		  function(idx) {
		    args <- as.list(params[idx,])
		    formals(MLFC_metrics) <- args
		    # the formals(), the list of arguments which control how you can call the function.
		    rtn <- replicate(1, MLFC_metrics(), simplify = FALSE)
		    # sims.per.seed: the number of replications to perform
		    # mediation_sim(): the expression that should be run repeatedly
		    # simplify, which contrlr the type of output the results of expr are saved into. 
		    # Use simplify = FALSE to get vectors saved into a list instead of in an array
		    #cat(paste0(idx, " ", "\n"),file = stdout())
		    return(rtn)
		  }
		)
		stopCluster(this.cluster)

		for( j in 1:length(MLFC_results)){
		  final_MLFC_results <- rbind(final_MLFC_results, MLFC_results[[j]][[1]])
		}
}


liptak_2 <- function(p) {
  w <- c(1, p[3])
  p <- p[1:2]
  if (missing(w)) {
    w <- rep(1, length(p)) / length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- -qnorm(p) 
  Z  <- sum(w*Zi) / sqrt(sum(w^2))
  p.val <- pnorm(-Z)
  return(p.val)
}

pvalues_cis_df$weight_2 <- 1

final_MLFC_results_cis <- final_MLFC_results[order(final_MLFC_results$MLFC_cis),]
final_MLFC_results_cis <- final_MLFC_results_cis[!duplicated(final_MLFC_results_cis$Mutation_type),]
final_MLFC_results_trans <- final_MLFC_results[order(final_MLFC_results$MLFC_trans),]
final_MLFC_results_trans <- final_MLFC_results_trans[!duplicated(final_MLFC_results_trans$Mutation_type),]
cat(paste0("Weights used for Liptak's p-value combination are chosen based on minimum MLFC metric.\n\n"), file = stdout())
cat(paste0("After calculation, the weights for different mutation types are shown below:\n"),file=stdout())
for(i in 1 : nrow(final_MLFC_results_cis)){
  pvalues_cis_df$weight_2[pvalues_cis_df$Mut_type == as.character(final_MLFC_results_cis$Mutation_type[i])] <- as.numeric(final_MLFC_results_cis$weight_2[i])
  cat(paste0("weights [CSEA, mediation] for cis ", final_MLFC_results_cis$Mutation_type[i], " mutations are [1, ", final_MLFC_results_cis$weight_2[i], "].\n"), file = stdout())
}

pvalues_trans_df$weight_2 <- 1

for(i in 1 : nrow(final_MLFC_results_trans)){
  pvalues_trans_df$weight_2[pvalues_trans_df$Mut_type == as.character(final_MLFC_results_trans$Mutation_type[i])] <- as.numeric(final_MLFC_results_trans$weight_2[i])
  cat(paste0("weights [CSEA, mediation] for trans ", final_MLFC_results_trans$Mutation_type[i], " mutations are [1, ", final_MLFC_results_trans$weight_2[i], "].\n"), file = stdout())
}

if(nrow(pvalues_cis_df) > 0){
		pvalues_cis_df$combined_P_vals <- apply(pvalues_cis_df[, c("p_CSEA", "p_cis", "weight_2")], 1 , liptak_2)
		pvalues_cis_df$combined_P_vals_adjusted <- p.adjust(pvalues_cis_df$combined_P_vals, method = "BH")

    ID_sig_cis <- pvalues_cis_df$ID[pvalues_cis_df$combined_P_vals_adjusted < arguments$FDR_cutoff_cis & pvalues_cis_df$p_CSEA < arguments$CSEA_cis_cutoff & pvalues_cis_df$p_cis < arguments$mediation_cis_cutoff] # 
    cat(paste0("Number of final cis genetic alterations is: ", length(ID_sig_cis), "\n"), file = stdout())
}

if(nrow(pvalues_trans_df) > 0){
    pvalues_trans_df$combined_P_vals <- apply(pvalues_trans_df[, c("p_CSEA", "p_trans", "weight_2")], 1, liptak_2)
		pvalues_trans_df$combined_P_vals_adjusted <- p.adjust(pvalues_trans_df$combined_P_vals, method = "BH")

    ID_sig_trans <- pvalues_trans_df$ID[pvalues_trans_df$combined_P_vals_adjusted < arguments$FDR_cutoff_trans & pvalues_trans_df$p_CSEA < arguments$CSEA_trans_cutoff & pvalues_trans_df$p_trans < arguments$mediation_trans_cutoff]
    cat(paste0("Number of final trans genetic alterations is: ", length(ID_sig_trans), "\n"), file = stdout())
}

merged_CSEA_Mediation_cis_sig <- merged_CSEA_Mediation[merged_CSEA_Mediation$ID %in% ID_sig_cis,]
merged_CSEA_Mediation_trans_sig <- merged_CSEA_Mediation[merged_CSEA_Mediation$ID %in% ID_sig_trans,]

for(i in 1 : length(mut_types)){
  cat(paste0("# of cis ", mut_types[i], " mutations is ", length(which(merged_CSEA_Mediation_cis_sig$Mutation_type == mut_types[i])), ".\n"), file = stdout())
  cat(paste0("# of trans ", mut_types[i], " mutations is ", length(which(merged_CSEA_Mediation_trans_sig$Mutation_type == mut_types[i])), ".\n"), file = stdout())
}

if(nrow(merged_CSEA_Mediation_cis_sig) > 0){
    merged_CSEA_Mediation_cis_sig <- merge(merged_CSEA_Mediation_cis_sig, pvalues_cis_df[, c("ID", "combined_P_vals", "combined_P_vals_adjusted")], by = 'ID', all.x = TRUE)
    merged_CSEA_Mediation_cis_sig$Effect_on_Cell_Proliferation[merged_CSEA_Mediation_cis_sig$ES > 0] <- "Increase_dependency"
    merged_CSEA_Mediation_cis_sig$Effect_on_Cell_Proliferation[merged_CSEA_Mediation_cis_sig$ES <= 0] <- "Decrease_dependency"
    colnames(merged_CSEA_Mediation_cis_sig) <- c("ID", "Line_in_PPI", "Gene_with_mutation", "Target_gene", "Mutation_type", "CSEA_Regulation", "CSEA_ES", "CSEA_NES", "CSEA_P_val", "CSEA_Leading_edge", "CSEA_P_val_adjusted", "gamma_1", "gamma_2", "gamma_3", "alpha_1", "alpha_2", "cis_effect", "P_val_gamma_1", "P_val_gamma_2", "P_val_gamma_3", "P_val_alpha_1", "P_val_alpha_2", "P_val_cis_effect", "mediation_effect", "P_val_mediation_effect", "combined_P_vals", "combined_P_vals_adjusted", "Effect_on_Cell_Proliferation")
  
    write.table(merged_CSEA_Mediation_cis_sig[, c("Line_in_PPI", "Gene_with_mutation", "Target_gene", "Mutation_type", "CSEA_ES", "CSEA_P_val", "CSEA_Leading_edge", "cis_effect", "P_val_cis_effect", "mediation_effect", "P_val_mediation_effect", "combined_P_vals", "combined_P_vals_adjusted", "Effect_on_Cell_Proliferation")], file = paste0("./",arguments$folder,"/", arguments$name, "cis_sig_genetic_alteration_target_results.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}else{
    write.table(paste0(c("ID", "Line_in_PPI", "Gene_with_mutation", "Target_gene", "Mutation_type", "CSEA_Regulation", "CSEA_ES", "CSEA_NES", "CSEA_P_val", "CSEA_Leading_edge", "CSEA_P_val_adjusted", "gamma_1", "gamma_2", "gamma_3", "alpha_1", "alpha_2", "cis_effect", "P_val_gamma_1", "P_val_gamma_2", "P_val_gamma_3", "P_val_alpha_1", "P_val_alpha_2", "P_val_cis_effect", "mediation_effect", "P_val_mediation_effect", "combined_P_vals", "combined_P_vals_adjusted", "Effect_on_Cell_Proliferation"), collapse="\t"), file = paste0("./",arguments$folder,"/", arguments$name, "cis_sig_genetic_alteration_target_results.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

if(nrow(merged_CSEA_Mediation_trans_sig) > 0){
    merged_CSEA_Mediation_trans_sig <- merge(merged_CSEA_Mediation_trans_sig,pvalues_trans_df[, c("ID", "combined_P_vals", "combined_P_vals_adjusted")], by='ID', all.x=T)
    merged_CSEA_Mediation_trans_sig$Effect_on_Cell_Proliferation[merged_CSEA_Mediation_trans_sig$ES > 0] <- "Increase_dependency"
    merged_CSEA_Mediation_trans_sig$Effect_on_Cell_Proliferation[merged_CSEA_Mediation_trans_sig$ES <= 0] <- "Decrease_dependency"
    
    colnames(merged_CSEA_Mediation_trans_sig) <- c("ID", "Line_in_PPI", "Gene_with_mutation", "Target_gene", "Mutation_type", "CSEA_Regulation", "CSEA_ES", "CSEA_NES", "CSEA_P_val", "CSEA_Leading_edge", "CSEA_P_val_adjusted", "gamma_1", "gamma_2", "gamma_3", "alpha_1", "alpha_2", "cis_effect", "P_val_gamma_1", "P_val_gamma_2", "P_val_gamma_3", "P_val_alpha_1", "P_val_alpha_2", "P_val_cis_effect", "mediation_effect", "P_val_mediation_effect", "combined_P_vals", "combined_P_vals_adjusted", "Effect_on_Cell_Proliferation")
  
    merged_CSEA_Mediation_trans_sig$cis_effect <- NA
    merged_CSEA_Mediation_trans_sig$P_val_cis_effect <- NA
    write.table(merged_CSEA_Mediation_trans_sig[, c("Line_in_PPI", "Gene_with_mutation", "Target_gene", "Mutation_type", "CSEA_ES", "CSEA_P_val", "CSEA_Leading_edge", "cis_effect", "P_val_cis_effect", "mediation_effect", "P_val_mediation_effect", "combined_P_vals", "combined_P_vals_adjusted", "Effect_on_Cell_Proliferation")], file = paste0("./",arguments$folder,"/", arguments$name, "trans_sig_genetic_alteration_target_results.tsv"), sep = "\t", row.names = FALSE, col.names=TRUE, quote=FALSE)
}else{
    write.table(paste0(c("ID", "Line_in_PPI", "Gene_with_mutation", "Target_gene", "Mutation_type", "CSEA_Regulation", "CSEA_ES", "CSEA_NES", "CSEA_P_val", "CSEA_Leading_edge", "CSEA_P_val_adjusted", "gamma_1", "gamma_2", "gamma_3", "alpha_1", "alpha_2", "cis_effect", "P_val_gamma_1", "P_val_gamma_2", "P_val_gamma_3", "P_val_alpha_1", "P_val_alpha_2", "P_val_cis_effect", "mediation_effect", "P_val_mediation_effect", "combined_P_vals", "combined_P_vals_adjusted", "Effect_on_Cell_Proliferation"), collapse = "\t"), file = paste0("./",arguments$folder,"/", arguments$name, "trans_sig_genetic_alteration_target_results.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

cat(paste0("Number of genes with cis mutations: ", length(unique(merged_CSEA_Mediation_cis_sig$Gene_with_mutation)), "; Number of genes with trans mutations: ", length(unique(merged_CSEA_Mediation_trans_sig$Gene_with_mutation)), "; Number of genes with cis/trans mutations:", length(unique(c(merged_CSEA_Mediation_cis_sig$Gene_with_mutation, merged_CSEA_Mediation_trans_sig$Gene_with_mutation))), "\n"), file = stdout())

if(arguments$output_mutation_site_centric_table){

    cat(paste0("Generating the genetic alteration tables.\n"), file = stdout())
    setkey(mutations_dt,gene_name, cell_line)
    
    CSEA_Mediation_cis_sig_anno_table <- NULL
    if(nrow(merged_CSEA_Mediation_cis_sig) > 0){
      for(index in 1:nrow(merged_CSEA_Mediation_cis_sig)){
        temp_line <- merged_CSEA_Mediation_cis_sig[index,]
        cls <- strsplit(temp_line$CSEA_Leading_edge, split=",")[[1]]
        if(index %% 100 == 0){
          cat(paste0(index, " / ", nrow(merged_CSEA_Mediation_cis_sig), "\n"), file = stdout())
        }
        CSEA_Mediation_cis_sig_anno_table_tmp <- NULL
        for(index2 in 1:length(cls)){
          dat_tmp <- mutations_dt[.(temp_line$Gene_with_mutation, cls[index2])]
      
          for(index3 in 1:nrow(dat_tmp)){
            if(nchar(dat_tmp[index3, "var_class"]) < 1){
              next
            }
            CSEA_Mediation_cis_sig_anno_table_tmp <- rbind(CSEA_Mediation_cis_sig_anno_table_tmp, data.frame("Line_in_PPI" = temp_line$Line_in_PPI, dat_tmp[index3,], "CSEA_ES" = temp_line$CSEA_ES, "CSEA_p_value" = temp_line$CSEA_P_val, "Mediation_p_value_cis" = temp_line$P_val_cis_effect, "Mediation_p_value_trans" = temp_line$P_val_mediation_effect, "Combined_p_value" = temp_line$combined_P_vals, "Combined_p_value_adjusted" = temp_line$combined_P_vals_adjusted, "Effect_on_Cell_Proliferation" = temp_line$Effect_on_Cell_Proliferation))
      
          }
        }
        CSEA_Mediation_cis_sig_anno_table <- rbind(CSEA_Mediation_cis_sig_anno_table, CSEA_Mediation_cis_sig_anno_table_tmp)
      }
      write.table(CSEA_Mediation_cis_sig_anno_table, file = paste0("./",arguments$folder,"/", arguments$name, "cis_sig_genetic_alteration_results.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }else{
      write.table("None", file = paste0("./",arguments$folder,"/", arguments$name, "cis_sig_genetic_alteration_results.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}

if(arguments$output_mutation_site_centric_table){
      
    CSEA_Mediation_trans_sig_anno_table <- NULL
    
    if(nrow(merged_CSEA_Mediation_trans_sig) > 0){
        for(index in 1 : nrow(merged_CSEA_Mediation_trans_sig)){
          temp_line <- merged_CSEA_Mediation_trans_sig[index,]
          cls <- strsplit(temp_line$CSEA_Leading_edge, split=",")[[1]]
          if(index %% 100 == 0){
            cat(paste0(index, " / ", nrow(merged_CSEA_Mediation_trans_sig), "\n"),file = stdout())
          }
          CSEA_Mediation_trans_sig_anno_table_tmp <- NULL
          for(index2 in 1:length(cls)){
            dat_tmp <- mutations_dt[.(temp_line$Gene_with_mutation, cls[index2])]
        
            for(index3 in 1:nrow(dat_tmp)){
              if(nchar(dat_tmp[index3, "var_class"]) < 1){
                next
              }
              CSEA_Mediation_trans_sig_anno_table_tmp <- rbind(CSEA_Mediation_trans_sig_anno_table_tmp, data.frame("Line_in_PPI" = temp_line$Line_in_PPI, dat_tmp[index3,], "CSEA_ES" = temp_line$CSEA_ES, "CSEA_p_value" = temp_line$CSEA_P_val, "Mediation_p_value_cis" = temp_line$P_val_cis_effect, "Mediation_p_value_trans" = temp_line$P_val_mediation_effect, "Combined_p_value" = temp_line$combined_P_vals, "Combined_p_value_adjusted" = temp_line$combined_P_vals_adjusted, "Effect_on_Cell_Proliferation" = temp_line$Effect_on_Cell_Proliferation))
        
            }
          }
          CSEA_Mediation_trans_sig_anno_table <- rbind(CSEA_Mediation_trans_sig_anno_table, CSEA_Mediation_trans_sig_anno_table_tmp)
        }
        write.table(CSEA_Mediation_trans_sig_anno_table, file = paste0("./",arguments$folder,"/", arguments$name, "trans_sig_genetic_alteration_results.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }else{
        write.table("None",file = paste0("./",arguments$folder,"/", arguments$name, "trans_sig_genetic_alteration_results.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}

cat("Complete.\n", file = stdout())