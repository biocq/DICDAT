# DICDAT

Command lines:

Rscript DICDAT.R --help

Usage: DICDAT.R [options]


Options:
	
	--folder=FOLDER
		Project folder to be used for storing data used and output by this program.

	--name=NAME
		Task name prefix.

	--stage1=STAGE1
		Whether perform CSEA analysis.

	--stage2=STAGE2
		Whether perform mediation analysis.

	--ga=GA
		Tabular file of genetic alteration data (i.e. mutations).

	--exp=EXP
		Tabular file of gene expression data.

	--dep=DEP
		Tabular file of gene dependency data [Optional]

	--cna=CNA
		Tabular file of CNV data [Optional]

	--cl=CL
		Whether the data is from DepMap project. [Optional]

	--map=MAP
		Tabular file of gene sample type mapper: map user-provided to in-house tissue types [Optional]

	--ppi=PPI
		Tabular file of two columns of genes, i.e. two-column tab delimited undirected PPI interactions. [Optional]

	--start_PPI=START_PPI
		Starting line of --ppi file. Set to -1 if no subsetting is needed [Optional]

	--end_PPI=END_PPI
		Ending line of --ppi file. Set to -1 if no subsetting is needed [Optional]

	--directed=DIRECTED
		If TRUE, only check the regulation of the gene in the first column upon the gene in the second column. [Optional]

	--reversedirected=REVERSEDIRECTED
		If TRUE, only check the regulation of the gene in the second column upon the gene in the first column. [Optional]

	--sample=SAMPLE
		Random sampling how many lines from the provided gene file (indicated with --genes). Set to -1 if no sampling is needed [Optional]

	--genes=GENES
		Genes to be subsetted, please also set --sample to a positive integer. [Optional]

	--core=CORE
		Cpu cores used to run this script concurrently. [Optional]

	--seed=SEED
		Seed used in set.seed(). Set to -1 to use rownumbers of PPI file as seeds. [Optional]

	--npermut=NPERMUT
		Number of permutations for fgsea. [Optional]

	--nproc=NPROC
		Number of processes to be used for fgsea. [Optional]

	--min_size=MIN_SIZE
		Min cell line number used in fgsea. [Optional]

	--max_size=MAX_SIZE
		Max cell line number used in fgsea. [Optional]

	--precise_mediation_pvalue=PRECISE_MEDIATION_PVALUE
		Whether more precise p-values are needed (TRUE: 100000 permutations; FALSE: 1000 permutations). [Optional]

	--remove_cis_confounding=REMOVE_CIS_CONFOUNDING
		Whether to remove trans records with potential cis regulation. [Optional]

	--cna_adjust=CNA_ADJUST
		Whether remove CNA confounding effect. Specify --cna also when using this argument [Optional]

	--cna_discretization=CNA_DISCRETIZATION
		Whether discretize the CNA values. Specify --cna and --cna_adjust also when using this argument [Optional]

	--cna_permutation_size=CNA_PERMUTATION_SIZE
		Calculate empirical p-values based on the predefined permutation size. Specify --cna and --cna_adjust also when using this argument [Optional]

	--cna_adjust_Pvalue_cutoff=CNA_ADJUST_PVALUE_CUTOFF
		Remove CNA confounding effect according to the P-value cutoff. Specify --cna and --cna_adjust also when using this argument [Optional]

	--cna_amplification_cutoff=CNA_AMPLIFICATION_CUTOFF
		CN > this cutoff is considered as amplication. Specify --cna and --cna_adjust also when using this argument [Optional]

	--cna_deletion_cutoff=CNA_DELETION_CUTOFF
		CN < this cutoff is considered as deletion. Specify --cna and --cna_adjust also when using this argument [Optional]

	--remove_contradiction=REMOVE_CONTRADICTION
		Whether remove contradiction between CSEA and mediation analysis. [Optional]

	--FDR_cutoff_cis=FDR_CUTOFF_CIS
		The FDR cutoff for defining cis genetic alterations. [Optional]

	--FDR_cutoff_trans=FDR_CUTOFF_TRANS
		The FDR cutoff for defining trans genetic alterations. [Optional]

	--CSEA_cis_cutoff=CSEA_CIS_CUTOFF
		The CSEA p-value cutoff for removing insignificant genetic alterations before p-value combination. [Optional]

	--mediation_cis_cutoff=MEDIATION_CIS_CUTOFF
		The mediation p-value cutoff for removing insignificant genetic alterations before p-value combination. [Optional]

	--CSEA_trans_cutoff=CSEA_TRANS_CUTOFF
		The CSEA p-value cutoff for removing insignificant genetic alterations before p-value combination. [Optional]

	--mediation_trans_cutoff=MEDIATION_TRANS_CUTOFF
		The mediation p-value cutoff for removing insignificant genetic alterations before p-value combination. [Optional]

	--output_mutation_site_centric_table=OUTPUT_MUTATION_SITE_CENTRIC_TABLE
		Output mutation information tables [Optional]

	-h, --help
		Show this help message and exit

Examples:

Download and unzip necessary data from https://figshare.com/projects/DICDAT_A_Data-driven_Method_Uncovers_Infrequent_Cancer-driving_Alterations_and_Associated_Targets_from_Multi-sample_Omics_Data/131000 and put into the right folder. 
You can also run data_prepare.Rmd if you are interested in processing raw data by yourself, however, it may require a lot of RAMs.


1. Run DepMap CCLE data (Line 1 to 20)
Rscript DICDAT.R --cl=TRUE --core=1 --folder=ex --name=ex_  --start_PPI=1 --end_PPI=20 --precise_mediation_pvalue=TRUE --cna_adjust=TRUE --ppi=./data/PPI_filtered.txt --exp=./data/expr_processed.tsv --ga=./data/mutation_processed.tsv --dep=./data/dependency_processed.tsv --cna=./data/CNV_processed.tsv --remove_cis_confounding=FALSE --FDR_cutoff_cis=0.1 --FDR_cutoff_trans=0.1 --cna_adjust=TRUE --seed=1 --folder=ex --output_mutation_site_centric_table=FALSE --precise_mediation_pvalue=TRUE --remove_contradiction=TRUE

2. Run TCGA data (Line 1 to 20)
Rscript DICDAT.R --cl=FALSE --core=1 --folder=ex --name=ex_  --start_PPI=1 --end_PPI=20 --precise_mediation_pvalue=TRUE --cna_adjust=TRUE --ppi=./data/PPI_filtered.txt --exp=./TCGA/TCGA_expr_processed.tsv --ga=./TCGA/TCGA_mutation_processed.tsv --dep=./data/dependency_processed.tsv --cna=./TCGA/TCGA_CNV_processed.tsv --map=./TCGA/TCGA_sample_mapping.tsv --remove_cis_confounding=FALSE --FDR_cutoff_cis=0.1 --FDR_cutoff_trans=0.1 --cna_adjust=TRUE --seed=1 --folder=ex --output_mutation_site_centric_table=FALSE --precise_mediation_pvalue=TRUE --remove_contradiction=TRUE
