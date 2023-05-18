suppressMessages(suppressWarnings(
    library("DESeq2", character.only=T, warn.conflicts = F, quietly = T)))

###########################################################################

# args from command line:
args<-commandArgs(TRUE)

# the raw counts for the TCGA cohort, in TSV format
RAW_COUNT_MATRIX <- args[1]

# the full annotation matrix:
FULL_ANN_MATRIX <- args[2]

# the gene we are using to segment the samples
SELECTED_GENE <- args[3]

# the quantile for the lower and upper cutoffs
LOW_Q <- as.numeric(args[4])
HIGH_Q <- as.numeric(args[5])

# An annotation file giving the samples that are high or low
# on the distribution
ANN_OUTPUT <- args[6]

###########################################################################

# Don't mangle the column names quite yet. Keep them as-is and then we can convert them after we
# create a column name mapping
counts <- read.table(RAW_COUNT_MATRIX, 
                    sep='\t', 
                    stringsAsFactors=F, 
                    row.names=1, 
                    header=T, 
                    check.names=F)
orig_colnames = colnames(counts)

# cast the column names to "proper" names for R. This way we don't run into trouble with
# any called functions
new_colnames = make.names(orig_colnames)
colnames(counts) = new_colnames

# Create a map of the column names so we can re-assign them at the end
colname_mapping = setNames(orig_colnames, new_colnames)

# create a "dummy" annotation matrix which we will later fill in
# with high or low expression status
ann_df <- data.frame(expression_state=rep(NA, length(new_colnames)))
rownames(ann_df) <- new_colnames

# Now run the normalization so we can find the "extremes"
# for the gene of interest:
deseq_sf = estimateSizeFactorsForMatrix(counts)
norm_mtx = sweep(counts,2,deseq_sf,'/')

# extract the expressions for the gene of interest and find
# the cutoffs for the desired quantiles
gene_array <- as.matrix(norm_mtx[SELECTED_GENE,])
cutoffs <- unname(
                  quantile(
                           gene_array, 
                           na.rm=T, 
                           probs=c(LOW_Q, HIGH_Q)
                  )
           )
low_cutoff <- cutoffs[1]
high_cutoff <- cutoffs[2]

# now get the samples corresponding to the low and high-expression groups
low_samples <- names(gene_array[,gene_array < low_cutoff])
high_samples <- names(gene_array[,gene_array > high_cutoff])

# complete the annotation dataframe so we can have for later
# in case we need it:
ann_df[low_samples, 'expression_state'] <- 'low'
ann_df[high_samples, 'expression_state'] <- 'high'

# the "un-interesting" samples (those with NAs) can be dropped
ann_df <- na.omit(ann_df)

# now load the full annotation matrix to get the other metadata for survival:
keep_cols =  c('days_to_death','vital_status')
full_annotations = read.csv(FULL_ANN_MATRIX, row.names=1)
full_annotations = full_annotations[,keep_cols]
rownames(full_annotations) = make.names(rownames(full_annotations))

# merge with the dataframe describing the expression status (high or low)
joined_ann = merge(ann_df, full_annotations, by=0)
rownames(joined_ann) <- colname_mapping[joined_ann[,'Row.names']]
joined_ann <- subset(joined_ann, select=c('days_to_death','vital_status','expression_state'))

write.table(joined_ann, ANN_OUTPUT, sep='\t', quote=F)