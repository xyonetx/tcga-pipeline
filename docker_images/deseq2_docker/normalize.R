suppressMessages(suppressWarnings(library('DESeq2', character.only=T, warn.conflicts = F, quietly = T)))

# args from command line:
args<-commandArgs(TRUE)

# path to the count matrix
COUNTS_FILE <- args[1]

# path to the output (normalized) expression
OUTPUT_FILE_PATH <- args[2]

# Don't mangle the column names quite yet. Keep them as-is and then we can convert them after we
# create a column name mapping
counts <- read.table(COUNTS_FILE, sep='\t', stringsAsFactors=F, row.names=1, header=T, check.names=F)
orig_colnames = colnames(counts)

# cast the column names to "proper" names for R. This way we don't run into trouble with
# any called functions
new_colnames = make.names(orig_colnames)
colnames(counts) = new_colnames

# Create a map of the column names so we can re-assign them at the end
colname_mapping = setNames(orig_colnames, new_colnames)

# Now run the normalization:
deseq_sf = estimateSizeFactorsForMatrix(counts)
norm_mtx = sweep(counts,2,deseq_sf,'/')

# map back to the original column/sample names
colnames(norm_mtx)<-colname_mapping[colnames(norm_mtx)]

write.table(norm_mtx, file=OUTPUT_FILE_PATH, sep='\t', quote=FALSE)
