library(getopt)
library(ggplot2)
library(clusterProfiler)
library(R.utils)

# R.utils::setOption("clusterProfiler.download.method","auto")
R.utils::setOption("clusterProfiler.download.method","wget")

spec <- matrix(
  c("input",  "i", 2, "character", "This is input file",
    "output", "o", 2, "character",  "This is output file",
    "help",   "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
if( !is.null(opt$help) || is.null(opt$input) || is.null(opt$output) ){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}


df <- read.csv(opt$input, sep="\t")
head(df$GeneID)
#gene <- bitr(df$GeneID,fromType = 'TAIR',toType = 'ENTREZID',OrgDb = "org.At.tair.db")
#dim(gene)
#head(gene)
#keytypes(org.At.tair.db)
#search_kegg_organism('ath', by='kegg_code')
kegg_analysis <- enrichKEGG(df$GeneID, organism="ath", pvalueCutoff = 1, pAdjustMethod = 'BH')
head(kegg_analysis)
#dotplot(kegg_analysis, title="Enrichment KEGG_dot")
ggsave(opt$output, width = 8, height = 6, units='in', dpi=300)
