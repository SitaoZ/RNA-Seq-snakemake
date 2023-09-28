library(getopt)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

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

go_analysis <- enrichGO(df$GeneID, OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pvalueCutoff = 0.05)
head(go_analysis)
# barplot(go_analysis, drop=TRUE, showCategory=12)
dotplot(go_analysis, showCategory = 10, font.size=14)
ggsave(opt$output, width = 8, height = 6, units = 'in', dpi=300)
