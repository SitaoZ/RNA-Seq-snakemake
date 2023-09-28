library(pheatmap)
library(ggplot2)
library(getopt)
library(dplyr)

spec <- matrix(
  c("input",  "i", 2, "character", "This is input file",
    "diff", "d", 2, "character", "This is diff gene file",
    "output", "o", 2, "character",  "This is output file",
    "help",   "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
if( !is.null(opt$help) || is.null(opt$input) || is.null(opt$output) || is.null(opt$diff)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

exp <- read.csv(opt$input, sep="\t")
diff <- read.csv(opt$diff, sep="\t")

DGE <- filter(exp, exp$GeneID %in% diff$GeneID)
DGE <- DGE[apply(DGE!=0, 1, all),]
rownames(DGE) <- DGE$GeneID
DGE$GeneID <- NULL
DGE_log <- log(DGE)
pplot <- pheatmap(DGE_log, scale='row')
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(pplot, opt$output)

