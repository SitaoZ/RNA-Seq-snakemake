library(getopt)
library(ggplot2)

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

ggplot(data = df, aes(x = log2FoldChange.Treat.Ctrl., y = -log10(Padj),
                      col=Up.Down.Regulation)) +
  geom_point(size=0.8) +
  theme_classic() +
  scale_color_manual(values=c(Up = "#F8766D", Down = "#00BA38"))
ggsave(opt$output, width = 8, height = 6,units = 'in',dpi=300)
