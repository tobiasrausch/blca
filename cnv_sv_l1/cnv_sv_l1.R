library(ggplot2)
library(reshape2)
library(ggsignif)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=axisTtlFontSize))

# Read data
x=read.table("summary.table", header=T)
x = x[x$l1sv != "None",]
x = x[x$l1loose != "None",]
x$l1sv = as.numeric(x$l1sv)
x$l1loose = as.numeric(x$l1loose)
x = x[(x$l1sv > 0) | (x$l1loose > 0),]

# Filter by histology (min. 5 samples)
min_samples = 5
histology_counts = table(x$histology)
histologies_to_keep = names(histology_counts[histology_counts >= min_samples])
x_filtered = x[x$histology %in% histologies_to_keep, ]

## Plotting
df=melt(x_filtered[,c("sample","histology","l1sv","l1loose")], id.vars=c("sample","histology"))

png("cnv_hist.png", width=1200, height=600)
p = ggplot(data=df, aes(x=variable, y=value))
p = p + geom_boxplot(aes(group=variable))
p = p + facet_wrap(~histology, scales="free_y")
p = p + geom_signif(comparisons = list(c("l1sv", "l1loose")), vjust=1.5, 
                   test.args = list(alternative = "two.sided", paired = TRUE, exact=FALSE))
p = p + scienceTheme
p = p + ylab("Fraction of somatic L1s\nat CNV breakpoints")
p = p + xlab("CNV with (without) SV call")
p = p + scale_x_discrete(labels = c("CNV with SV", "CNV without SV"))
p
dev.off()

png("cnv.png", width=800, height=600)
p = ggplot(data=df, aes(x=variable, y=value))
p = p + geom_boxplot(aes(group=variable), alpha=0.5)
p = p + geom_line(aes(group=sample), alpha=0.2, color="gray50")
p = p + geom_point(alpha=0.4, size=2)
p = p + geom_signif(comparisons = list(c("l1sv", "l1loose")), map_signif_level = TRUE, textsize = 6, test.args = list(alternative = "two.sided", paired = TRUE, exact=FALSE))
p = p + scienceTheme
p = p + ylab("Fraction of somatic L1s\nat CNV breakpoints")
p = p + xlab("CNV with (without) SV call")
p = p + scale_x_discrete(labels = c("CNV with SV", "CNV without SV"))
p
dev.off()

# Save filtered data for inspection
write.table(x_filtered, "summary.table.filtered", sep="\t", row.names=F, quote=F)
