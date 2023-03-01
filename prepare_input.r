#!/usr/bin/R

## Prepare the input files with ROIs for DNA FISH, RNA FISH or both

require(data.table)
require(ggplot2)

## Read ROI

roi_data = rbindlist(list(
    fread("data/rois/df_DNA_FISH.txt", col.names=c("oligos","chrom", "DNA_start", "DNA_end", "type", "ref", "length")),
    fread("data/rois/df_RNA_FISH.txt")[, .(
        chrom=gtf.chr, Gene_start=gtf.start, Gene_end=gtf.end, Gene_strand=gtf.strand, Gene_name=gtf.name, Gene_id=ids, ref=gtf.ref, length=length, oligos=oligos)],
    fread("data/rois/df_DNA_RNA_FISH.txt")[, .(
        chrom=DNA.chr, DNA_start=DNA.start, DNA_end=DNA.end, window=DNA.windowID, type=DNA.group,
        Gene_start=RNA.start, Gene_end=RNA.end, Gene_strand=RNA.strand, Gene_name=RNA.name, Gene_id=RNA.ensembleID, length=length)]
), use.names=T, fill=T)
#head(roi_data)

roi_data[, design_type := "DNA+RNA"]
roi_data[is.na(Gene_start), design_type := "DNA"]
roi_data[is.na(DNA_start), design_type := "RNA"]


## Plot ROIs for DNA+RNA FISH

#options(repr.plot.width=24, repr.plot.height=12)
#ggplot(roi_data[design_type == "DNA+RNA"]) +
#    geom_segment(aes(x=DNA_start/2e6, xend=DNA_end/2e6, color="DNA"), y=1, yend=1) +
#    geom_segment(aes(x=Gene_start/2e6, xend=Gene_end/2e6, color="Gene"), y=2, yend=2) +
#    facet_wrap(~sprintf("%s - %s", Gene_name, Gene_id), scales="free") + ylim(0, 3) + scale_color_brewer(palette="Set1") +
#    theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="bottom") +
#   labs(x="Genomic coordinate (Mb)", color="Target", title="mm10", subtitle="Chromosome 5 regions of interest")


## HUSH status

hush_threshold = 8e5

# DNA FISH only
unhushed = roi_data[design_type == "DNA" & DNA_start%%1e6 > hush_threshold]
unhushed
dim(roi_data[design_type == "DNA"])
round(nrow(unhushed)/nrow(roi_data[design_type == "DNA"])*100, 2)

# RNA FISH only
unhushed = roi_data[design_type == "RNA" & Gene_start%%1e6 > hush_threshold]
unhushed
dim(roi_data[design_type == "RNA"])
round(nrow(unhushed)/nrow(roi_data[design_type == "RNA"])*100, 2)

# DNA+RNA FISH
unhushed = roi_data[design_type == "DNA+RNA" & DNA_start%%1e6 > hush_threshold]
unhushed
dim(roi_data[design_type == "DNA+RNA"])
round(nrow(unhushed)/nrow(roi_data[design_type == "DNA+RNA"])*100, 2)

# Reduce ROI
roi_data[design_type != "RNA", Window_start := min(DNA_start, Gene_start, na.rm=T), by=c("DNA_start", "DNA_end")]
roi_data[design_type != "RNA", Window_end := max(DNA_end, Gene_end, na.rm=T), by=c("DNA_start", "DNA_end")]
roi_data[design_type == "RNA", Window_start := min(DNA_start, Gene_start, na.rm=T), by=c("Gene_start", "Gene_end")]
roi_data[design_type == "RNA", Window_end := max(DNA_end, Gene_end, na.rm=T), by=c("Gene_start", "Gene_end")]
head(roi_data)

windows = unique(roi_data[, .(Window_start, Window_end)])
windows[, window_id := .I]
setkeyv(windows, c("Window_start", "Window_end"))

setkeyv(roi_data, c("Window_start", "Window_end"))
roi_data2 = windows[roi_data]

# Manual corrections
# remove DNA only FISH regions that are already present in the DNA+RNA FISH (done manually on input)
roi_data2[, .N, by=window_id][N > 1]
roi_data2[!is.na(Gene_name), .N, by=Gene_name][N>1]
roi_data2[, max(window_id)]

setorder(roi_data2,window_id)#,chrom,Window_start,Window_end)

# Write to file
fwrite(roi_data2, "data/rois/all_regions.tsv", sep="\t")