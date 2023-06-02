#!/usr/bin/env Rscript

# this script extracts annotations for partial/incomplete ORFs from prodigal results

# configure environment
library(tidyverse)
library(Biostrings)
library(ape)
library(optparse)

# options
option_list = list(
  make_option("--full_gff", type="character", default="", help="full gff provided by prokka [default= %default]", dest="full_gff"),
  make_option("--prodigal_gff", type="character", default="", help="provided by prodigal - provides partial gene calls [default= %default]", dest="prodigal_gff"),
  make_option("--assembly", type="character", default="", help="fna file from prokka - contains assembled contigs [default= %default]", dest="assembly")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load assembly data
assembly <- readDNAStringSet(opt$assembly)

# setting up output files
output_fasta <- paste0(gsub(".gff$", "", opt$full_gff), "_partials.fa")
output_gff_table <- paste0(gsub(".gff$", "", opt$full_gff), "_partials.gff_table")

# get the fasta portion of the gff
fasta_portion <- system(paste0("grep -n \"^>\" ", opt$full_gff), intern = TRUE)

# remove the sequence
system(paste0("head -n ", as.numeric(gsub(":.+", "", fasta_portion[1]))-2, " ", opt$full_gff, " > ", paste0(gsub(".gff$", "", opt$full_gff), "_noseq.gff")))

# read the reduced gff
gff <- read.gff(paste0(gsub(".gff$", "", opt$full_gff), "_noseq.gff")) %>%
  filter(type == "CDS") %>%
  mutate(locus_number = gsub(";.+", "", gsub("^.+?_", "", attributes)))

# get the partial entries
gff_partial <- read.gff(opt$prodigal_gff) %>%
  select(seqid, type, start, end, strand, phase, attributes) %>%
  mutate(partial = gsub(";.*", "", gsub(".*;partial=", "", attributes)))

gff_partial <- gff_partial %>%
  filter(partial %in% c("10", "01"))

system(paste0("rm ", output_fasta))

# process the partial entries
if (nrow(gff_partial) == 0) {
  write.table("", output_fasta)
} else {
  max_locus_number <- max(as.numeric(gff$locus_number))
  locus_tag <- gsub("^ID=", "", gsub("_.+", "", gff[1,]$attributes))
  for (i in 1:nrow(gff_partial)) {
    locus_number <- max_locus_number + i
    locus_number_format <- formatC(locus_number, width = max(nchar(gff$locus_number)), format = "d", flag = "0")
    rec <- subseq(assembly[gff_partial[i,]$seqid], start = gff_partial[i,]$start, gff_partial[i,]$end)
    names(rec) <- paste(names(rec), paste(locus_tag, locus_number_format, sep = "_"), sep = "|")
    if (gff_partial[i,]$strand == "-") {
      rec <- reverseComplement(rec)
    }
    writeXStringSet(rec, output_fasta, append = TRUE)
    gff <- gff %>%
      bind_rows(data.frame(seqid = gff_partial[i,]$seqid,
                           source = "manual_addition",
                           type = "CDS",
                           start = gff_partial[i,]$start,
                           end = gff_partial[i,]$end,
                           score = NA,
                           strand = gff_partial[i,]$strand,
                           phase = gff_partial[i,]$phase,
                           attributes = paste("ID", paste(locus_tag, locus_number_format, sep = "_"), sep = "="),
                           locus_number = locus_number_format))
    
  }  
}

# write to file
write.table(gff, file = output_gff_table, row.names = FALSE, sep = "\t")

