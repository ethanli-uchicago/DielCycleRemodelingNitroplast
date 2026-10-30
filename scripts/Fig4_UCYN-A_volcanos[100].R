# Reproduce the main-text and supplementary UCYN-A proteomics volcano figures.
#
# INPUTS: use locally downloaded, unmodified files from the Coale et al. (2024)
# Dryad deposit https://doi.org/10.5061/dryad.2z34tmptf (version 301306,
# updated 24 June 2024). 
# Put these THREE files together in a folder named input_data (or choose another
# folder with the first command-line argument / input_dir setting below):
#   ADK1075_ProteinQuantifications.csv   -- original protein quantities
#   ADK1075_proteomics_DB_2.fasta        -- sequences for theoretical peptide counts
#   UCYNA2_genome_descriptions.csv      -- genome annotations for candidate selection
#
# USAGE (relative paths are resolved against your working directory):
#   Rscript Fig4_UCYN-A_volcanos.R
#   Rscript Fig4_UCYN-A_volcanos.R input_data results kofam_reference
# In RStudio, edit input_dir/output_dir below, set your working directory, then
# source this script.
# Also download ko_list.gz and profiles.tar.gz from:
# https://www.genome.jp/ftp/db/kofam/ into kofam_reference.
# This pinned version uses the reference retrieved 21 October 2026.
# Install HMMER 3.4 (hmmsearch on PATH, or in WSL on Windows).
# The R script extracts references, annotates Dryad sequences, and makes both plots.
# Allow several GB of disk space and time for the full-profile sequence search.
#
# DEPENDENCIES (install once if needed; not installed automatically):
# install.packages(c('dplyr','tidyr','readr','tibble','ggplot2','ggrepel',
#                    'scales','statmod','digest','openxlsx','BiocManager'))
# BiocManager::install('limma', ask=FALSE, update=FALSE)
#
# OUTPUTS: Main_text_Fig4g.pdf/.png and SI_all_UCYNA.pdf/.png in output_dir;
# Supplementary_Data_S1.xlsx and Supplementary_Data_S2.xlsx, plus
# matching plain CSV tables in supplementary_data/S1 and supplementary_data/S2;
# complete source tables, testability audits and checks in output_dir/source_data;
# sequence-based KO evidence, reference hashes and raw searches in KO_annotation.
# Existing output files with these names will be replaced on rerun.
#


input_dir <- 'input_data'
output_dir <- 'volcano_results'
reference_dir <- 'kofam_reference'
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>=1) input_dir <- args[1]
if(length(args)>=2) output_dir <- args[2]
if(length(args)>=3) reference_dir <- args[3]
if(length(args)>3) stop('Usage: Rscript Fig4_UCYN-A_volcanos.R [input-folder] [output-folder] [KOfam-reference-folder]')
needed <- c('dplyr','tidyr','readr','tibble','ggplot2','ggrepel','scales','limma','statmod','digest','openxlsx')
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Missing R packages: ',paste(missing,collapse=', '),'. See installation comments above.')
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(tibble);library(ggplot2)})
required_files <- c('ADK1075_ProteinQuantifications.csv','ADK1075_proteomics_DB_2.fasta','UCYNA2_genome_descriptions.csv')
absent <- required_files[!file.exists(file.path(input_dir,required_files))]
if(length(absent)) stop('Missing Dryad input files: ',paste(absent,collapse=', '))
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
outdir <- file.path(output_dir,'source_data')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

# 1. Read the original quantities; fit one statistical model; calculate abundance
#
quant_file <- file.path(input_dir, "ADK1075_ProteinQuantifications.csv")
fasta_file <- file.path(input_dir, "ADK1075_proteomics_DB_2.fasta")
stopifnot(digest::digest(file = quant_file, algo = "sha256") == "825acccd992863a6b0df14b0633f1d3f4c66cccb2fa4a23bfe179ab10818b3c1")
stopifnot(digest::digest(file = fasta_file, algo = "sha256") == "76c51de9db4d8e19fc36d797e70227eb78cdc4d5ad5e8da10a79dfccb04c985a")

# Dryad has three header rows. Derive sample identity from its first two rows
headers <- read.csv(quant_file, header = FALSE, nrows = 3, stringsAsFactors = FALSE, check.names = FALSE)
raw <- read.csv(quant_file, skip = 3, header = FALSE, stringsAsFactors = FALSE,
                na.strings = c("NaN", "NA", ""), check.names = FALSE)
stopifnot(ncol(raw)==24, headers[1,24]=="Sequence source", !anyDuplicated(raw[[1]]))
quantity_idx <- which(grepl("PG.Quantity", as.character(headers[3,])))
stopifnot(length(quantity_idx)==12)
sample_key <- tibble(column_index = quantity_idx,
  published_column = as.character(headers[3,quantity_idx]),
  condition = as.character(headers[1,quantity_idx]),
  sample_type = as.character(headers[2,quantity_idx])) %>%
  group_by(condition, sample_type) %>% mutate(replicate = row_number()) %>% ungroup() %>%
  mutate(sample = paste0(if_else(sample_type=="Isolated UCYNA", "Diazoplast", "Whole"),
                        replicate + if_else(condition=="Night", 3L, 0L)),
         used = sample_type=="Isolated UCYNA")
stopifnot(identical(sample_key$sample, c("Diazoplast1","Diazoplast2","Diazoplast3","Whole1","Whole2","Whole3",
                                     "Diazoplast4","Diazoplast5","Diazoplast6","Whole4","Whole5","Whole6")))
write_csv(sample_key, file.path(outdir,"Sample_key.csv"))
quant <- tibble(PG.ProteinAccessions=raw[[1]], sequence_source=raw[[24]],
                Origin=if_else(raw[[24]]=="UCYN-A genome", "UCYNA", raw[[24]]))
for (i in seq_len(nrow(sample_key))) quant[[sample_key$sample[i]]] <- as.numeric(raw[[sample_key$column_index[i]]])
write_csv(quant, file.path(outdir,"Parsed_Dryad_quantities.csv"))
day_samples <- sample_key$sample[sample_key$used & sample_key$condition=="Day"]
night_samples <- sample_key$sample[sample_key$used & sample_key$condition=="Night"]
stopifnot(length(day_samples)==3, length(night_samples)==3)

# Refit the statistical model 

mat <- as.matrix(quant[,c(day_samples,night_samples)])
rownames(mat) <- quant$PG.ProteinAccessions
# All observed quantities in the pinned Dryad samples are strictly positive.
# Log-transform directly; retain missing values without imputation.
# Stop on unexpected zeros.
stopifnot(!any(is.infinite(mat)), !any(mat<=0,na.rm=TRUE))
log2int <- log2(mat)
sample_medians <- apply(log2int,2,median,na.rm=TRUE)
norm <- sweep(log2int,2,sample_medians,"-")
keep <- rowSums(!is.na(norm[,day_samples,drop=FALSE]))>=2 & rowSums(!is.na(norm[,night_samples,drop=FALSE]))>=2
condition <- factor(c(rep("Day",3),rep("Night",3)),levels=c("Night","Day"))
design <- model.matrix(~0+condition); colnames(design) <- levels(condition)
fit <- limma::lmFit(norm[keep,,drop=FALSE],design)
fit <- limma::eBayes(limma::contrasts.fit(fit,limma::makeContrasts(DayVsNight=Day-Night,levels=design)),robust=TRUE,trend=TRUE)
fc <- as.numeric(fit$coefficients[,1]); se <- as.numeric(fit$stdev.unscaled[,1]*sqrt(fit$s2.post))
critical <- qt(.975,df=fit$df.total)
limma <- tibble(PG.ProteinAccessions=rownames(fit$coefficients), log2FC_DayMinusNight=fc,
 SE=se, CI_lo=fc-critical*se,CI_hi=fc+critical*se,p_value=as.numeric(fit$p.value[,1]),
 FDR=p.adjust(p_value,"BH"),overall_mean=rowMeans(norm[keep,,drop=FALSE],na.rm=TRUE)) %>%
 left_join(quant %>% select(PG.ProteinAccessions,Origin,sequence_source),by="PG.ProteinAccessions")
write_csv(limma,file.path(outdir,"DayVsNight_limma_results_INTENSITY.csv"))
write_csv(tibble(sample=colnames(mat),log2_median=sample_medians,pseudocount=0),file.path(outdir,"Normalization_parameters.csv"))
write_csv(tibble(PG.ProteinAccessions=rownames(mat),n_day=rowSums(!is.na(mat[,day_samples])),
 n_night=rowSums(!is.na(mat[,night_samples])),tested=keep),file.path(outdir,"Statistical_filter_audit.csv"))

tryptic_peptide_count <- function(sequence, min_length = 6L, max_length = 30L) {
  sequence <- gsub("[^A-Z]", "", toupper(sequence))
  aa <- strsplit(sequence, "", fixed = TRUE)[[1]]
  if (!length(aa)) return(NA_integer_)
  cleavage <- which(aa %in% c("K", "R") & c(aa[-1] != "P", TRUE))
  ends <- sort(unique(c(cleavage, length(aa))))
  starts <- c(1L, head(ends, -1L) + 1L)
  peptides <- mapply(substr, sequence, starts, ends, USE.NAMES = FALSE)
  sum(nchar(peptides) >= min_length & nchar(peptides) <= max_length)
}

read_fasta <- function(path) {
  x <- readLines(path, warn = FALSE)
  header_idx <- which(startsWith(x, ">"))
  header_end <- c(header_idx[-1] - 1L, length(x))
  ids <- sub("^>([^[:space:]]+).*", "\\1", x[header_idx])
  sequences <- vapply(seq_along(header_idx), function(i) {
    paste0(x[(header_idx[i] + 1L):header_end[i]], collapse = "")
  }, character(1))
  tibble(accession = ids, sequence = sequences)
}


fasta <- read_fasta(fasta_file) %>% mutate(observable_tryptic_peptides=vapply(sequence,tryptic_peptide_count,integer(1)))
stopifnot(!anyDuplicated(fasta$accession))
quant_long <- quant %>%
  filter(Origin == "UCYNA") %>%
  select(PG.ProteinAccessions, all_of(c(day_samples, night_samples))) %>%
  inner_join(fasta %>% select(accession, observable_tryptic_peptides),
             by = c("PG.ProteinAccessions" = "accession")) %>%
  filter(is.finite(observable_tryptic_peptides), observable_tryptic_peptides > 0) %>%
  pivot_longer(all_of(c(day_samples, night_samples)), names_to = "sample", values_to = "PG_Quantity") %>%
  filter(is.finite(PG_Quantity), PG_Quantity > 0) %>%
  mutate(
    iBAQ_like = PG_Quantity / observable_tryptic_peptides,
    condition = if_else(sample %in% day_samples, "Day", "Night"),
    replicate = if_else(condition == "Day", match(sample, day_samples), match(sample, night_samples)),
    condition = factor(condition, levels = c("Night", "Day"))
  ) %>%
  group_by(sample) %>%
  mutate(riBAQ_percent_UCYNA = 100 * iBAQ_like / sum(iBAQ_like, na.rm = TRUE)) %>%
  ungroup()

# Annotation rules are fixed independently of abundance and differential statistics.
annotation_file <- file.path(input_dir,'UCYNA2_genome_descriptions.csv')
stopifnot(digest::digest(file=annotation_file,algo='sha256')=='3ccd99474659ccc7323043a42ea7331689f7406f86e9af24630e2cc059a5cda7')
ann <- read_csv(annotation_file,show_col_types=FALSE) %>% rename(accession=Accession,description=Description) %>%
 mutate(description=trimws(sub(' \\[.*$', '', description)))
stopifnot(nrow(ann)==1186,!anyDuplicated(ann$accession))
expected_chunk_hashes <- c('7c45b48ebf7f0f6489bc76b204507c32a650204593755e0e9bfb0552eda8ecfa','f312e337be6c3e6b1b2d7a279c61404eb04657148a8b85f944842aad17723dbf','bb64a5361753e9b198caeb6cfb4f1f7013576aa0ba7805139e0547564a69a0bd','df38373982ef97467a656e707504bc7c6a91edd9be818b519d458692b0f8ba30','d1ce8233a917406a7611833d70365464afd6e7ac458b83bcb4d75c1149357604','0070971197da7a802c148e81a6b7a2d949466cef2ffca8a4442dcfb5787f6ee7','901a9299b51021bc23c529eb57c39218a2f59c060b0750882621092effd8d6c2','6f775620b558a65414d3c387b3e73b6537568b334eb7370a1602af6c5d10f431')
# Sequence-based KO annotation. The only study sequences are from Dryad.
# KOfam is a public reference database
ko_dir <- file.path(output_dir,'KO_annotation')
dir.create(ko_dir,recursive=TRUE,showWarnings=FALSE)
sha <- function(p) digest::digest(file=p,algo='sha256')
ko_list_file <- file.path(reference_dir,'ko_list.gz')
archive <- file.path(reference_dir,'profiles.tar.gz')
if(!all(file.exists(c(ko_list_file,archive)))) stop('Put ko_list.gz and profiles.tar.gz from https://www.genome.jp/ftp/db/kofam/ in reference_dir (third argument).')
stopifnot(sha(archive)=='6f5e69cd2c7811ad0eb6f4dc5cec603f4c6c7034c4321ba21ac84b68a7825a59')
stopifnot(sha(ko_list_file)=='4131a8aebc7e68f9a03939bdd575d3cee58c2b42135b1893828270bb179fca01')
con <- gzfile(ko_list_file,'rt')
ko_ref <- read.delim(con,check.names=FALSE,quote='',stringsAsFactors=FALSE); close(con)
ko_ref$threshold <- suppressWarnings(as.numeric(ko_ref$threshold))
query <- fasta %>% semi_join(ann,by='accession') %>% arrange(accession)
stopifnot(nrow(query)==1186)
query_file <- file.path(ko_dir,'UCYNA_from_Dryad.fasta')
writeLines(as.vector(rbind(paste0('>',query$accession),query$sequence)),query_file,useBytes=TRUE)
# Normalize line endings to LF, including on Windows
con <- file(query_file,'wb'); writeBin(charToRaw(paste0(paste(as.vector(rbind(paste0('>',query$accession),query$sequence)),collapse='\n'),'\n')),con); close(con)
chunk_dir <- file.path(reference_dir,'chunks')
chunks <- file.path(chunk_dir,sprintf('chunk_%02d.hmm',0:7))

if(!all(file.exists(chunks))) {
  members <- utils::untar(archive,list=TRUE)
  members <- members[grepl('\\.hmm$',members)]
  if(any(grepl('(^/|(^|/)\\.\\.(/|$)|^[A-Za-z]:)',members))) stop('Unsafe reference archive paths')
  stopifnot(length(members)==27930)
  extracted <- file.path(reference_dir,'extracted')
  dir.create(extracted,recursive=TRUE,showWarnings=FALSE)
  if(utils::untar(archive,files=members,exdir=extracted)!=0) stop('Reference extraction failed')
  dir.create(chunk_dir,recursive=TRUE,showWarnings=FALSE)
  handles <- lapply(chunks,file,open='wb')
  for(i in seq_along(members)) {
    p <- file.path(extracted,members[i])
    writeBin(readBin(p,'raw',n=file.info(p)$size),handles[[(i-1)%%8+1]])
  }
  invisible(lapply(handles,close))
}
chunk_hashes <- vapply(chunks,sha,character(1))
stopifnot(identical(unname(chunk_hashes),expected_chunk_hashes))
# Native HMMER on PATH is preferred. Windows can use HMMER installed in WSL.
use_wsl <- !nzchar(Sys.which('hmmsearch')) && .Platform$OS.type=='windows' && nzchar(Sys.which('wsl.exe'))
if(!use_wsl && !nzchar(Sys.which('hmmsearch'))) stop('Install HMMER 3.4 and place hmmsearch on PATH (or install it in WSL on Windows).')
native_path <- function(p) {
  p <- normalizePath(p,winslash='/',mustWork=TRUE)
  if(use_wsl) trimws(system2('wsl.exe',c('-e','wslpath','-a',shQuote(p)),stdout=TRUE)) else p
}
hmmer_help <- if(use_wsl) system2('wsl.exe',c('-e','hmmsearch','-h'),stdout=TRUE) else system2('hmmsearch','-h',stdout=TRUE)
hmmer_version <- grep('HMMER ',hmmer_help,value=TRUE)[1]
if(!grepl('HMMER 3.4 ',hmmer_version,fixed=TRUE)) stop('This pinned workflow requires HMMER 3.4.')
query_hash <- sha(query_file)
raw_dir <- file.path(ko_dir,'hmmer'); dir.create(raw_dir,showWarnings=FALSE)
tables <- file.path(raw_dir,sprintf('chunk_%02d.tbl',0:7))
for(i in seq_along(chunks)) {
  key <- paste(query_hash,chunk_hashes[i],hmmer_version,'--cpu 1 -E 0.01 --domE 0.01',sep='|')
  stamp <- paste0(tables[i],'.cache_key')
  valid <- file.exists(stamp) && file.exists(tables[i]) && identical(readLines(stamp,warn=FALSE),c(key,sha(tables[i]))) &&
    any(grepl('^# \\[ok\\]',tail(readLines(tables[i],warn=FALSE),10)))
  if(!valid) {
    message('Searching KOfam chunk ',i,' / ',length(chunks),' ...')
    # Output files must exist before path conversion on Windows.
    file.create(tables[i])
    dom <- sub('\\.tbl$','.domtbl',tables[i]); file.create(dom)
    log <- sub('\\.tbl$','.log',tables[i])
    opts <- c('--noali','--cpu','1','-E','0.01','--domE','0.01',
      '--tblout',shQuote(native_path(tables[i])),'--domtblout',shQuote(native_path(dom)),
      shQuote(native_path(chunks[i])),shQuote(native_path(query_file)))
    status <- if(use_wsl) system2('wsl.exe',c('-e','hmmsearch',opts),stdout=log,stderr=log) else system2('hmmsearch',opts,stdout=log,stderr=log)
    if(status!=0 || !any(grepl('^# \\[ok\\]',tail(readLines(tables[i],warn=FALSE),10)))) stop('Incomplete HMMER search: ',i)
    writeLines(c(key,sha(tables[i])),stamp)
  }
}
parse_hits <- function(p) {
  rows <- readLines(p,warn=FALSE); rows <- rows[nzchar(rows)&!startsWith(rows,'#')]
  parts <- strsplit(trimws(rows),'[[:space:]]+')
  if(!length(parts)) return(tibble())
  tibble(accession=vapply(parts,`[`,character(1),1),KO=vapply(parts,`[`,character(1),3),
    full_evalue=as.numeric(vapply(parts,`[`,character(1),5)),full_score=as.numeric(vapply(parts,`[`,character(1),6)),
    domain_evalue=as.numeric(vapply(parts,`[`,character(1),8)),domain_score=as.numeric(vapply(parts,`[`,character(1),9)))
}
hits <- bind_rows(lapply(tables,parse_hits)) %>%
 left_join(as_tibble(ko_ref) %>% select(KO=knum,threshold,score_type,definition),by='KO') %>%
 mutate(score=if_else(score_type=='domain',domain_score,full_score),
        evalue=if_else(score_type=='domain',domain_evalue,full_evalue),
        passes=!is.na(threshold)&score>threshold&evalue<=.01)
write_csv(hits,file.path(ko_dir,'All_reported_KO_hits.csv'))
passed <- hits %>% filter(passes) %>% arrange(accession,KO)
write_csv(passed,file.path(ko_dir,'Threshold_passing_KO_hits.csv'))
# Exact KO membership map, fixed before examining abundance or differential tests.
# This covers the thylakoid functional scope.
ko_groups <- list(
 'Photosystem I'=c(sprintf('K%05d',2689:2702),'K08905','K14332'),
 'Cytochrome b6f'=c('K02634','K02635','K02636','K02637','K02640','K02642','K02643','K03689'),
 'ATP synthase'=c(sprintf('K%05d',2108:2116),'K01549'),
 'NDH complex'=c(sprintf('K%05d',330:343),sprintf('K%05d',5572:5585),'K00345','K23051','K23052'),
 'Terminal oxidase'=c(sprintf('K%05d',2274:2277),sprintf('K%05d',404:407)),
 'Electron carriers'=c('K02638','K02641','K00528','K21567','K08906'),
 'Insertion / proteostasis'=c('K03217','K03798','K03070','K03072','K03073','K03074','K03075','K03076','K03116','K03117','K03118'),
 'Cofactor delivery'=c('K02199','K07399'),
 'Other redox candidates'=c('K02639','K00204','K00205','K17892','K03839'),
 'Pigment / quinone metabolism'=c('K00218','K04034','K04035','K04036','K04037','K04038','K04039','K10960','K03183','K21231','K17872'))
ko_map <- bind_rows(lapply(names(ko_groups),function(g) tibble(KO=ko_groups[[g]],KO_category=g))) %>%
 left_join(as_tibble(ko_ref) %>% select(KO=knum,definition),by='KO')
stopifnot(!anyDuplicated(ko_map$KO),!anyNA(ko_map$definition))
write_csv(ko_map,file.path(ko_dir,'KO_functional_scope.csv'))
# Conservative rule: only a unique threshold-passing KO can add a new candidate.
# Multiple passing KOs remain in the audit.
ko_summary <- passed %>% group_by(accession) %>% summarise(
 KO_count=n_distinct(KO),KO_hits=paste(unique(KO),collapse=';'),
 KO_definition=if(n_distinct(KO)==1) first(definition) else NA_character_,.groups='drop') %>%
 mutate(KO_unique=if_else(KO_count==1,KO_hits,NA_character_)) %>%
 left_join(ko_map %>% select(KO_unique=KO,KO_category),by='KO_unique')
ann <- ann %>% left_join(ko_summary,by='accession') %>% mutate(KO_count=coalesce(KO_count,0L),
 KO_status=case_when(KO_count==0 ~ 'No threshold-passing KO',KO_count==1 ~ 'Unique threshold-passing KO',TRUE ~ 'Multiple threshold-passing KOs; unresolved'))
write_csv(tibble(reference=c('profiles.tar.gz','ko_list.gz','Dryad-derived query',basename(chunks)),
 sha256=c(sha(archive),sha(ko_list_file),query_hash,chunk_hashes)),file.path(ko_dir,'Reference_checksums.csv'))
writeLines(c(hmmer_version,'Reference: https://www.genome.jp/ftp/db/kofam/',
 'Reference downloaded 2026-10-21; full collection of 27930 profiles.',
 'Threshold: score > KO-specific threshold; score type full/domain; corresponding E <= 0.01.',
 'KO prediction is functional evidence, not proof of thylakoid localization.'),file.path(ko_dir,'Annotation_method.txt'))

ann <- ann %>% mutate(category=case_when(
 grepl('photosystem I assembly',description,ignore.case=TRUE) ~ 'PSI assembly',
 grepl('photosystem I ',description,ignore.case=TRUE) ~ 'Photosystem I',
 grepl('cytochrome b6|^cytochrome f$',description,ignore.case=TRUE) ~ 'Cytochrome b6f',
 grepl('ATP synthase',description,ignore.case=TRUE) ~ 'ATP synthase',
 grepl('^NAD.*quinone oxidoreductase',description) ~ 'NDH complex',
 grepl('cytochrome c oxidase',description) ~ 'Terminal oxidase',
 grepl('^plastocyanin$|^ferredoxin--NADP reductase$',description) ~ 'Electron carriers',
 grepl('FtsH|YidC|protein insertion efficiency|translocase|preprotein translocase|SecDF|SecF',description) ~ 'Insertion / proteostasis',
 grepl('cytochrome c biogenesis|PacS',description) ~ 'Cofactor delivery',
 grepl('^ferredoxin($|-)|^NADH dehydrogenase$',description) ~ 'Other redox candidates',
 grepl('protochlorophyllide|Ycf54|demethylphylloquinone|demethylmenaquinone',description) ~ 'Pigment / quinone metabolism',
 TRUE ~ 'Outside defined scope'),
 description_category=category,
 KO_addition=category=='Outside defined scope' & !is.na(KO_category),
 category=if_else(KO_addition,KO_category,category),
 selection_evidence=case_when(KO_addition ~ 'Unique sequence-derived KO',category!='Outside defined scope' & !is.na(KO_category) ~ 'Dryad description and sequence-derived KO',category!='Outside defined scope' ~ 'Dryad description',TRUE ~ 'Outside scope'),
 KO_category_disagreement=!is.na(KO_category) & description_category!='Outside defined scope' & KO_category!=description_category,
 tier=case_when(category %in% c('Photosystem I','Cytochrome b6f','ATP synthase','NDH complex','Terminal oxidase') ~ 'Complex subunit',
 category %in% c('PSI assembly','Electron carriers') ~ 'Specific association',
 category!='Outside defined scope' ~ 'Broader candidate',TRUE ~ 'Not selected'),
 localization_caveat=case_when(tier=='Complex subunit' ~ 'Complex membership inferred from annotation; UCYN-A membrane location not directly measured; not every subunit spans membrane',
 tier=='Specific association' ~ 'Associated function; may be soluble or peripheral, not evidence of integral membrane localization',
 tier=='Broader candidate' ~ 'Possible relationship only; may act in other membranes or soluble pathways',
 TRUE ~ 'No match to the stated functional rules; does not establish absence from thylakoids'))
abundance <- quant_long %>% group_by(PG.ProteinAccessions) %>% summarise(n_observed=n(),
 median_percent=median(riBAQ_percent_UCYNA), median_raw_quantity=median(PG_Quantity),
 observable_tryptic_peptides=first(observable_tryptic_peptides),.groups='drop')
audit <- ann %>% left_join(abundance,by=c('accession'='PG.ProteinAccessions')) %>%
 mutate(n_observed=coalesce(n_observed,0L), main_eligible=tier %in% c('Complex subunit','Specific association') & n_observed>=4)
ranking <- audit %>% filter(main_eligible) %>% arrange(desc(median_percent),accession) %>% mutate(abundance_rank=row_number()) %>%
 mutate(raw_quantity_rank=rank(-median_raw_quantity,ties.method='first'), main_selected=abundance_rank<=15)
audit <- audit %>% left_join(ranking %>% select(accession,abundance_rank,raw_quantity_rank,main_selected),by='accession') %>%
 mutate(main_selected=coalesce(main_selected,FALSE),SI_selected=tier!='Not selected' & n_observed>0,
 selection_reason=case_when(main_selected ~ 'Top 15 eligible pooled peptide-adjusted abundances',
 main_eligible ~ 'Eligible; below main-text abundance cutoff',
 tier=='Not selected' ~ 'No functional-rule match',n_observed==0 ~ 'No usable uniquely mapped public-fraction abundance',
 tier=='Broader candidate' ~ 'Broader functional candidate; SI only',TRUE ~ 'Fewer than four usable public samples; SI only'))
# NifH is an abundant nitrogenase subunit, appropriate for a diazotroph proteome benchmark
audit <- audit %>% mutate(benchmark=accession=='BDA40030.1',
 main_display=main_selected | benchmark, SI_selected=SI_selected | benchmark,
 category=if_else(benchmark,'Nitrogen fixation benchmark',category),
 tier=if_else(benchmark,'Benchmark',tier),
 selection_reason=if_else(benchmark,'User-requested NifH abundance benchmark; outside thylakoid selection',selection_reason),
 localization_caveat=if_else(benchmark,'Nitrogenase benchmark; not a thylakoid-localization claim',localization_caveat))
category_order <- c('Photosystem I','Cytochrome b6f','ATP synthase','Electron carriers',
 'NDH complex','Terminal oxidase','PSI assembly','Nitrogen fixation benchmark',
 'Insertion / proteostasis','Cofactor delivery','Other redox candidates','Pigment / quinone metabolism')
# Statistics are joined only after the selection and ranking are complete.
audit <- audit %>% left_join(limma %>% select(-Origin,-sequence_source),by=c('accession'='PG.ProteinAccessions')) %>%
 mutate(main_display=SI_selected & is.finite(FDR) & is.finite(log2FC_DayMinusNight),
 selection_reason=case_when(benchmark ~ 'NifH benchmark',main_display ~ paste('Testable candidate:',selection_evidence),SI_selected ~ 'Detected candidate; insufficient data for test',tier!='Not selected' ~ 'Functional candidate without usable public abundance',TRUE ~ 'Outside defined functional scope'))
# Display aliases do not affect inclusion. Other labels retain the public annotation.
aliases <- c('BDA40030.1'='NifH','BDA40052.1'='PsaF','BDA39698.1'='PsaC','BDA39703.1'='PsaD','BDA39708.1'='ATP synthase beta',
 'BDA40097.1'='PsaK','BDA39955.1'='FNR (PetH)','BDA39883.1'='PetC','BDA39584.1'='ATP synthase alpha',
 'BDA39736.1'='PsaE','BDA39585.1'='ATP synthase delta','BDA40126.1'='PsaA','BDA40127.1'='PsaB',
 'BDA39882.1'='PetA','BDA39171.1'='PetB','BDA39172.1'='PetD','BDA39678.1'='Plastocyanin (PetE1)',
 'BDA39583.1'='ATP synthase gamma','BDA40160.1'='PsaN','BDA39967.1'='PsaL')
short_label <- function(x) {
 x <- sub('^putative ','',x)
 x <- sub('NAD\\(P\\)H-quinone oxidoreductase','NDH',x)
 x <- sub('NADH-quinone oxidoreductase','NDH',x)
 x <- sub('photosystem I assembly protein','PSI assembly',x)
 x <- sub('ATP-dependent zinc metalloprotease ','',x)
 x <- sub('light-independent protochlorophyllide reductase','Protochlorophyllide reductase',x)
 x <- sub('Sec-independent protein translocase protein ','',x)
 x <- sub('protein translocase subunit ','',x)
 x <- sub('ferredoxin-thioredoxin reductase','Fd-thioredoxin reductase',x)
 x <- sub('iron-sulfur ATP-binding protein','Fe-S ATP-binding',x)
 x
}
audit$label <- short_label(audit$description)
# Preserve public Dryad descriptions; label new KO-only candidates as predictions.
ko_new <- audit$KO_addition
audit$label[ko_new] <- paste0(audit$KO_unique[ko_new], ' (pred.)')
ko_labels <- c(K10960='Geranylgeranyl reductase',K04035='Chlorophyll cyclase',K02639='Fd',K02696='PSI subunit VIII',K00345='NDH-S',K02116='ATP-I',K02642='b6f subunit 6',K03075='SecG',K02275='CoxII')
known <- ko_new & audit$KO_unique %in% names(ko_labels)
audit$label[known] <- paste0(unname(ko_labels[audit$KO_unique[known]]),' (pred.)')
has_alias <- audit$accession %in% names(aliases)
audit$label[has_alias] <- unname(aliases[audit$accession[has_alias]])
atp <- grepl('ATP synthase',audit$description)
audit$label[atp] <- sub('ATP synthase (subunit )?', 'ATP-',audit$description[atp])
audit$label[atp] <- sub(' chain$', '',audit$label[atp])
for(pair in list(c('alpha','\u03b1'),c('beta','\u03b2'),c('gamma','\u03b3'),c('delta','\u03b4'),c('epsilon','\u03b5'))) {
 audit$label[atp] <- sub(pair[1],pair[2],audit$label[atp])
}

write_csv(audit,file.path(outdir,'Genome_wide_selection_audit.csv'))
write_csv(audit %>% filter(main_eligible) %>% arrange(abundance_rank),file.path(outdir,'Abundance_ranking_and_sensitivity.csv'))
write_csv(quant_long,file.path(outdir,'Supplementary_Data_all_quantified_UCYNA_iBAQ_like.csv'))
write_csv(fasta %>% semi_join(ann,by='accession') %>% select(-sequence),file.path(outdir,'Genome_observable_peptides.csv'))
source <- quant_long %>% inner_join(audit,by=c('PG.ProteinAccessions'='accession')) %>%
 left_join(sample_key %>% select(sample,published_column),by='sample')
write_csv(source %>% filter(main_display),file.path(outdir,'Source_Data_Fig4g.csv'))
write_csv(source %>% filter(SI_selected),file.path(outdir,'Source_Data_expanded_SI.csv'))
condition_summary <- source %>% filter(SI_selected) %>% group_by(PG.ProteinAccessions,condition) %>%
 summarise(n=n(),median_percent=median(riBAQ_percent_UCYNA),min_percent=min(riBAQ_percent_UCYNA),max_percent=max(riBAQ_percent_UCYNA),.groups='drop')
write_csv(condition_summary,file.path(outdir,'Condition_abundance_summary.csv'))


# 2. Define curated candidates and their labels; both plots display all testable UCYN-A proteins.
mean_abundance <- quant_long %>% group_by(PG.ProteinAccessions) %>%
 summarise(mean_percent=mean(riBAQ_percent_UCYNA),mean_raw_quantity=mean(PG_Quantity),n_for_mean=n(),.groups='drop')
volcano <- audit %>% filter(SI_selected) %>% left_join(mean_abundance,by=c('accession'='PG.ProteinAccessions')) %>% mutate(
 tested=is.finite(FDR)&is.finite(log2FC_DayMinusNight),
 minus_log10_FDR=if_else(tested,-log10(pmax(FDR,.Machine$double.xmin)),NA_real_),
 significance=case_when(tested & FDR<.05 & log2FC_DayMinusNight<0 ~ 'Night enriched',
 tested & FDR<.05 & log2FC_DayMinusNight>0 ~ 'Day enriched',TRUE ~ 'FDR >= 0.05'),
 point_type=if_else(benchmark,'NifH benchmark','Thylakoid-related candidate'))
# Label significant candidates and NifH in the main-text plot.
# Abbreviations affect display only. Full annotations/categories remain in CSV.
volcano$volcano_label <- volcano$label
volcano$volcano_label <- sub('ferredoxin-thioredoxin reductase','Fd-thioredoxin reductase',volcano$volcano_label)
volcano$volcano_label <- sub('Protochlorophyllide reductase subunit ','Chl reductase ',volcano$volcano_label)
duplicates <- duplicated(volcano$volcano_label)|duplicated(volcano$volcano_label,fromLast=TRUE)
volcano$volcano_label[duplicates] <- paste0(volcano$volcano_label[duplicates],' [',sub('BDA','',volcano$accession[duplicates]),']')
symbol_map <- c('BDA39955.1'='PetH','BDA39678.1'='PetE1',
 'BDA40040.1'='PetF','BDA39441.1'='Fd')
for(a in names(symbol_map)) volcano$volcano_label[volcano$accession==a] <- symbol_map[a]
volcano <- volcano %>% mutate(volcano_main_display=tested)
write_csv(volcano,file.path(outdir,'Source_Data_volcano_all_candidates.csv'))
write_csv(volcano %>% filter(tested),file.path(outdir,'Source_Data_curated_candidates.csv'))
write_csv(volcano %>% filter(!tested),file.path(outdir,'SI_candidates_not_tested.csv'))
plotted <- volcano %>% filter(tested)
xmax <- ceiling(max(abs(plotted$log2FC_DayMinusNight))+.6)
ymax <- max(2,ceiling((max(plotted$minus_log10_FDR)+.45)*2)/2)
palette <- c('Night enriched'='#277DA8','Day enriched'='#C66B0A','FDR >= 0.05'='#A2A2A2')
make_volcano <- function(si=FALSE) {
 d <- plotted
 labels <- d %>% filter((curated_candidate & FDR<.05) | benchmark)
 ggplot(d,aes(log2FC_DayMinusNight,minus_log10_FDR)) +
 geom_hline(yintercept=-log10(.05),linetype='dashed',linewidth=.35,colour='grey60') +
 geom_vline(xintercept=0,linewidth=.35,colour='grey70') +
 geom_point(aes(colour=main_colour,size=mean_percent,alpha=main_emphasis),shape=16) +
 scale_alpha_manual(values=c('Background'=.45,'Highlighted'=.95),guide='none') +
 ggrepel::geom_text_repel(data=labels,aes(label=volcano_label),seed=41,
   size=if(si) 2.8 else 3,colour='#252525',max.overlaps=Inf,
   nudge_y=ifelse(labels$volcano_label %in% c('NDH subunit J','NDH-S (pred.)'),.22,0),
   box.padding=.55,point.padding=.3,min.segment.length=0,
   segment.colour='grey65',segment.size=.3,force=3,
   max.time=5,max.iter=20000,show.legend=FALSE) +
 annotate('text',x=-xmax*.54,y=ymax*.98,label='Higher at night',fontface='bold',size=3.3,hjust=.5) +
 annotate('text',x=xmax*.54,y=ymax*.98,label='Higher during day',fontface='bold',size=3.3,hjust=.5) +
 annotate('text',x=-xmax+.07,y=-log10(.05)+.06,label='FDR = 0.05',hjust=0,size=2.7,colour='grey40') +
 scale_colour_manual(values=main_palette,limits=names(main_palette),name=NULL,drop=FALSE) +
 scale_size_continuous(range=c(1.5,7),limits=c(0,max(plotted$mean_percent)),breaks=c(1,5,20),name='Mean abundance (% UCYN-A)') +
 scale_x_continuous(limits=c(-xmax,xmax),breaks=seq(-xmax,xmax,by=2),expand=expansion(mult=.02)) +
 scale_y_continuous(limits=c(-.15,ymax),breaks=seq(0,ymax,by=1),expand=expansion(mult=c(.02,.05))) +
 labs(x='log2(day / night)',y=expression(-log[10]('FDR')),
 title='UCYN-A proteins: thylakoid-related candidates in context',
 subtitle=paste0(nrow(d),' proteins | colour: significant curated candidates | size: mean abundance')) +
 theme_classic(base_size=10) + theme(axis.text=element_text(colour='black',size=9),
 plot.title=element_text(face='bold',size=12),plot.subtitle=element_text(size=9),
 legend.position='bottom',legend.box='vertical',legend.text=element_text(size=8),
 plot.margin=margin(8,10,6,8)) +
 guides(size=guide_legend(order=2,override.aes=list(colour='#777777')),colour=guide_legend(order=1,override.aes=list(size=3)))
}
stopifnot(sum(volcano$main_display)==sum(volcano$tested),sum(volcano$main_selected)==15,sum(volcano$benchmark)==1)
writeLines(c(paste('Main proteins:',nrow(plotted)),paste('Main candidate/reference proteins plotted:',nrow(plotted)),
 paste('SI candidates without estimable test:',sum(!volcano$tested)),
 'Colour requires BH FDR < 0.05 across the original 925 tests; no fold-change cutoff.',
 'NifH is a separate benchmark. Candidate selection uses public Dryad descriptions plus unique sequence-derived KOs.',
 'All quantities and statistics computed from the six public isolated-UCYN-A samples.'),file.path(outdir,'Volcano_validation.txt'))

# 3. Genome-wide SI: all 495 testable UCYN-A proteins.
means <- quant_long %>% group_by(PG.ProteinAccessions) %>%
 summarise(mean_percent=mean(riBAQ_percent_UCYNA),mean_raw_quantity=mean(PG_Quantity),n_for_mean=n(),.groups='drop')
all_genome <- audit %>% left_join(means,by=c('accession'='PG.ProteinAccessions')) %>% mutate(
 tested=is.finite(FDR)&is.finite(log2FC_DayMinusNight),
 minus_log10_FDR=if_else(tested,-log10(pmax(FDR,.Machine$double.xmin)),NA_real_),
 significance=case_when(tested & FDR<.05 & log2FC_DayMinusNight<0 ~ 'Night enriched',
 tested & FDR<.05 & log2FC_DayMinusNight>0 ~ 'Day enriched',TRUE ~ 'FDR >= 0.05'),
 plot_status=case_when(tested ~ 'Plotted',n_observed>0 ~ 'Not tested: insufficient observations per condition',TRUE ~ 'Not tested: no usable exact-accession abundance'))
d <- all_genome %>% filter(tested)
stopifnot(nrow(d)==sum(limma$Origin=='UCYNA'),all(is.finite(d$mean_percent)))
# Select labels only; every tested genome protein is plotted irrespective of FDR.
label_ids <- d %>% filter(FDR<.05) %>% group_by(significance) %>%
 arrange(FDR,accession,.by_group=TRUE) %>% slice_head(n=8) %>% pull(accession)
label_ids <- union(label_ids,d$accession[(d$SI_selected & d$FDR<.05)|d$benchmark])
d <- d %>% mutate(label_shown=accession %in% label_ids,
 volcano_label=if_else(nchar(label)<=18 & !grepl('hypothetical|putative',label),label,accession))
aliases <- c('BDA39955.1'='PetH','BDA39678.1'='PetE1','BDA40040.1'='PetF','BDA39441.1'='Fd')
for(a in names(aliases)) d$volcano_label[d$accession==a] <- aliases[a]
# Readable labels derived only from the supplied public genome descriptions.
# Do not assign a function to proteins annotated as hypothetical.
annotation_labels <- c('BDA40254.1'='Alanine dehydrogenase',
 'BDA39958.1'='Triosephosphate\nisomerase',
 'BDA39939.1'='Protochlorophyllide\nreductase B',
 'BDA39693.1'='Ribosomal L11',
 'BDA39619.1'='DNA-binding HU',
 'BDA39575.1'='UDP-GlcNAc\nacyltransferase')
for(a in names(annotation_labels)) d$volcano_label[d$accession==a] <- annotation_labels[a]
hypothetical <- grepl('^hypothetical protein',d$description,ignore.case=TRUE)
d$volcano_label[hypothetical] <- paste0('Hypothetical\n',d$accession[hypothetical])
d$volcano_label[d$KO_addition] <- d$label[d$KO_addition]
dup <- duplicated(d$volcano_label)|duplicated(d$volcano_label,fromLast=TRUE)
d$volcano_label[dup] <- d$accession[dup]
write_csv(d,file.path(outdir,'Source_Data_volcano_all_UCYNA.csv'))
write_csv(all_genome,file.path(outdir,'All_genome_testability_audit.csv'))
write_csv(all_genome %>% filter(!tested),file.path(outdir,'Genome_proteins_not_tested.csv'))
write_csv(quant_long,file.path(outdir,'Source_Data_all_UCYNA_by_sample.csv'))
# Record group identifiers that cannot be treated as one exact genome accession.
write_csv(quant %>% filter(Origin=='UCYNA') %>% anti_join(ann,by=c('PG.ProteinAccessions'='accession')),
 file.path(outdir,'Unmapped_or_ambiguous_UCYNA_groups.csv'))
xmax <- ceiling(max(abs(d$log2FC_DayMinusNight)))
ymax <- ceiling((max(d$minus_log10_FDR)+.45)*2)/2
palette <- c('Night enriched'='#277DA8','Day enriched'='#C66B0A','FDR >= 0.05'='#A2A2A2')
# Draw low-significance observations first, so significant circles remain visible.
d <- d %>% arrange(FDR<.05,desc(mean_percent))
fig <- ggplot(d,aes(log2FC_DayMinusNight,minus_log10_FDR)) +
 geom_hline(yintercept=-log10(.05),linetype='dashed',linewidth=.35,colour='grey60') +
 geom_vline(xintercept=0,linewidth=.35,colour='grey70') +
 geom_point(aes(colour=significance,size=mean_percent),alpha=.85) +
 ggrepel::geom_text_repel(aes(label=if_else(label_shown,volcano_label,'')),
  seed=41,size=2.6,colour='#252525',max.overlaps=Inf,box.padding=.45,point.padding=.3,
  min.segment.length=0,segment.colour='grey65',segment.size=.3,force=3,
  max.time=10,max.iter=30000,show.legend=FALSE) +
 annotate('text',x=-xmax*.52,y=ymax*.98,label='Higher at night',fontface='bold',size=3.3) +
 annotate('text',x=xmax*.52,y=ymax*.98,label='Higher during day',fontface='bold',size=3.3) +
 annotate('text',x=-xmax+.08,y=-log10(.05)+.07,label='FDR = 0.05',hjust=0,size=2.7,colour='grey40') +
 scale_colour_manual(values=palette,limits=names(palette),name=NULL,drop=FALSE) +
 scale_size_continuous(range=c(1.5,7),limits=c(0,max(d$mean_percent)),breaks=c(1,5,20),name='Mean abundance (% UCYN-A)') +
 scale_x_continuous(limits=c(-xmax,xmax),breaks=seq(-xmax,xmax,by=2),expand=expansion(mult=.02)) +
 scale_y_continuous(limits=c(-.15,ymax),breaks=seq(0,ymax,by=1),expand=expansion(mult=c(.02,.05))) +
 labs(x='log2(day / night)',y=expression(-log[10]('FDR')),
 title='All testable UCYN-A-encoded proteins',
 subtitle=paste0(nrow(d),' proteins | ',sum(d$FDR<.05),' at FDR < 0.05 | circle size: mean abundance')) +
 theme_classic(base_size=10) + theme(axis.text=element_text(colour='black',size=9),
 plot.title=element_text(face='bold',size=12),plot.subtitle=element_text(size=9),
 legend.position='bottom',legend.box='vertical',legend.text=element_text(size=8),plot.margin=margin(8,10,6,8)) +
 guides(size=guide_legend(order=2,override.aes=list(colour='#777777')),colour=guide_legend(order=1,override.aes=list(size=3)))
ggsave(file.path(output_dir,'SI_all_UCYNA.pdf'),fig,width=210,height=180,units='mm',device=cairo_pdf)
ggsave(file.path(output_dir,'SI_all_UCYNA.png'),fig,width=210,height=180,units='mm',dpi=300,type='cairo')
writeLines(c(paste('Plotted:',nrow(d)),paste('FDR < .05:',sum(d$FDR<.05)),
 paste('Night significant:',sum(d$significance=='Night enriched')),paste('Day significant:',sum(d$significance=='Day enriched')),
 paste('Genome proteins not tested:',sum(!all_genome$tested)),paste('Labels:',sum(d$label_shown)),
 'BH correction remains across 925 tests. No thylakoid, abundance or significance filter applied to plotted membership.'),file.path(outdir,'All_UCYNA_validation.txt'))


# Main text: the same genome-wide point universe, highlighting only significant
# curated functional candidates. Grey includes significant out-of-scope proteins.
main_palette <- c('Night-enriched candidate'='#277DA8',
 'Day-enriched candidate'='#C66B0A','Other / not significant'='#AAAAAA')
plotted <- d %>% mutate(curated_candidate=accession %in% volcano$accession[!volcano$benchmark],
 main_colour=case_when(curated_candidate & FDR<.05 & log2FC_DayMinusNight<0 ~ 'Night-enriched candidate',
 curated_candidate & FDR<.05 & log2FC_DayMinusNight>0 ~ 'Day-enriched candidate',TRUE ~ 'Other / not significant'),
 main_emphasis=if_else(curated_candidate & FDR<.05,'Highlighted','Background'),
 main_label_shown=(curated_candidate & FDR<.05)|benchmark,
 volcano_main_display=TRUE) %>% arrange(main_emphasis=='Highlighted',desc(mean_percent))
label_match <- match(plotted$accession,volcano$accession)
plotted$volcano_label[!is.na(label_match)] <- volcano$volcano_label[label_match[!is.na(label_match)]]
write_csv(plotted,file.path(outdir,'Source_Data_volcano_main.csv'))
write_csv(plotted,file.path(outdir,'Source_Data_Fig4g.csv'))
xmax <- ceiling(max(abs(plotted$log2FC_DayMinusNight)))
ymax <- ceiling((max(plotted$minus_log10_FDR)+.45)*2)/2
for(si in FALSE) {
 figure <- make_volcano(si)
 name <- 'Main_text_Fig4g'
 ggsave(file.path(output_dir,paste0(name,'.pdf')),figure,width=if(si) 190 else 180,height=if(si) 155 else 150,units='mm',device=cairo_pdf)
 ggsave(file.path(output_dir,paste0(name,'.png')),figure,width=if(si) 190 else 180,height=if(si) 155 else 150,units='mm',dpi=300,type='cairo')
}

stopifnot(nrow(plotted)==495,sum(plotted$curated_candidate)==63,
 sum(plotted$main_emphasis=='Highlighted')==12,sum(plotted$main_label_shown)==13)
writeLines(c('Main: 495 testable UCYN-A proteins; curated set: 63 candidates plus NifH.',
 'Main colour: 12 significant curated candidates (5 night, 7 day).',
 'Grey includes nonsignificant proteins AND significant proteins outside the curated set.',
 'SI colour: all 69 significant genome proteins (29 night, 40 day).',
 'Same 925-test BH correction, coordinates and abundance size mapping in both plots.'),
 file.path(outdir,'Volcano_validation.txt'))

# Record package versions.
write_csv(tibble(package=needed,version=vapply(needed,function(x) as.character(packageVersion(x)),character(1))),
 file.path(outdir,'Package_versions.csv'))
writeLines(c(paste('R version:',getRversion()),
 paste0('Main: ',nrow(plotted),' testable UCYN-A proteins with curated highlighting; SI: ',nrow(d),' testable UCYN-A proteins.'),
 'BH correction: 925 tests; SI significant: 69 (40 day, 29 night).',
 'Circle size: mean observed peptide-adjusted relative abundance, with a minimum visible size.'),
 file.path(outdir,'Reproduction_summary.txt'))
# 3. Supplementary workbooks and matching CSVs. 
# of the R analysis. Blank cells mean missing or not applicable.

read_output <- function(name) readr::read_csv(file.path(outdir,name),show_col_types=FALSE)
figure_data <- read_output('Source_Data_volcano_main.csv') %>% transmute(
 `Protein ID`=accession, `Protein description`=description, `Protein function`=category,
 `Day vs night difference (log2)`=log2FC_DayMinusNight,
 `Adjusted P value`=FDR, `Mean abundance (%)`=mean_percent,
 `Samples used for mean`=n_for_mean,
 `Thylakoid-related candidate`=curated_candidate, `NifH reference`=benchmark,
 `Main figure colour group`=main_colour)
tests <- read_output('DayVsNight_limma_results_INTENSITY.csv') %>% transmute(
 `Protein ID`=PG.ProteinAccessions, `Protein origin`=Origin,
 `Day vs night difference (log2)`=log2FC_DayMinusNight,
 `Standard error`=SE, `95% confidence interval lower`=CI_lo,
 `95% confidence interval upper`=CI_hi, `P value`=p_value, `Adjusted P value`=FDR)
quantities <- quant %>% select(`Protein ID`=PG.ProteinAccessions,`Protein origin`=Origin,
 `Day sample 1`=Diazoplast1,`Day sample 2`=Diazoplast2,`Day sample 3`=Diazoplast3,
 `Night sample 1`=Diazoplast4,`Night sample 2`=Diazoplast5,`Night sample 3`=Diazoplast6)
abundance <- read_output('Source_Data_all_UCYNA_by_sample.csv')
abundance_summary <- abundance %>% group_by(PG.ProteinAccessions) %>% summarise(
 `Mean abundance (%)`=mean(riBAQ_percent_UCYNA),`Samples used for mean`=n(),.groups='drop') %>%
 rename(`Protein ID`=PG.ProteinAccessions)
abundance_samples <- abundance %>% transmute(`Protein ID`=PG.ProteinAccessions,
 `Time of day`=condition,`Sample number`=replicate,`Measured quantity`=PG_Quantity,
 `Predicted peptide count`=observable_tryptic_peptides,
 `Quantity divided by peptide count`=iBAQ_like,`Relative abundance (%)`=riBAQ_percent_UCYNA)
s1 <- list(Figure_values=figure_data,Statistical_results=tests,Original_quantities=quantities,
 Abundance_by_sample=abundance_samples,Mean_abundance=abundance_summary)
genome <- read_output('All_genome_testability_audit.csv')
annotations <- genome %>% transmute(`Protein ID`=accession,`Protein description`=description,
 `Protein function`=category,`Accepted KO matches`=KO_hits,`Predicted functions`=KO_definition,
 `KO match status`=KO_status,`Reason for selection`=selection_evidence,
 `Thylakoid-related candidate`=tier!='Not selected' & !benchmark,
 `NifH reference`=benchmark,`Shown in figures`=tested)
matches <- readr::read_csv(file.path(output_dir,'KO_annotation','Threshold_passing_KO_hits.csv'),show_col_types=FALSE) %>%
 transmute(`Protein ID`=accession,`KO ID`=KO,`Predicted function`=definition,
 `Match score`=score,`Required score`=threshold,`Score type`=score_type,`Match E value`=evalue)
category_map <- readr::read_csv(file.path(output_dir,'KO_annotation','KO_functional_scope.csv'),show_col_types=FALSE) %>%
 transmute(`KO ID`=KO,`Protein function`=KO_category,`KO description`=definition)
s2 <- list(Protein_descriptions=annotations,Accepted_KO_matches=matches,KO_categories=category_map)
write_plain_workbook <- function(sheets,path) {
 # Write data and widths directly, without modifying the XLSX archive.
 wb <- openxlsx::createWorkbook()
 for(name in names(sheets)) {
  x <- sheets[[name]]
  openxlsx::addWorksheet(wb,name)
  openxlsx::writeData(wb,name,as.data.frame(x),headerStyle=NULL,keepNA=FALSE)
  widths <- vapply(seq_along(x),function(j) {
   v <- as.character(x[[j]]);v <- v[!is.na(v)]
   min(100,max(18,nchar(names(x)[j])+2,if(length(v)) max(nchar(v))+2 else 0))
  },numeric(1))
  openxlsx::setColWidths(wb,name,cols=seq_along(x),widths=widths)
 }
 openxlsx::saveWorkbook(wb,path,overwrite=TRUE)
}
stopifnot(nrow(figure_data)==495,nrow(tests)==925,nrow(abundance_summary)==605,nrow(annotations)==1186)
used_packages <- sort(loadedNamespaces())
software <- bind_rows(
 tibble(Software='R',Version=as.character(getRversion()),Type='Analysis environment'),
 tibble(Software='HMMER hmmsearch',Version=trimws(sub('^#\\s*','',hmmer_version)),Type='Sequence search'),
 tibble(Software=used_packages,
  Version=vapply(used_packages,function(x) as.character(packageVersion(x)),character(1)),
  Type=ifelse(used_packages %in% needed,'R package','R dependency or base package')))
native_versions <- c(extSoftVersion()[setdiff(names(extSoftVersion()),'BLAS')],grDevices::grSoftVersion())
software <- bind_rows(software,tibble(Software=names(native_versions),Version=unname(native_versions),Type='R supporting library'))
if(use_wsl) {
 kernel_version <- trimws(paste(system2('wsl.exe',c('-e','uname','-r'),stdout=TRUE),collapse=' '))
 software <- bind_rows(software,tibble(Software='WSL Linux kernel',Version=kernel_version,Type='HMMER environment'))
}
software <- software %>% filter(!is.na(Version),nzchar(Version)) %>% distinct(Software,Version,.keep_all=TRUE)
s1$Software <- software
write_csv(software,file.path(outdir,'Software_versions.csv'))
for (id in c('S1','S2')) {
 sheets <- if(id=='S1') s1 else s2
 write_plain_workbook(sheets,file.path(output_dir,paste0('Supplementary_Data_',id,'.xlsx')))
 csv_dir <- file.path(output_dir,'supplementary_data',id)
 dir.create(csv_dir,recursive=TRUE,showWarnings=FALSE)
 for (sheet in names(sheets)) write_csv(sheets[[sheet]],file.path(csv_dir,paste0(sheet,'.csv')),na='')
}
cat('Completed both volcano figures, source tables and plain supplementary datasets S1 and S2.\n')

