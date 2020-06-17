#!/usr/bin/Rscript

library(ggplot2)
library(ape)
library(plyr)
library(reshape2)
library(cluster)
library(RColorBrewer)
library(phyloseq)
library(grid)
library(gridExtra)
library(gplots)
library(vegan)
library(irr)
library(useful)
library(pscl)
library(parallel)
library(poLCA)
library(igraph)
library(randomForest)
library(ROCR)
library(metagenomeSeq)
library(stringi)
library(Rtsne)
library(ggfortify)
library(mixOmics)

source("mcc.R")
source("utils.R")
data_dir <- "data/"
out_dir <- "data/" # change appropriately

distance_metrics <- c("bray", "jaccard", "jsd")
alpha_metrics <- c("Chao1", "Shannon", "Simpson", "Observed")

#########################################################################################################
## load mapping files
mapping <- read.table(sprintf("%s/NAFLDomics_Mapping.combined.txt", data_dir), header=T, as.is=T, sep="\t", comment.char="", row.names=1)
metadata <- read.delim(sprintf("%s/metadata.052218.txt", data_dir), header=T, sep="\t", comment.char=""); colnames(metadata)[1] <- "StudyID"
# impute missing BMI values
metadata$BMI....kg.m.2. <- as.numeric(as.character(metadata$BMI....kg.m.2.))
metadata[which(is.na(metadata$BMI....kg.m.2.)), "BMI....kg.m.2."] <- mean(metadata$BMI....kg.m.2., na.rm=T)
mapping$Cohort <- metadata$Cohort[match(mapping$StudyID, metadata$StudyID)]
mapping$BMI <- metadata[match(mapping$StudyID, metadata$StudyID), "BMI....kg.m.2."]
mapping$Age <- metadata[match(mapping$StudyID, metadata$StudyID), "Age.at.Date.of.Visit"]
mapping$Sex <- metadata[match(mapping$StudyID, metadata$StudyID), "Sex"]

## parse metabolon data and Z-transform
df.metabolon <- list()
metabolon_map <- data.frame()
metabolon_sortorder <- list()
for (st in c("plasma", "fecal")) {
	df.metabolon[[st]] <- list()
	metabolon <- read.table(sprintf("%s/metabolites_%s.txt", data_dir, st), header=T, as.is=T, sep="\t", quote="", comment.char="")
	metabolon_map <- rbind(metabolon_map, metabolon[,c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY")])
	metabolite_levels <- c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY")
	sel <- setdiff(colnames(metabolon), metabolite_levels)
	for (mlevel in metabolite_levels) {
		tmp <- metabolon[,c(mlevel, sel)]
		agg <- aggregate(as.formula(sprintf(". ~ %s", mlevel)), tmp, sum); rownames(agg) <- agg[,mlevel]; agg <- agg[,-1]
		agg <- apply(agg, 1, function(x) (x-mean(x))/sd(x))
		agg <- agg[, setdiff(1:ncol(agg), which(is.na(apply(agg, 2, sd))))] # remove entries with zero variation
		# fix metabolite names as necessary
		colnames(agg) <- gsub("'", "prime", colnames(agg))
		df.metabolon[[st]][[length(df.metabolon[[st]])+1]] <- agg
	}
	names(df.metabolon[[st]]) <- metabolite_levels
	
	tmp <- read.table(sprintf("%s/metabolite_map.%s.txt", data_dir, st), header=T, as.is=T, sep="\t", quote="")
	tmp$BIOCHEMICAL <- gsub("'", "prime", tmp$BIOCHEMICAL)
	metabolon_sortorder[[st]] <- tmp
}
names(df.metabolon) <- c("plasma", "fecal")
metabolon_map <- unique(metabolon_map); rownames(metabolon_map) <- metabolon_map$BIOCHEMICAL
cols.superpathway <- c(brewer.pal(length(unique(metabolon_map$SUPER.PATHWAY)), "Set1"), "#bbbbbb"); names(cols.superpathway) <- c(unique(metabolon_map$SUPER.PATHWAY), "NOT_METABOLITE")

# load from DADA2
seqtab <- readRDS(sprintf("%s/merged_seqtab.120717.rds", data_dir))
rownames(seqtab) <- gsub("_F_filt.fastq.gz", "", rownames(seqtab))
blast_fn <- sprintf("%s/BLAST_results.parsed.txt", data_dir)

# load BLAST taxa
taxa <- read.table(blast_fn, header=F, as.is=T, sep="\t", row.names=1)
inds_with_taxa <- as.numeric(gsub("seq", "", rownames(taxa))) # some sequences had no BLAST hits, so exclude these from the phyloseq object
seqtab <- seqtab[, inds_with_taxa]
blast_taxa <- do.call(rbind, lapply(taxa$V2, function(x) unlist(stri_split_fixed(x, ";"))))
rownames(blast_taxa) <- colnames(seqtab); colnames(blast_taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
ps <- phyloseq(otu_table(t(seqtab), taxa_are_rows=TRUE), sample_data(mapping), tax_table(blast_taxa))
set.seed(prod(dim(seqtab)))

# set up metadata
metadata_variables <- read.table(sprintf("%s/metadata_types.txt", data_dir), header=T, as.is=T, sep="\t", row.names=1)
sel <- intersect(rownames(metadata_variables), colnames(mapping))
metadata_variables <- metadata_variables[sel,, drop=F]
mapping.sel <- mapping[rownames(sample_data(ps)), sel]
## fix column types
for (mvar in rownames(metadata_variables)) {
	if (metadata_variables[mvar, "type"] == "factor") {
		mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
		if (metadata_variables[mvar, "baseline"] != "") {
			mapping.sel[,mvar] <- relevel(mapping.sel[,mvar], metadata_variables[mvar, "baseline"])
		}
	} else if (metadata_variables[mvar, "type"] == "ordered") {
		if (metadata_variables[mvar, "baseline"] == "") {
			lvls <- unique(mapping.sel[,mvar])
		} else {
			lvls <- unlist(strsplit(metadata_variables[mvar, "baseline"], ","))
		}
		mapping.sel[,mvar] <- ordered(mapping.sel[,mvar], levels=lvls)
	} else if (metadata_variables[mvar, "type"] == "numeric") {
		mapping.sel[,mvar] <- as.numeric(as.character(mapping.sel[,mvar]))
	}
}
sample_data(ps) <- mapping.sel

color_table <- read.table(sprintf("%s/taxa_coloring.Genus.050818.txt", data_dir), header=T, as.is=T, sep="\t", comment.char="")
coloring <- color_table$Color
names(coloring) <- color_table$Genus
ordering <- rev(names(coloring))
cols <- colorRampPalette(c("white", "red"), space = "rgb")

ordering.genus <- color_table
ordering.genus$Phylum <- factor(ordering.genus$Phylum)
ordering.genus$Class <- factor(ordering.genus$Class)
inds=order(ordering.genus$Phylum, ordering.genus$Class, ordering.genus$Genus)
ordering.genus <- ordering.genus[inds,]
ordering.genus$Genus <- factor(ordering.genus$Genus, levels=unique(ordering.genus$Genus))
cols <- colorRampPalette(c("white", "red"), space = "rgb")

cols.sig <- c("black", "red", "grey"); names(cols.sig) <- c("ns", "sig", "NA")
cols.cohort <- c("#888888", "orange", "red", "#a400a4"); names(cols.cohort) <- c("NMLBMI", "ObeseNML", "NAFLD", "NASH")

out_pdf <- sprintf("%s/phyloseq_output.%s.%s.pdf", out_dir, "NAFLDomics", format(Sys.Date(), "%m%d%y"))
pdf(out_pdf, width=12)

##########################################################################################################################
## 16S analysis

## PCoA pre-decontam
ps.relative <- transform_sample_counts(ps, function(x) x / sum(x) )
for (distance_metric in distance_metrics) {
	ordi <- ordinate(ps.relative, method = "PCoA", distance = distance_metric)
	p <- plot_ordination(ps.relative, ordi, "samples", color = "SampleType2") + geom_point(size = 2) + scale_fill_brewer(type = "qual", palette = "Set1") + scale_colour_brewer(type = "qual", palette = "Set1") + theme_classic() + ggtitle(sprintf("Overall PCoA - %s", distance_metric)) + geom_text(aes(label=StudyID), hjust=-0.4, vjust=0.1, size=3)
	print(p)
}

## SV filtering by %blank
ps.sel <- ps
otu_table <- as.matrix(as.data.frame(otu_table(ps.sel)))
negids <- rownames(subset(mapping, SampleType2 %in% c("BufferControl", "PCRBlank", "PCRWater")))
rs <- rowSums(otu_table[,negids])
rs.true <- rowSums(otu_table(ps)[, setdiff(sample_names(ps), negids)])
pct_blank <- 100* (rs / (rs + rs.true))
hist(pct_blank, breaks=100)
otus_to_remove <- names(which(pct_blank > 10))
# read counts after removing contaminant SVs
ps.pctblank <- prune_taxa(setdiff(taxa_names(ps), otus_to_remove), ps)
df <- melt(sample_sums(ps.pctblank)); df$SampleID <- rownames(df); df <- df[order(df$value),]; df$SampleID <- factor(df$SampleID, levels=df$SampleID)
df$SampleType2 <- ifelse(mapping.sel[as.character(df$SampleID), "SampleType2"] %in% c("BufferControl", "PCRBlank", "PCRWater"), "Blank", "Sample")
p <- ggplot(df, aes(x=SampleID, y=value, fill=SampleType2)) + geom_bar(stat="identity") + theme_classic() + ggtitle("Read counts after contaminant SV removal (%% blank)") + scale_fill_manual(values=c("red", "black"))
print(p)
# PCoA after %blank removal
ps.pctblank <- subset_samples(ps.pctblank, sample_sums(ps.pctblank) > 0)
ps.relative <- transform_sample_counts(ps.pctblank, function(x) x / sum(x) )
for (distance_metric in distance_metrics) {
	ordi <- ordinate(ps.relative, method = "PCoA", distance = distance_metric)
	p <- plot_ordination(ps.relative, ordi, "samples", color = "SampleType2") + geom_point(size = 2) + scale_fill_brewer(type = "qual", palette = "Set1") + scale_colour_brewer(type = "qual", palette = "Set1") + theme_classic() + ggtitle(sprintf("Overall PCoA - %s", distance_metric)) + geom_text(aes(label=StudyID), hjust=-0.4, vjust=0.1, size=3)
	print(p)
}

## final setup of phyloseq objects
ps <- subset_samples(ps, SampleType2=="OralSwab" | SampleType2=="RectalSwab")
ps.relative <- transform_sample_counts(ps, function(x) x / sum(x) )
ps.rarefied <- rarefy_even_depth(ps, sample.size=11712, rngseed=sum(dim(mapping)))
ps <- prune_samples(sample_names(ps.rarefied), ps)
ps.relative <- prune_samples(sample_names(ps.rarefied), ps.relative)

## PCoA (overall)
for (distance_metric in distance_metrics) {
    ordi <- ordinate(ps.relative, method = "PCoA", distance = distance_metric)
    p <- plot_ordination(ps.relative, ordi, "samples", color = "SampleType") + geom_point(size = 2) + scale_fill_brewer(type = "qual", palette = "Set1") + scale_colour_brewer(type = "qual", palette = "Set1") + theme_classic() + ggtitle(sprintf("Overall PCoA - %s", distance_metric)) + geom_text(aes(label=StudyID), hjust=-0.4, vjust=0.1, size=3)
    print(p)
  }

## taxa barplots
for (psvar in c("OralSwab", "RectalSwab")){
  ps.sel <- subset_samples(ps.relative, SampleType==psvar)
  otu.filt <- normalizeByCols(as.data.frame(otu_table(ps.sel)))
  otu.filt$Family <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level="Family")
  agg <- aggregate(. ~ Family, otu.filt, sum)
  families <- agg$Family
  agg <- agg[,-1]
  agg <- normalizeByCols(agg)
  inds_to_grey <- which(rowMeans(agg)<0.005 & !(apply(agg, 1, function(x) any(x>0.05))))
  families[inds_to_grey] <- "Other"
  agg$Family <- families
  df <- melt(agg, variable.name="SampleID")
  df2 <- aggregate(as.formula("value ~ Family + SampleID"), df, sum)
  p <- ggplot(df2, aes_string(x="SampleID", y="value", fill="Family", order="Family")) + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.0)) + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring) + ggtitle(sprintf("Taxa barplots overall - %s", psvar)) + guides(col = guide_legend(ncol = 3))
  print(p)
  
  for (cohort in levels(as(sample_data(ps.sel), "data.frame")$Cohort)){
    ps.sel2 <- subset_samples(ps.sel, Cohort==cohort)
    otu.filt <- normalizeByCols(as.data.frame(otu_table(ps.sel2)))
    otu.filt$Family <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel2), level="Family")
    agg <- aggregate(. ~ Family, otu.filt, sum)
    families <- agg$Family
    agg <- agg[,-1]
    agg <- normalizeByCols(agg)
    inds_to_grey <- which(rowMeans(agg)<0.005 & !(apply(agg, 1, function(x) any(x>0.05))))
    families[inds_to_grey] <- "Other"
    agg$Family <- families
    df <- melt(agg, variable.name="SampleID")
    df2 <- aggregate(as.formula("value ~ Family + SampleID"), df, sum)
    p <- ggplot(df2, aes_string(x="SampleID", y="value", fill="Family", order="Family")) + geom_bar(stat="identity", position="stack") + ylim(c(-0.1, 1.0)) + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring) + ggtitle(sprintf("Taxa barplots %s - %s", cohort, psvar)) + guides(col = guide_legend(ncol = 3))
    print(p)
    
    ps.selr <- subset_samples(ps.rarefied, SampleType==psvar & Cohort==cohort)
    adiv <- estimate_richness(ps.selr)
    rownames(adiv) <- gsub("^X", "", rownames(adiv))
    adiv <- adiv[levels(df2$SampleID),]
    adiv$SampleID <- rownames(adiv)
    adiv$SampleID <- factor(adiv$SampleID, levels=levels(df2$SampleID))
    
    for (am in alpha_metrics){
      p <- ggplot(adiv, aes_string(x="SampleID", y=am)) + geom_bar(stat="identity", position="identity") + theme_classic() + ggtitle(sprintf("Alpha Diversity %s - %s (%s)", am, cohort, psvar)) + theme( axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8))
      print(p)
    }
  }
}

## PCoA and PERMANOVA (oral/rectal separately)
for (psvar in c("OralSwab", "RectalSwab")){
  ps.sel <- subset_samples(ps.relative, SampleType==psvar)
  for (var in c("Cohort")){
    for (distance_metric in distance_metrics) {
      ordi <- ordinate(ps.sel, method = "PCoA", distance = distance_metric)
        p <- plot_ordination(ps.sel, ordi, "samples", color = var) + geom_point(size = 2) + scale_fill_brewer(type = "qual", palette = "Set1") + scale_colour_brewer(type = "qual", palette = "Set1") + theme_classic() + ggtitle(sprintf("By %s - %s (%s)", var, distance_metric, psvar))
        print(p)
      }
    }
}
sink(sprintf("%s/PERMANOVA.16S.txt", out_dir))
for (psvar in c("OralSwab", "RectalSwab")){
  ps.sel <- subset_samples(ps.relative, SampleType==psvar)
  df <- as(sample_data(ps.sel), "data.frame")
  dm <- list()
  for (distance_metric in distance_metrics) {
    dm[[length(dm)+1]] <- as.matrix(distance(ps.sel, method=distance_metric))
  }
  names(dm) <- distance_metrics
  print(sprintf("PERMANOVA manually selected variables (%s)", psvar))
  for (distance_metric in distance_metrics) {
    print(distance_metric)
    form <- as.formula("as.dist(dm[[distance_metric]]) ~  Cohort")
    res <- adonis(form , data=df, permutations=999)
    print(res)
  }
}
sink()
  
## ZINB
for (psvar in c("OralSwab", "RectalSwab")){
  ps.sel <- subset_samples(ps.rarefied, SampleType==psvar)
  for (level in c("Species", "Genus")){
    otu.filt <- as.data.frame(otu_table(ps.sel))
    otu.filt.abundance <- normalizeByCols(otu.filt)
    otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
    otu.filt.abundance[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
    agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
    agg.abundance <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt.abundance, sum)
    lvl <- agg[[level]]
    lvl.abundance <- agg.abundance[[level]]
    agg <- agg[,-1]
    agg.abundance <- agg.abundance[,-1]
    rownames(agg) <- lvl
    rownames(agg.abundance) <- lvl.abundance
    agg.abundance <- agg.abundance[which(rowSums(agg.abundance > 0.01) >= 3),]
    agg <- agg[rownames(agg.abundance),]
    agg[[level]] <- rownames(agg)
    out <- mclapply(agg[[level]], function(f) {
      df <- melt(agg[f,]); colnames(df) <- sprintf(c("%s", "SampleID", "value"), level); df$SampleID <- as.character(df$SampleID)
      tmp <- {}
      print(f)
      for (mvar in c("Cohort")) {
        #print(sprintf("%s %s", f, mvar))
        df2 <- df
        df2[, mvar] <- mapping.sel[df2$SampleID, mvar]
        tt <- try(m <- zeroinfl(as.formula(sprintf("value ~ %s | 1", mvar)), data = df2, dist = "negbin", EM = FALSE, maxit=100), silent=T) # using EM=TRUE causes certain models to hang...
        if (class(tt) == "zeroinfl" && tt$converged) {
          coef <- summary(m)$coefficients$count # rows are [(Intercept), comparisons, Log(theta)], columns are [Estimate, SE, Z, pval]
          tmp <- rbind(tmp, cbind(f, mvar, class(mapping.sel[,mvar]), "ZINB", rownames(coef), coef))
        } else {
          tt <- try(m <- glm.nb(as.formula(sprintf("value ~ %s", mvar)), data = df2), silent=T)
          if (class(tt)[1] == "negbin") { 
            coef <- summary(m)$coefficients # rows are [(Intercept), comparisons], columns are [Estimate, SE, Z, pval]
            tmp <- rbind(tmp, cbind(f, mvar, class(mapping.sel[,mvar]), "NB", rownames(coef), coef))
          }
        }
        #print(sprintf("finished %s %s", f, mvar))
      }
      print(sprintf("finished %s", f))
      tmp
    }, mc.cores=16)
    res <- as.data.frame(do.call(rbind, out))
    colnames(res) <- sprintf(c("%s", "metadata_variable", "class", "method", "coefficient", "Estimate", "SE", "Z", "pvalue"), level)
    res <- subset(res, !coefficient %in% c("(Intercept)", "Log(theta)")) # remove intercept and log-theta terms before FDR correction
    res$pvalue <- as.numeric(as.character(res$pvalue))
    ## p-value adjust
    res$padj <- p.adjust(res$pvalue, method="fdr")
    res <- res[order(res$padj, decreasing=F),]
    resSig <- subset(res, padj<0.05)
    write.table(res, file=sprintf("%s/ZINB_Cohort.%s.%s.txt", out_dir, psvar, level), quote=F, sep="\t", row.names=F, col.names=T)
  }
}

# randomForest classification of Cohort (one-vs-one setup, X vs NMLBMI baseline)
for (psvar in c("OralSwab", "RectalSwab")){
	ps.sel <- subset_samples(ps.relative, SampleType2==psvar)
	for (level in c("Species", "Genus")){
	  otu.filt <- as.data.frame(otu_table(ps.sel))
	  otu.filt.abundance <- normalizeByCols(otu.filt)
	  otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
	  otu.filt.abundance[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
	  # rename Prevotella_6, etc -> Prevotella
	  otu.filt[[level]] <- gsub("_\\d$", "", otu.filt[[level]])
	  otu.filt.abundance[[level]] <- gsub("_\\d$", "", otu.filt.abundance[[level]])
	  agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
	  agg.abundance <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt.abundance, sum)
	  lvl <- agg[[level]]
	  lvl.abundance <- agg.abundance[[level]]
	  agg <- agg[,-1]
	  agg.abundance <- agg.abundance[,-1]
	  rownames(agg) <- lvl
	  rownames(agg.abundance) <- lvl.abundance
	  agg.abundance <- agg.abundance[which(rowSums(agg.abundance > 0.01) >= ceiling(ncol(agg.abundance)/10)),]
	  agg <- agg[rownames(agg.abundance),]
	  res.mean <- {}; res.sd <- {}
  	for (mvar_level in setdiff(levels(mapping.sel[, "Cohort"]), "NMLBMI")) {
			sel <- sample_names(subset_samples(ps.relative, SampleType2==psvar & Cohort %in% c(mvar_level, "NMLBMI")))
			data.sel <- as.data.frame(t(agg[,sel]))
			data.sel$BMI <- as(sample_data(ps.sel), "data.frame")[rownames(data.sel), "BMI"]
			data.sel <- as.matrix(data.sel)
			response <- droplevels(mapping.sel[sel, "Cohort"]); names(response) <- sel

			## after running for the first time, COMMENT OUT THIS BLOCK ##
			num_iter <- 100
			ncores <- 20
			out <- mclapply(1:num_iter, function (dummy) {
					importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
				}, mc.cores=ncores )
			collated.importance <- do.call(cbind, out)
			out <- mclapply(1:num_iter, function (dummy) {
					rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
				}, mc.cores=ncores )
			collated.cv <- do.call(cbind, out)

			write.table(collated.importance, file=sprintf("%s/randomForest.%s_%s.%s.importance.txt", out_dir, psvar, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=F)
			write.table(collated.cv, file=sprintf("%s/randomForest.%s_%s.%s.cv.txt", out_dir, psvar, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=F)
			## END BLOCK TO COMMENT ##
			
			collated.importance <- read.table(sprintf("%s/randomForest.%s_%s.%s.importance.txt", out_dir, psvar, mvar_level, level), header=F, as.is=T, sep="\t", row.names=1)
			collated.cv <- read.table(sprintf("%s/randomForest.%s_%s.%s.cv.txt", out_dir, psvar, mvar_level, level), header=F, as.is=T, sep="\t", row.names=1)
			importance.mean <- rowMeans(collated.importance)
			importance.sd <- unlist(apply(collated.importance, 1, sd))
			cv.mean <- rowMeans(collated.cv)
			cv.sd <- unlist(apply(collated.cv, 1, sd))
			res.mean <- rbind(res.mean, importance.mean)
			res.sd <- rbind(res.sd, importance.sd)
			inds <- order(importance.mean, decreasing=T)
			inds <- inds[1:as.numeric(names(which.min(cv.mean)))] # edit as appropriate
			
			## after running for the first time, COMMENT OUT THIS BLOCK ##
			# using a sparse model with N predictors
			sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
			save(sparseRF, file=sprintf("%s/randomForest.%s_%s.%s.model", out_dir, psvar, mvar_level, level))
			# ROC and AUC of final sparseRF model
			pred <- predict(sparseRF, type="prob")
			pred2 <- prediction(pred[,2], ordered(response))
			perf <- performance(pred2, "tpr", "fpr")
			perf.auc <- performance(pred2, "auc")
			plot(perf, main=sprintf("ROC %s %s %s (sparseRF final model)", psvar, mvar_level, level)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", unlist(perf.auc@y.values)))		
			# 5x cross-validation ROC and AUC
			nperm <- 100
			out <- mclapply(1:nperm, function (dummy) {
				tmp <- data.frame(SampleID=rownames(data.sel), Cohort=droplevels(mapping.sel[rownames(data.sel), "Cohort"]))
				tmp$SampleID <- as.character(tmp$SampleID)
				train.sel <- unlist(lapply(levels(tmp$Cohort), function(lvl) {
					tmp2 <- subset(tmp, Cohort==lvl)
					sample(tmp2$SampleID, size=round(nrow(tmp2)*0.8))
				}))
#				train.sel <- sample(rownames(data.sel), size=round(nrow(data.sel)*0.8))
				test.sel <- setdiff(rownames(data.sel), train.sel)
				rf <- randomForest(x=data.sel[train.sel, names(importance.mean[inds])], y=response[train.sel], ntree=10000, importance=T)
				pred <- predict(rf, data.sel[test.sel, names(importance.mean[inds])], type="prob")
				pred[,2]
			}, mc.cores=ncores )
			pred <- prediction(out, lapply(out, function(x) response[names(x)]) )
			perf <- performance(pred, "tpr", "fpr")
			perf.auc <- performance(pred, "auc"); auc.mean <- mean(unlist(perf.auc@y.values))
#			plot(perf, main="ROC plot")
			plot(perf,col="grey82",lty=3, main=sprintf("ROC %s %s %s (5x cross-validation)", psvar, mvar_level, level))
			plot(perf, lwd=3, avg="vertical", spread.estimate="boxplot", add=TRUE) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", auc.mean))		
			## END BLOCK TO COMMENT ##
			
			load(sprintf("%s/randomForest.%s_%s.%s.model", out_dir, psvar, mvar_level, level))
			# plotting - per-group sparse model
			df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
			colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
			print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s", psvar, mvar_level, level)))
			df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
			print(ggplot(df, aes(x=OTU, y=importance, label=OTU)) + geom_bar(position=position_dodge(), stat="identity", fill="#999999") + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory %s", psvar, mvar_level, level)))
		}
		# plotting - all groups importance values
		rownames(res.mean) <- setdiff(levels(mapping.sel[, "Cohort"]), "NMLBMI"); rownames(res.sd) <- rownames(res.mean)
		df <- res.mean; df[which(df<0.001, arr.ind=T)] <- 0
		heatmap.2(df, Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(10,10), cexCol=1, cexRow=1, main=sprintf("RF importance values (%s, %s)", psvar, level))
		
	}
}

## Bacteroidetes/Firmicutes	ratio
level <- "Phylum"
ps.sel <- subset_samples(ps.relative, SampleType=="RectalSwab")
otu.filt <- as.data.frame(otu_table(ps.sel))
otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
lvl <- agg[[level]]
agg <- agg[,-1]
rownames(agg) <- lvl
agg <- agg[which(rowSums(agg > 0.01) >= ceiling(ncol(agg)/10)),]
agg <- normalizeByCols(agg)
df <- melt(agg["Bacteroidetes",] / agg["Firmicutes",]); colnames(df) <- c("SampleID", "BFRatio"); df$SampleID <- as.character(df$SampleID); df$BFLogRatio <- log10(df$BFRatio)
df$Cohort <- mapping[df$SampleID, "Cohort"]; df$Cohort <- factor(as.character(df$Cohort), levels=c("NMLBMI", "ObeseNML", "NAFLD", "NASH"))
df <- subset(df, !is.na(df$Cohort))
test <- kruskal.test(BFLogRatio ~ Cohort, df)
p <- ggplot(df, aes(x=Cohort, y=BFLogRatio)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("BFLogRatio by Cohort (KW p=%.4g)", test$p.value))
print(p)





###############################################################################################################
### METABOLITES

## linear regression of METABOLITE data ~ Cohort + BMI + Age + Sex
for (st in c("fecal", "plasma")) {
	for (mlevel in metabolite_levels) {
		data <- df.metabolon[[st]][[mlevel]]
		mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% rownames(data))
		out <- mclapply(colnames(data), function(metabolite) {
			df <- melt(data[,metabolite])	
			df[, "Cohort"] <- unlist(sample_data(ps)[match(rownames(df), sample_data(ps)$StudyID), "Cohort"])
			df[, "BMI"] <- unlist(sample_data(ps)[match(rownames(df), sample_data(ps)$StudyID), "BMI"])
			df[, "Age"] <- unlist(sample_data(ps)[match(rownames(df), sample_data(ps)$StudyID), "Age"])
			df[, "Sex"] <- unlist(sample_data(ps)[match(rownames(df), sample_data(ps)$StudyID), "Sex"])
			tt <- try(m <- lm(as.formula("value ~ Cohort + BMI + Age + Sex"), data = df), silent=T)
			if (class(tt) == "lm") {
		    coef <- summary(m)$coefficients # rows are [(Intercept), comparisons], columns are [Estimate, SE, Z, pval]
		    tmp <- cbind(metabolite, "Linear regression", rownames(coef), coef)
			}
			print(sprintf("finished %s", metabolite))
			tmp
		}, mc.cores=16)
		res <- as.data.frame(do.call(rbind, out))
		colnames(res) <- c("metabolite", "method", "coefficient", "Estimate", "SE", "Z", "pvalue")
		res <- subset(res, !coefficient %in% c("(Intercept)")) # remove intercept terms before FDR correction
		res$pvalue <- as.numeric(as.character(res$pvalue))
		## p-value adjust
		res$padj <- p.adjust(res$pvalue, method="fdr")
		res <- res[order(res$padj, decreasing=F),]
		resSig <- subset(res, padj<0.05)
		write.table(res, file=sprintf("%s/METABOLITE_linear_regression.%s.%s.txt", out_dir, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
	}
}

## ordination (t-SNE, PCA) and PERMANOVA
set.seed(112917)
for (st in c("fecal", "plasma")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% rownames(data))
	agg <- aggregate(Cohort ~ StudyID, mapping.sel, unique); agg2 <- aggregate(cbind(BMI, Age) ~ StudyID, mapping.sel, unique); agg3 <- aggregate(Sex ~ StudyID, mapping.sel, unique)
	mapping.sel <- merge(merge(agg, agg2, by="StudyID"), agg3, by="StudyID"); rownames(mapping.sel) <- mapping.sel$StudyID
	# PCA
	pca <- prcomp(data, center=F, scale=F)
	eigs <- pca$sdev^2
	pvar <- 100*(eigs / sum(eigs))
	df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], SampleID=rownames(pca$x))
	for (mvar in rownames(subset(metadata_variables, useForPERMANOVA=="yes"))) {
		df[, mvar] <- mapping.sel[rownames(df), mvar]
		if (mvar == "Cohort") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s PCoA (Euclidean distance)", st, mlevel)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_manual(values=cols.cohort) + stat_ellipse(type="t")
		} else if (metadata_variables[mvar, "type"] == "numeric") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s PCoA (Euclidean distance)", st, mlevel)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_gradient(low="grey", high="red")
		}
		else {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s PCoA (Euclidean distance)", st, mlevel)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + stat_ellipse(type="t")
		}
		print(p)
	}
	# PERMANOVA
	res <- adonis(data ~ BMI + Age + Sex + Cohort, data=mapping.sel, permutations=999, method="eu")
	sink(sprintf("%s/METABOLITE_PERMANOVA.%s.%s.txt", out_dir, st, mlevel))
	print(res)
	sink()
	# t-SNE
	tsne.out <- Rtsne(data, perplexity=20)
	df <- data.frame(PC1 = tsne.out$Y[,1], PC2 = tsne.out$Y[,2], SampleID=rownames(data)); rownames(df) <- df$SampleID
	for (mvar in rownames(subset(metadata_variables, useForPERMANOVA=="yes"))) {
		df[, mvar] <- mapping.sel[rownames(df), mvar]
		if (mvar == "Cohort") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s tSNE", st, mlevel)) + scale_color_manual(values=cols.cohort) + stat_ellipse(type="t")
		} else if (metadata_variables[mvar, "type"] == "numeric") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s tSNE", st, mlevel)) + scale_color_gradient(low="grey", high="red")
		}
		else {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s tSNE", st, mlevel)) + stat_ellipse(type="t")
		}
		print(p)
	}
}


## randomForest classification of Cohort (one-vs-one setup, X vs NMLBMI baseline); using METABOLITE data [plasma, fecal]
set.seed(112717)	
# MAIN LOOP for random forest (through metabolite levels)
for (st in c("fecal", "plasma")) {
	for (mlevel in metabolite_levels) {
		data <- df.metabolon[[st]][[mlevel]]
		mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% rownames(data))
		# random forest classification
		res.mean <- {}; res.sd <- {}; res.sig <- list()
		for (mvar_level in setdiff(levels(mapping.sel[, "Cohort"]), "NMLBMI")) {
			sel <- intersect(rownames(data), unique(subset(sample_data(ps), Cohort %in% c(mvar_level, "NMLBMI"))$StudyID))
			data.sel <- as.data.frame(data[sel,])
			response <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Cohort"]); names(response) <- sel
			# add BMI, Age, Sex as covariates
			data.sel$BMI <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "BMI"])
			data.sel$Age <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Age"])
			data.sel$Sex <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Sex"])
			agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Sex")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
			
			## after running for the first time, COMMENT OUT THIS BLOCK ##
			num_iter <- 100
			ncores <- 20
			out <- mclapply(1:num_iter, function (dummy) {
					importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
				}, mc.cores=ncores )
			collated.importance <- do.call(cbind, out)
			out <- mclapply(1:num_iter, function (dummy) {
					rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
				}, mc.cores=ncores )
			collated.cv <- do.call(cbind, out)

			write.table(collated.importance, file=sprintf("%s/randomForest_METABOLITE.%s_%s.%s.importance.txt", out_dir, mvar_level, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
			write.table(collated.cv, file=sprintf("%s/randomForest_METABOLITE.%s_%s.%s.cv.txt", out_dir, mvar_level, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
			## END BLOCK TO COMMENT ##

			collated.importance <- read.table(sprintf("%s/randomForest_METABOLITE.%s_%s.%s.importance.txt", out_dir, mvar_level, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
			collated.cv <- read.table(sprintf("%s/randomForest_METABOLITE.%s_%s.%s.cv.txt", out_dir, mvar_level, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
			importance.mean <- rowMeans(collated.importance)
			importance.sd <- unlist(apply(collated.importance, 1, sd))
			cv.mean <- rowMeans(collated.cv)
			cv.sd <- unlist(apply(collated.cv, 1, sd))
			res.mean <- rbind(res.mean, importance.mean)
			res.sd <- rbind(res.sd, importance.sd)
			inds <- order(importance.mean, decreasing=T)
			inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.001)))] # edit as appropriate
			res.sig[[length(res.sig)+1]] <- names(importance.mean[inds])
			write.table(melt(importance.mean[inds]), file=sprintf("%s/randomForest_METABOLITE.%s_%s.%s.features.txt", out_dir, mvar_level, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

			## after running for the first time, COMMENT OUT THIS BLOCK ##
			# using a sparse model with N predictors
			sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
			save(sparseRF, file=sprintf("%s/randomForest_METABOLITE.%s_%s.%s.model", out_dir, mvar_level, mlevel, st))
			load(sprintf("%s/randomForest_METABOLITE.%s_%s.%s.model", out_dir, mvar_level, mlevel, st))
			# ROC and AUC of final sparseRF model
			pred <- predict(sparseRF, type="prob")
			pred2 <- prediction(pred[,2], ordered(response))
			perf <- performance(pred2, "tpr", "fpr")
			perf.auc <- performance(pred2, "auc")
			pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
			confusion_matrix <- table(pred_df[, c("true", "predicted")])
			accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
			vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
			mccvalue <- mcc(vec.pred, vec.true)
			plot(perf, main=sprintf("ROC %s %s %s (sparseRF final model)", mvar_level, mlevel, st)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g, accuracy=%.2f%%, MCC=%.4g", unlist(perf.auc@y.values), accuracy, mccvalue))			
			## END BLOCK TO COMMENT ##
			
			# plotting - per-group sparse model
			df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
			colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
			print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s", mvar_level, mlevel, st)))
			# plotting - per-group variables
			df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
			df$metabolite_name <- as.character(df$OTU)
			res <- read.table(sprintf("%s/METABOLITE_linear_regression.%s.%s.txt", out_dir, st, mlevel), header=T, sep="\t", as.is=T, quote="")
			res <- subset(res, coefficient == paste("Cohort", mvar_level, sep="")); rownames(res) <- res$metabolite
			df$Z <- res[as.character(df$OTU), "Z"]; df$padj <- res[as.character(df$OTU), "padj"]
			if (mlevel == "BIOCHEMICAL") {
				df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
				df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
				df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g; %s; %s)", df$OTU, df$Z, df$padj, df$subpathway, df$superpathway)
			} else if (mlevel == "SUB.PATHWAY") {
				df$subpathway <- df$metabolite_name
				df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
				df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g)", df$OTU, df$Z, df$padj)
			} else {
				df$superpathway <- df$metabolite_name
				df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g)", df$OTU, df$Z, df$padj)
			}
			df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
			df$sig <- ifelse(df$padj < 0.1, "sig", "ns"); df$sig[which(is.na(df$sig))] <- "ns"
			p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string, color=sig), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory %s", mvar_level, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
			print(p)
			# shading rectangles of importance
			df.rect <- df
			df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
			df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
			p <- ggplot(df.rect, aes(x=x, y=OTU, color=sig, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s explanatory %s", mvar_level, mlevel, st)) + scale_fill_gradient(low="white", high="black") + scale_color_manual(values=cols.sig)
			print(p)
			# violin plots of relative abundance
			agg.melt <- agg.melt.stored
			agg.melt$Cohort <- mapping.sel[match(agg.melt$SampleID, mapping.sel$StudyID), "Cohort"]
			agg.melt <- subset(agg.melt, metabolite %in% rownames(df) & Cohort %in% c(mvar_level, "NMLBMI"))
			agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
			p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mvar_level, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)		
			
			## custom plot to show value for a single misclassified subject (L099)
			agg.melt$L099 <- ifelse(agg.melt$SampleID=="L099", "L099", "not L099")
			p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=L099)) + geom_violin() + geom_point(size=1) + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mlevel, st)) + coord_flip() + scale_color_manual(values=c("red", "grey"))
			print(p)
			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Cohort)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
			agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
			p <- ggplot(agg.melt2, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Cohort)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
			print(p)
			## END ADDED	
		}
		# plotting - all groups importance values
		rownames(res.mean) <- setdiff(levels(mapping.sel[, "Cohort"]), "NMLBMI"); rownames(res.sd) <- rownames(res.mean)
		df <- res.mean[, unique(unlist(res.sig))]; df[which(df<0.001, arr.ind=T)] <- 0
		df.sig <- matrix("", nrow=nrow(df), ncol=ncol(df)); rownames(df.sig) <- rownames(df); colnames(df.sig) <- colnames(df)
		for (i in 1:length(res.sig)) {
			df.sig[i, res.sig[[i]]] <- "*"
		}
		heatmap.2(df, Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(10,10), cexCol=0.6, cexRow=0.6, cellnote=df.sig, notecex=0.8, notecol="red", main=sprintf("RF importance values (%s, %s)", st, mlevel))
	}
}

## randomForest classification of Cohort (one-vs-one setup, X vs ObeseNML baseline); using METABOLITE data [plasma, fecal]
set.seed(041918)	
# MAIN LOOP for random forest (through metabolite levels)
for (st in c("fecal", "plasma")) {
	for (mlevel in metabolite_levels) {
		data <- df.metabolon[[st]][[mlevel]]
		mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% rownames(data))
		# random forest classification
		res.mean <- {}; res.sd <- {}; res.sig <- list()
		for (mvar_level in setdiff(levels(mapping.sel[, "Cohort"]), c("NMLBMI", "ObeseNML"))) {
			sel <- intersect(rownames(data), unique(subset(sample_data(ps), Cohort %in% c(mvar_level, "ObeseNML"))$StudyID))
			data.sel <- as.data.frame(data[sel,])
			response <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Cohort"]); names(response) <- sel
			# add BMI, Age, Sex as covariates
			data.sel$BMI <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "BMI"])
			data.sel$Age <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Age"])
			data.sel$Sex <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Sex"])
			agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Sex")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
			
			## after running for the first time, COMMENT OUT THIS BLOCK ##
			num_iter <- 100
			ncores <- 20
			out <- mclapply(1:num_iter, function (dummy) {
					importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
				}, mc.cores=ncores )
			collated.importance <- do.call(cbind, out)
			out <- mclapply(1:num_iter, function (dummy) {
					rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
				}, mc.cores=ncores )
			collated.cv <- do.call(cbind, out)

			write.table(collated.importance, file=sprintf("%s/randomForest_METABOLITE_ObeseNML.%s_%s.%s.importance.txt", out_dir, mvar_level, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
			write.table(collated.cv, file=sprintf("%s/randomForest_METABOLITE_ObeseNML.%s_%s.%s.cv.txt", out_dir, mvar_level, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
			## END BLOCK TO COMMENT ##

			collated.importance <- read.table(sprintf("%s/randomForest_METABOLITE_ObeseNML.%s_%s.%s.importance.txt", out_dir, mvar_level, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
			collated.cv <- read.table(sprintf("%s/randomForest_METABOLITE_ObeseNML.%s_%s.%s.cv.txt", out_dir, mvar_level, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
			importance.mean <- rowMeans(collated.importance)
			importance.sd <- unlist(apply(collated.importance, 1, sd))
			cv.mean <- rowMeans(collated.cv)
			cv.sd <- unlist(apply(collated.cv, 1, sd))
			res.mean <- rbind(res.mean, importance.mean)
			res.sd <- rbind(res.sd, importance.sd)
			inds <- order(importance.mean, decreasing=T)
			inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.001)))] # edit as appropriate
			res.sig[[length(res.sig)+1]] <- names(importance.mean[inds])
			write.table(melt(importance.mean[inds]), file=sprintf("%s/randomForest_METABOLITE_ObeseNML.%s_%s.%s.features.txt", out_dir, mvar_level, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

			## after running for the first time, COMMENT OUT THIS BLOCK ##
			# using a sparse model with N predictors
			sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, ntree=10000, importance=T)
			save(sparseRF, file=sprintf("%s/randomForest_METABOLITE_ObeseNML.%s_%s.%s.model", out_dir, mvar_level, mlevel, st))
			load(sprintf("%s/randomForest_METABOLITE_ObeseNML.%s_%s.%s.model", out_dir, mvar_level, mlevel, st))
			# ROC and AUC of final sparseRF model
			pred <- predict(sparseRF, type="prob")
			pred2 <- prediction(pred[,2], ordered(response))
			perf <- performance(pred2, "tpr", "fpr")
			perf.auc <- performance(pred2, "auc")
			pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
			confusion_matrix <- table(pred_df[, c("true", "predicted")])
			accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
			vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
			mccvalue <- mcc(vec.pred, vec.true)
			plot(perf, main=sprintf("ROC %s %s %s (sparseRF final model ObeseNML baseline)", mvar_level, mlevel, st)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g, accuracy=%.2f%%, MCC=%.4g", unlist(perf.auc@y.values), accuracy, mccvalue))			
			## END BLOCK TO COMMENT ##
			
			# plotting - per-group sparse model
			df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
			colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
			print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s", mvar_level, mlevel, st)))
			# plotting - per-group variables
			df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
			df$metabolite_name <- as.character(df$OTU)
			res <- read.table(sprintf("%s/METABOLITE_linear_regression.%s.%s.txt", out_dir, st, mlevel), header=T, sep="\t", as.is=T, quote="")
			res <- subset(res, coefficient == paste("Cohort", mvar_level, sep="")); rownames(res) <- res$metabolite
			df$Z <- res[as.character(df$OTU), "Z"]; df$padj <- res[as.character(df$OTU), "padj"]
			if (mlevel == "BIOCHEMICAL") {
				df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
				df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
				df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g; %s; %s)", df$OTU, df$Z, df$padj, df$subpathway, df$superpathway)
			} else if (mlevel == "SUB.PATHWAY") {
				df$subpathway <- df$metabolite_name
				df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
				df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g)", df$OTU, df$Z, df$padj)
			} else {
				df$superpathway <- df$metabolite_name
				df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g)", df$OTU, df$Z, df$padj)
			}
			df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
			df$sig <- ifelse(df$padj < 0.1, "sig", "ns"); df$sig[which(is.na(df$sig))] <- "ns"
			p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string, color=sig), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory %s", mvar_level, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
			print(p)
			# shading rectangles of importance
			df.rect <- df
			df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
			df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
			p <- ggplot(df.rect, aes(x=x, y=OTU, color=sig, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s explanatory %s", mvar_level, mlevel, st)) + scale_fill_gradient(low="white", high="black") + scale_color_manual(values=cols.sig)
			print(p)
			# violin plots of relative abundance
			agg.melt <- agg.melt.stored
			agg.melt$Cohort <- mapping.sel[match(agg.melt$SampleID, mapping.sel$StudyID), "Cohort"]
			agg.melt <- subset(agg.melt, metabolite %in% rownames(df) & Cohort %in% c(mvar_level, "ObeseNML"))
			agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
			if (nrow(agg.melt) > 0) {
				p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mvar_level, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
				agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
				p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Cohort)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mvar_level, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
				agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
				agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
				p <- ggplot(agg.melt2, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mvar_level, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
				p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Cohort)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mvar_level, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
			}
		}
		# plotting - all groups importance values
		rownames(res.mean) <- setdiff(levels(mapping.sel[, "Cohort"]), c("NMLBMI", "ObeseNML")); rownames(res.sd) <- rownames(res.mean)
		df <- res.mean[, unique(unlist(res.sig))]; df[which(df<0.001, arr.ind=T)] <- 0
		df.sig <- matrix("", nrow=nrow(df), ncol=ncol(df)); rownames(df.sig) <- rownames(df); colnames(df.sig) <- colnames(df)
		for (i in 1:length(res.sig)) {
			df.sig[i, res.sig[[i]]] <- "*"
		}
		heatmap.2(df, Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(10,10), cexCol=0.6, cexRow=0.6, cellnote=df.sig, notecex=0.8, notecol="red", main=sprintf("RF importance values (%s, %s)", st, mlevel))
	}
}

## randomForest classification of Cohort (multiclass); using METABOLITE data [plasma, fecal]
set.seed(112817)	
# MAIN LOOP for random forest (through metabolite levels)
for (st in c("fecal", "plasma")) {
	for (mlevel in metabolite_levels) {
		data <- df.metabolon[[st]][[mlevel]]
		mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% rownames(data))
		data.sel <- as.data.frame(data)
		response <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Cohort"]); names(response) <- rownames(data.sel)
		# add BMI, Age, Sex as covariates
		data.sel$BMI <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "BMI"])
		data.sel$Age <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Age"])
		data.sel$Sex <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Sex"])
		agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Sex")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
		
		## after running for the first time, COMMENT OUT THIS BLOCK ##
		num_iter <- 100
		ncores <- 20
		out <- mclapply(1:num_iter, function (dummy) {
				importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
			}, mc.cores=ncores )
		collated.importance <- do.call(cbind, out)
		out <- mclapply(1:num_iter, function (dummy) {
				rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
			}, mc.cores=ncores )
		collated.cv <- do.call(cbind, out)

		write.table(collated.importance, file=sprintf("%s/randomForest_METABOLITE.%s_%s.%s.importance.txt", out_dir, "multiclass", mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
		write.table(collated.cv, file=sprintf("%s/randomForest_METABOLITE.%s_%s.%s.cv.txt", out_dir, "multiclass", mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
		## END BLOCK TO COMMENT ##

		collated.importance <- read.table(sprintf("%s/randomForest_METABOLITE.%s_%s.%s.importance.txt", out_dir, "multiclass", mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
		collated.cv <- read.table(sprintf("%s/randomForest_METABOLITE.%s_%s.%s.cv.txt", out_dir, "multiclass", mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
		importance.mean <- rowMeans(collated.importance)
		importance.sd <- unlist(apply(collated.importance, 1, sd))
		cv.mean <- rowMeans(collated.cv)
		cv.sd <- unlist(apply(collated.cv, 1, sd))
		inds <- order(importance.mean, decreasing=T)
		inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.001)))] # edit as appropriate
		write.table(melt(importance.mean[inds]), file=sprintf("%s/randomForest_METABOLITE.%s_%s.%s.features.txt", out_dir, "multiclass", mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

		## after running for the first time, COMMENT OUT THIS BLOCK ##
		# using a sparse model with N predictors
		sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
		save(sparseRF, file=sprintf("%s/randomForest_METABOLITE.%s_%s.%s.model", out_dir, "multiclass", mlevel, st))
		load(sprintf("%s/randomForest_METABOLITE.%s_%s.%s.model", out_dir, "multiclass", mlevel, st))
		# accuracy of final sparseRF model
		pred <- predict(sparseRF, type="prob")
		pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
		pred_df_out <- merge(pred_df, data.sel, by="row.names")
		write.table(pred_df_out, file=sprintf("%s/randomForest_METABOLITE.%s_%s.%s.predictions.txt", out_dir, "multiclass", mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
		confusion_matrix <- table(pred_df[, c("true", "predicted")])
		class_errors <- unlist(lapply(levels(mapping.sel$Cohort), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$Cohort)
		accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
		vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
		mccvalue <- mcc(vec.pred, vec.true)
		df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
		p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", "multiclass", mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
		print(p)
		write.table(confusion_matrix, file=sprintf("%s/randomForest_METABOLITE.%s_%s.%s.confusion_matrix.txt", out_dir,  "multiclass", mlevel, st), quote=F, sep="\t", row.names=T, col.names=T)
		## END BLOCK TO COMMENT ##
		
		# plotting - per-group sparse model
		df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
		colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
		print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s", "multiclass", mlevel, st)))
		# plotting - per-group variables
		df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
		df$metabolite_name <- as.character(df$OTU)
		res <- read.table(sprintf("%s/METABOLITE_linear_regression.%s.%s.txt", out_dir, st, mlevel), header=T, sep="\t", as.is=T, quote="")
		res <- res[grep("Cohort", res$coefficient),]
		res <- aggregate(padj ~ metabolite, res, min); rownames(res) <- res$metabolite
		df$padj <- res[as.character(df$OTU), "padj"]
		df$sig <- ifelse(df$padj < 0.1, "sig", "ns"); df$sig[which(is.na(df$sig))] <- "ns"
		if (mlevel == "BIOCHEMICAL") {
			df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
			df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
			df$OTU_string <- sprintf("%s (min padj=%.4g; %s; %s)", df$OTU, df$padj, df$subpathway, df$superpathway)
		} else if (mlevel == "SUB.PATHWAY") {
			df$subpathway <- df$metabolite_name
			df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
			df$OTU_string <- sprintf("%s (min padj=%.4g)", df$OTU, df$padj)
		} else {
			df$superpathway <- df$metabolite_name
			df$OTU_string <- sprintf("%s (min padj=%.4g)", df$OTU, df$padj)
		}
		df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
		p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string, color=sig), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory %s", "multiclass", mlevel, st)) + scale_fill_manual(values=cols.superpathway) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
		print(p) 
		# shading rectangles of importance values
		df.rect <- df
		df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
		df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
		p <- ggplot(df.rect, aes(x=x, y=OTU, color=sig, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s explanatory %s", "multiclass", mlevel, st)) + scale_fill_gradient(low="white", high="black") + scale_color_manual(values=cols.sig)
		print(p)
		# violin plots of relative abundance
		agg.melt <- agg.melt.stored
		agg.melt$Cohort <- mapping.sel[match(agg.melt$SampleID, mapping.sel$StudyID), "Cohort"]
		agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$Cohort))
		agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
		agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
		agg.melt$L099 <- ifelse(agg.melt$SampleID=="L099", "L099", "not L099")
		p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=L099)) + geom_violin() + geom_point(size=1) + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mlevel, st)) + coord_flip() + scale_color_manual(values=c("red", "grey"))
		print(p)
		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Cohort)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
		agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
		p <- ggplot(agg.melt2, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Cohort)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mlevel, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
		# violin plots stratified by prediction/truth
		for (met in levels(agg.melt$metabolite)) {
			tmp <- subset(agg.melt, metabolite==met)
			p <- ggplot(tmp, aes(x=Cohort, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s)", met, mlevel, st)) + scale_fill_manual(values=cols.cohort)
			print(p)			
			# color scheme - manual
			p <- ggplot(tmp, aes(x=Cohort, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s)", met, mlevel, st)) + scale_fill_manual(values=cols.cohort)
			print(p)
			p <- ggplot(tmp, aes(x=Cohort, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s)", met, mlevel, st)) + scale_color_manual(values=cols.cohort)
			print(p)
			
		}
	}
}


## randomForest classification of Cohort (one-vs-one setup, NASH vs NAFL baseline); using METABOLITE data [plasma, fecal]
set.seed(091818)	
# MAIN LOOP for random forest (through metabolite levels)
for (st in c("fecal", "plasma")) {
	for (mlevel in metabolite_levels) {
		data <- df.metabolon[[st]][[mlevel]]
		mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% rownames(data))
		# random forest classification
		res.mean <- {}; res.sd <- {}; res.sig <- list()
		for (mvar_level in setdiff(levels(mapping.sel[, "Cohort"]), c("NMLBMI", "ObeseNML", "NAFLD"))) {
			sel <- intersect(rownames(data), unique(subset(sample_data(ps), Cohort %in% c(mvar_level, "NAFLD"))$StudyID))
			data.sel <- as.data.frame(data[sel,])
			response <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Cohort"]); names(response) <- sel
			# add BMI, Age, Sex as covariates
			data.sel$BMI <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "BMI"])
			data.sel$Age <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Age"])
			data.sel$Sex <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Sex"])
			agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Sex")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
			
			## after running for the first time, COMMENT OUT THIS BLOCK ##
			num_iter <- 100
			ncores <- 20
			out <- mclapply(1:num_iter, function (dummy) {
					importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
				}, mc.cores=ncores )
			collated.importance <- do.call(cbind, out)
			out <- mclapply(1:num_iter, function (dummy) {
					rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
				}, mc.cores=ncores )
			collated.cv <- do.call(cbind, out)

			write.table(collated.importance, file=sprintf("%s/randomForest_METABOLITE_NAFLD.%s_%s.%s.importance.txt", out_dir, mvar_level, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
			write.table(collated.cv, file=sprintf("%s/randomForest_METABOLITE_NAFLD.%s_%s.%s.cv.txt", out_dir, mvar_level, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
			## END BLOCK TO COMMENT ##

			collated.importance <- read.table(sprintf("%s/randomForest_METABOLITE_NAFLD.%s_%s.%s.importance.txt", out_dir, mvar_level, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
			collated.cv <- read.table(sprintf("%s/randomForest_METABOLITE_NAFLD.%s_%s.%s.cv.txt", out_dir, mvar_level, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
			importance.mean <- rowMeans(collated.importance)
			importance.sd <- unlist(apply(collated.importance, 1, sd))
			cv.mean <- rowMeans(collated.cv)
			cv.sd <- unlist(apply(collated.cv, 1, sd))
			res.mean <- rbind(res.mean, importance.mean)
			res.sd <- rbind(res.sd, importance.sd)
			inds <- order(importance.mean, decreasing=T)
			inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.001)))] # edit as appropriate
			res.sig[[length(res.sig)+1]] <- names(importance.mean[inds])
			write.table(melt(importance.mean[inds]), file=sprintf("%s/randomForest_METABOLITE_NAFLD.%s_%s.%s.features.txt", out_dir, mvar_level, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

			## after running for the first time, COMMENT OUT THIS BLOCK ##
			# using a sparse model with N predictors
			sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, ntree=10000, importance=T)
			save(sparseRF, file=sprintf("%s/randomForest_METABOLITE_NAFLD.%s_%s.%s.model", out_dir, mvar_level, mlevel, st))
			load(sprintf("%s/randomForest_METABOLITE_NAFLD.%s_%s.%s.model", out_dir, mvar_level, mlevel, st))
			# ROC and AUC of final sparseRF model
			pred <- predict(sparseRF, type="prob")
			pred2 <- prediction(pred[,2], ordered(response))
			perf <- performance(pred2, "tpr", "fpr")
			perf.auc <- performance(pred2, "auc")
			pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
			confusion_matrix <- table(pred_df[, c("true", "predicted")])
			accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
			vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
			mccvalue <- mcc(vec.pred, vec.true)
			plot(perf, main=sprintf("ROC %s %s %s (sparseRF final model NAFLD baseline)", mvar_level, mlevel, st)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g, accuracy=%.2f%%, MCC=%.4g", unlist(perf.auc@y.values), accuracy, mccvalue))			
			## END BLOCK TO COMMENT ##
			
			# plotting - per-group sparse model
			df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
			colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
			print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s", mvar_level, mlevel, st)))
			# plotting - per-group variables
			df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
			df$metabolite_name <- as.character(df$OTU)
			res <- read.table(sprintf("%s/METABOLITE_linear_regression.%s.%s.txt", out_dir st, mlevel), header=T, sep="\t", as.is=T, quote="")
			res <- subset(res, coefficient == paste("Cohort", mvar_level, sep="")); rownames(res) <- res$metabolite
			df$Z <- res[as.character(df$OTU), "Z"]; df$padj <- res[as.character(df$OTU), "padj"]
			if (mlevel == "BIOCHEMICAL") {
				df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
				df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
				df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g; %s; %s)", df$OTU, df$Z, df$padj, df$subpathway, df$superpathway)
			} else if (mlevel == "SUB.PATHWAY") {
				df$subpathway <- df$metabolite_name
				df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
				df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g)", df$OTU, df$Z, df$padj)
			} else {
				df$superpathway <- df$metabolite_name
				df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g)", df$OTU, df$Z, df$padj)
			}
			df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
			df$sig <- ifelse(df$padj < 0.1, "sig", "ns"); df$sig[which(is.na(df$sig))] <- "ns"
			p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string, color=sig), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory %s", mvar_level, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
			print(p)
			# shading rectangles of importance
			df.rect <- df
			df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
			df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
			p <- ggplot(df.rect, aes(x=x, y=OTU, color=sig, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s explanatory %s", mvar_level, mlevel, st)) + scale_fill_gradient(low="white", high="black") + scale_color_manual(values=cols.sig)
			print(p)
			# violin plots of relative abundance
			agg.melt <- agg.melt.stored
			agg.melt$Cohort <- mapping.sel[match(agg.melt$SampleID, mapping.sel$StudyID), "Cohort"]
			agg.melt <- subset(agg.melt, metabolite %in% rownames(df) & Cohort %in% c(mvar_level, "NAFLD"))
			agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
			if (nrow(agg.melt) > 0) {
				p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mvar_level, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
				agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
				p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Cohort)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mvar_level, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
				agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
				agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
				p <- ggplot(agg.melt2, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mvar_level, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
				p <- ggplot(agg.melt2, aes(x=metabolite, y=value, color=Cohort)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", mvar_level, st)) + coord_flip() + scale_color_manual(values=cols.cohort)
				print(p)
			}
		}
	}
}



## PLS-DA classification of Cohort (multiclass); using plasma METABOLITE data
set.seed(112817)
nrepeats <- 500
ncomps <- 5

# no looping for PLS-DA/sPLS-DA as need to manually tune parameters
# plasma BIOCHEMICAL
st <- "plasma"; mlevel <- "BIOCHEMICAL"
data <- df.metabolon[[st]][[mlevel]]
mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% rownames(data))
data.sel <- as.data.frame(data)
response <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Cohort"]); names(response) <- rownames(data.sel)
data.sel$BMI <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "BMI"])
data.sel$Age <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Age"])
data.sel$Sex <- as.numeric(unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Sex"]))
agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Sex")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
data.sel <- as.matrix(data.sel)
# PLS-DA
plsda.res <- plsda(data.sel, response, ncomp = ncomps) # where ncomp is the number of components wanted
perf.plsda <- perf(plsda.res, validation = "Mfold", folds = 5, progressBar = FALSE, auc = TRUE, nrepeat = nrepeats)
plot(perf.plsda, sd = TRUE, legend.position = "horizontal")
plotIndiv(plsda.res, comp = 1:2, group = plsda.res$Y, ind.names = FALSE, ellipse = TRUE, legend = TRUE, title=sprintf("PLS-DA (%s, %s)", st, mlevel))
plotIndiv(plsda.res, comp = 1:2, group = plsda.res$Y, ind.names = FALSE, ellipse = TRUE, star = TRUE, legend = TRUE, title=sprintf("PLS-DA (%s, %s)", st, mlevel))
background <- background.predict(plsda.res, comp.predicted=2, dist = "max.dist") 
plotIndiv(plsda.res, comp = 1:2, group = response, ind.names = FALSE, title = sprintf("PLS-DA (%s, %s)", st, mlevel), legend = TRUE,  background = background)
# sPLS-DA
list.keepX <- c(1:10,  seq(from=20, to=100, by=10))
tune.splsda.multiclass <- tune.splsda(data.sel, response, ncomp = ncomps, validation = 'Mfold', folds = 5, progressBar = TRUE, dist = 'max.dist', measure = "BER", test.keepX = list.keepX, nrepeat = nrepeats, cpus = 16)
error <- tune.splsda.multiclass$error.rate
ncomp <- tune.splsda.multiclass$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp
select.keepX <- tune.splsda.multiclass$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX
plot(tune.splsda.multiclass)
splsda.final <- splsda(data.sel, response, ncomp = ncomp, keepX = select.keepX)
plotIndiv(splsda.final, ind.names = FALSE, legend=TRUE, ellipse = TRUE, title=sprintf("SPLS-DA final result (%s, %s)", st, mlevel))
plotIndiv(splsda.final, ind.names = FALSE, legend=TRUE, ellipse = TRUE, star = TRUE, title=sprintf("SPLS-DA final result (%s, %s)", st, mlevel))
background <- background.predict(splsda.final, comp.predicted=2, dist = "max.dist") 
plotIndiv(splsda.final, comp = 1:2, group = response, ind.names = FALSE, title = sprintf("SPLS-DA (%s, %s)", st, mlevel), legend = TRUE,  background = background)
perf.final <- perf(splsda.final, validation = "Mfold", folds = 5, dist = 'max.dist', nrepeat = nrepeats, progressBar = FALSE) 
plot(perf.final)
# compute mean+SD for performance metrics (using ncomp components as this is model with the selected number of components)
cmlist <- lapply(1:nrepeats, function(repi) get.confusion_matrix(truth=response, predicted=perf.final$class[["max.dist"]][,repi,ncomp]))
cmarr <- abind(cmlist, along=3)
cm.mean <- apply(cmarr, c(1,2), mean); cm.sd <- apply(cmarr, c(1,2), sd)
accuracy <- 100*apply(cmarr, 3, function(x) sum(diag(x))/sum(x))
mccvalues <- {}
for (i in 1:nrepeats) {
	predicted <- factor(perf.final$class[["max.dist"]][,i,ncomp][names(response)], levels=levels(response))
	mccvalue <- mcc(as.numeric(predicted)-1, as.numeric(response)-1)
	mccvalues <- c(mccvalues, mccvalue)
}
metrics <- data.frame(accuracy=accuracy, mcc=mccvalues)
cmtab <- matrix(paste(formatC(cm.mean, format="f", digits=3), formatC(cm.sd, format="f", digits=3), sep=" +/- "), nrow=nrow(cm.mean), dimnames=dimnames(cm.mean))
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix mean+/-SD (%s, %s)", st, mlevel)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(cmtab), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)
p <- ggplot(metrics, aes(accuracy)) + geom_histogram() + theme_classic() + ggtitle(sprintf("Accuracy (%s, %s) mean=%.4g SD=%.4g", st, mlevel, mean(metrics$accuracy), sd(metrics$accuracy))) + geom_vline(xintercept=mean(metrics$accuracy), col="red")
print(p)
p <- ggplot(metrics, aes(mcc)) + geom_histogram() + theme_classic() + ggtitle(sprintf("MCC (%s, %s) mean=%.4g SD=%.4g", st, mlevel, mean(metrics$mcc), sd(metrics$mcc))) + geom_vline(xintercept=mean(metrics$mcc), col="red")
print(p)
# arrow plots, loadings, feature stability
plotArrow(splsda.final, legend=T)
for (i in 1:ncomp) {
	plotLoadings(splsda.final, comp = i, title = sprintf("Loadings on comp %d (%s, %s)", i, st, mlevel), contrib = 'max', method = 'mean')
}
for (i in 1:ncomp) {
	inds <- match(selectVar(splsda.final, comp = i)$name, names(perf.final$features$stable[[i]]))
	freq <- as.numeric(perf.final$features$stable[[i]][inds])
	df <- data.frame(selectVar(splsda.final, comp = i)$value, freq); df <- df[order(df$freq, decreasing=F),]; df$feature <- factor(rownames(df), levels=rownames(df)) # feature stability
	p <- ggplot(df, aes(x=feature, y=freq)) + geom_bar(stat="identity") + coord_flip() + theme_classic() + ggtitle(sprintf("Feature stability (%s, %s)", st, mlevel))
	print(p)
}



###############################################################################################################
### METABOLITES (Biopsy subgroup)


## read in biopsy data
biopsy <- read.table(sprintf("%s/liver_biopsy.txt", data_dir), header=T, row.names=1, as.is=T, sep="\t")
colnames(biopsy) <- c("NAS_score", "Diagnosis", "Steatosis", "LobularInflammation", "PortalInflammation", "Ballooning", "Fibrosis", "Cirrhosis", "Metabolon")
biopsy$Diagnosis2 <- factor(ifelse(biopsy$NAS_score >= 5, "DefiniteNASH", "BorderlineNASH"))
biopsy$Fibrosis2 <- factor(ifelse(biopsy$Fibrosis >= 3, "Advanced", "Mild"), levels=c("Mild", "Advanced"))
biopsy$Steatosis2 <- factor(pmax(biopsy$Steatosis, 1)) # convert 0 -> 1
biopsy$LobularInflammation2 <- factor(pmax(biopsy$LobularInflammation, 1)) # convert 0 -> 1
biopsy$Ballooning2 <- factor(biopsy$Ballooning)
biopsy_variables <- c("Diagnosis2", "Steatosis2", "Fibrosis2", "LobularInflammation2", "Ballooning2")

## linear regression of METABOLITE data ~ biopsy_variable + BMI
for (st in c("fecal", "plasma")) {
	for (mlevel in metabolite_levels) {
		data <- df.metabolon[[st]][[mlevel]]
		sel <- intersect(rownames(biopsy), rownames(data))
		data <- data[sel,]
		mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% sel)
		outfull <- list()
		for (biopsy_variable in biopsy_variables) {
			out <- mclapply(colnames(data), function(metabolite) {
				df <- melt(data[,metabolite])	
				df[, "Cohort"] <- unlist(sample_data(ps)[match(rownames(df), sample_data(ps)$StudyID), "Cohort"])
				df[, "BMI"] <- unlist(sample_data(ps)[match(rownames(df), sample_data(ps)$StudyID), "BMI"])
				df[, "Age"] <- unlist(sample_data(ps)[match(rownames(df), sample_data(ps)$StudyID), "Age"])
				df[, "Sex"] <- unlist(sample_data(ps)[match(rownames(df), sample_data(ps)$StudyID), "Sex"])
				df <- merge(df, biopsy[, biopsy_variables], by="row.names")
				tt <- try(m <- lm(as.formula(sprintf("value ~ %s + BMI + Age + Sex", biopsy_variable)), data = df), silent=T)
				if (class(tt) == "lm") {
					coef <- summary(m)$coefficients # rows are [(Intercept), comparisons], columns are [Estimate, SE, Z, pval]
					tmp <- cbind(biopsy_variable, metabolite, "Linear regression", rownames(coef), coef)
				}
				print(sprintf("finished %s", metabolite))
				tmp
			}, mc.cores=16)
			outfull <- c(outfull, out)
		}
		res <- as.data.frame(do.call(rbind, outfull))
		colnames(res) <- c("biopsy_variable", "metabolite", "method", "coefficient", "Estimate", "SE", "Z", "pvalue")
		res <- subset(res, !coefficient %in% c("(Intercept)", "Age", "BMI", "SexM")) # remove intercept/control terms before FDR correction
		res$pvalue <- as.numeric(as.character(res$pvalue))
		## p-value adjust
		res$padj <- p.adjust(res$pvalue, method="fdr")
		res <- res[order(res$padj, decreasing=F),]
		resSig <- subset(res, padj<0.05)
		write.table(res, file=sprintf("%s/BIOPSY_METABOLITE_linear_regression.%s.%s.txt", out_dir, st, mlevel), quote=F, sep="\t", row.names=F, col.names=T)
	}
}

## ordination (t-SNE, PCA) and PERMANOVA
set.seed(120617)
for (st in c("fecal", "plasma")) {
	mlevel <- "BIOCHEMICAL"
	data <- df.metabolon[[st]][[mlevel]]
	sel <- intersect(rownames(biopsy), rownames(data))
	data <- data[sel,]
	mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% sel)
	agg <- aggregate(Cohort ~ StudyID, mapping.sel, unique); agg2 <- aggregate(cbind(BMI, Age) ~ StudyID, mapping.sel, unique); agg3 <- aggregate(Sex ~ StudyID, mapping.sel, unique)
	mapping.sel <- merge(merge(agg, agg2, by="StudyID"), agg3, by="StudyID"); rownames(mapping.sel) <- mapping.sel$StudyID
	mapping.sel <- merge(mapping.sel, biopsy[, biopsy_variables], by="row.names"); rownames(mapping.sel) <- mapping.sel$StudyID
	# PCA
	pca <- prcomp(data, center=F, scale=F)
	eigs <- pca$sdev^2
	pvar <- 100*(eigs / sum(eigs))
	df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], SampleID=rownames(pca$x))
	for (mvar in biopsy_variables) {
		df[, mvar] <- mapping.sel[rownames(df), mvar]
		if (class(mapping.sel[, mvar])[1] == "numeric") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s PCoA (Euclidean distance)", st, mlevel)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_gradient(low="grey", high="red")
		}
		else {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s PCoA (Euclidean distance)", st, mlevel)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2]))
		}
		print(p)
	}
	# PERMANOVA
	res <- adonis(as.formula(sprintf("data ~ BMI + Age + Sex + %s", paste(biopsy_variables, collapse=" + "))), data=mapping.sel, permutations=999, method="eu")
	sink(sprintf("%s/BIOPSY_METABOLITE_PERMANOVA.%s.%s.txt", out_dir, st, mlevel))
	print(res)
	# demographics
	for (bvar in c("NAS_score", biopsy_variables)) {
		tab <- table(biopsy[sel,bvar])
		print(bvar)
		print(tab)
	}
	sink()
	# t-SNE
	tsne.out <- Rtsne(data, perplexity=5)
	df <- data.frame(PC1 = tsne.out$Y[,1], PC2 = tsne.out$Y[,2], SampleID=rownames(data)); rownames(df) <- df$SampleID
	for (mvar in biopsy_variables) {
		df[, mvar] <- mapping.sel[rownames(df), mvar]
		if (class(mapping.sel[, mvar])[1] == "numeric") {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s tSNE", st, mlevel)) + scale_color_gradient(low="grey", high="red")
		}
		else {
			p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("%s %s tSNE", st, mlevel))
		}
		print(p)
	}
}

## randomForest classification of biopsy variables (multiclass); using METABOLITE data [plasma, fecal]
set.seed(120717)	
# MAIN LOOP for random forest (through metabolite levels)
for (st in c("fecal", "plasma")) {
	for (mlevel in metabolite_levels[1:2]) {
		data <- df.metabolon[[st]][[mlevel]]
		sel <- intersect(rownames(biopsy), rownames(data))
		data <- data[sel,]
		mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% rownames(data))
		data.sel <- as.data.frame(data)
		for (biopsy_variable in biopsy_variables) {
			response <- biopsy[sel, biopsy_variable]; names(response) <- sel
			# add BMI, Age, Sex as covariates
			data.sel$BMI <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "BMI"])
			data.sel$Age <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Age"])
			data.sel$Sex <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Sex"])
			agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Sex")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
		
			## after running for the first time, COMMENT OUT THIS BLOCK ##
			num_iter <- 100
			ncores <- 20
			out <- mclapply(1:num_iter, function (dummy) {
					importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
				}, mc.cores=ncores )
			collated.importance <- do.call(cbind, out)
			out <- mclapply(1:num_iter, function (dummy) {
					rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
				}, mc.cores=ncores )
			collated.cv <- do.call(cbind, out)

			write.table(collated.importance, file=sprintf("%s/randomForest_BIOPSY_METABOLITE.%s_%s.%s.importance.txt", out_dir, biopsy_variable, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
			write.table(collated.cv, file=sprintf("%s/randomForest_BIOPSY_METABOLITE.%s_%s.%s.cv.txt", out_dir, biopsy_variable, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)
			## END BLOCK TO COMMENT ##

			collated.importance <- read.table(sprintf("%s/randomForest_BIOPSY_METABOLITE.%s_%s.%s.importance.txt", out_dir, biopsy_variable, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
			collated.cv <- read.table(sprintf("%s/randomForest_BIOPSY_METABOLITE.%s_%s.%s.cv.txt", out_dir, biopsy_variable, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
			importance.mean <- rowMeans(collated.importance)
			importance.sd <- unlist(apply(collated.importance, 1, sd))
			cv.mean <- rowMeans(collated.cv)
			cv.sd <- unlist(apply(collated.cv, 1, sd))
			inds <- order(importance.mean, decreasing=T)
			inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.001)))] # edit as appropriate
			write.table(melt(importance.mean[inds]), file=sprintf("%s/randomForest_BIOPSY_METABOLITE.%s_%s.%s.features.txt", out_dir, biopsy_variable, mlevel, st), quote=F, sep="\t", row.names=T, col.names=F)

			## after running for the first time, COMMENT OUT THIS BLOCK ##
			# using a sparse model with N predictors
			sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, ntree=10000, importance=T)
			save(sparseRF, file=sprintf("%s/randomForest_BIOPSY_METABOLITE.%s_%s.%s.model", out_dir, biopsy_variable, mlevel, st))
			load(sprintf("%s/randomForest_BIOPSY_METABOLITE.%s_%s.%s.model", out_dir, biopsy_variable, mlevel, st))
			# accuracy of final sparseRF model
			pred <- predict(sparseRF, type="prob")
			pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
			pred_df_out <- merge(pred_df, data.sel, by="row.names")
			write.table(pred_df_out, file=sprintf("%s/randomForest_BIOPSY_METABOLITE.%s_%s.%s.predictions.txt", out_dir, biopsy_variable, mlevel, st), quote=F, sep="\t", row.names=F, col.names=T)
			confusion_matrix <- table(pred_df[, c("true", "predicted")])
			class_errors <- unlist(lapply(levels(biopsy[,biopsy_variable]), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(biopsy[,biopsy_variable])
			accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
			vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
			mccvalue <- mcc(vec.pred, vec.true)
			df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
			p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s, %s) (accuracy = %.2f%%, MCC = %.4f)", "multiclass", mlevel, st, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
			print(p)
			write.table(confusion_matrix, file=sprintf("%s/randomForest_BIOPSY_METABOLITE.%s_%s.%s.confusion_matrix.txt", out_dir, biopsy_variable, mlevel, st), quote=F, sep="\t", row.names=T, col.names=T)
			# ROC analysis
			if (nlevels(response)==2) {
				pred <- predict(sparseRF, type="prob")
				pred2 <- prediction(pred[,2], ordered(response))
				perf <- performance(pred2, "tpr", "fpr")
				perf.auc <- performance(pred2, "auc")
				print(plot(perf, main=sprintf("ROC %s %s %s (sparseRF final model)", biopsy_variable, mlevel, st)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", unlist(perf.auc@y.values))))
			}
			## END BLOCK TO COMMENT ##
		
			# plotting - per-group sparse model
			df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
			colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
			print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s %s", "multiclass", mlevel, st)))
			# plotting - per-group variables
			df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
			df$metabolite_name <- as.character(df$OTU)
			res <- read.table(sprintf("%s/BIOPSY_METABOLITE_linear_regression.%s.%s.txt", out_dir, st, mlevel), header=T, sep="\t", as.is=T, quote=""); colnames(res)[1] <- "BiopsyVariable"
			res <- subset(res, BiopsyVariable==biopsy_variable)
			res <- ddply(res, .(metabolite), function(x) {
				x[which.min(x[, "padj"]),]
			})
			rownames(res) <- res$metabolite
			df$Z <- res[as.character(df$OTU), "Z"]; df$padj <- res[as.character(df$OTU), "padj"]
			if (mlevel == "BIOCHEMICAL") {
				df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
				df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
				df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g; %s; %s)", df$OTU, df$Z, df$padj, df$subpathway, df$superpathway)
			} else if (mlevel == "SUB.PATHWAY") {
				df$subpathway <- df$metabolite_name
				df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
				df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g)", df$OTU, df$Z, df$padj)
			} else {
				df$superpathway <- df$metabolite_name
				df$OTU_string <- sprintf("%s (Z=%.4g, padj=%.4g)", df$OTU, df$Z, df$padj)
			}
			df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
			df$sig <- ifelse(df$padj < 0.1, "sig", "ns"); df$sig[which(is.na(df$sig))] <- "ns"
			p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string, color=sig), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory %s", biopsy_variable, mlevel, st)) + scale_fill_manual(values=cols.superpathway) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
			print(p)
			# violin plots of relative abundance
			agg.melt <- agg.melt.stored
			agg.melt$Truth <- biopsy[agg.melt$SampleID, biopsy_variable]
			agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$Truth))
			agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
			agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
			p <- ggplot(agg.melt, aes(x=Truth, y=value, color=Truth)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s, %s)", biopsy_variable, mlevel, st)) + coord_flip()
			print(p)
			# violin plots stratified by prediction/truth
			for (met in levels(agg.melt$metabolite)) {
				tmp <- subset(agg.melt, metabolite==met)
				# color scheme - Set1
				p <- ggplot(tmp, aes(x=Truth, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s)", met, mlevel, st)) + scale_fill_brewer(palette="Set1")
				print(p)
				p <- ggplot(tmp, aes(x=Truth, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s, %s)", met, mlevel, st)) + scale_color_brewer(palette="Set1")
				print(p)
			}
		}
	}
}

# heatmap of all biopsy classification models
for (cutoff in c(0.001, 0.002, 0.005)) {
	for (st in c("fecal", "plasma")) {
		for (mlevel in metabolite_levels[1:2]) {
			res <- {}
			for (biopsy_variable in biopsy_variables) {
				collated.importance <- read.table(sprintf("%s/randomForest_BIOPSY_METABOLITE.%s_%s.%s.importance.txt", out_dir, biopsy_variable, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
				collated.cv <- read.table(sprintf("%s/randomForest_BIOPSY_METABOLITE.%s_%s.%s.cv.txt", out_dir, biopsy_variable, mlevel, st), header=F, as.is=T, sep="\t", row.names=1)
				importance.mean <- rowMeans(collated.importance)
				importance.sd <- unlist(apply(collated.importance, 1, sd))
				cv.mean <- rowMeans(collated.cv)
				cv.sd <- unlist(apply(collated.cv, 1, sd))
				inds <- order(importance.mean, decreasing=T)
				inds <- inds[1:min(30, length(which(importance.mean[inds] > cutoff)))]
				tmp <- data.frame(metabolite=names(importance.mean[inds]), importance=importance.mean[inds], mvar=biopsy_variable)
				res <- rbind(res, tmp)
			}
			df <- dcast(res, metabolite ~ mvar, value.var="importance")
			rownames(df) <- df$metabolite; df <- as.matrix(df[,-1])
#			df <- log2(df)
			df[which(is.na(df), arr.ind=T)] <- 0 # NA -> 0
			res <- read.table(sprintf("%s/BIOPSY_METABOLITE_linear_regression.%s.%s.txt", out_dir, st, mlevel), header=T, sep="\t", as.is=T, quote=""); colnames(res)[1] <- "BiopsyVariable"
			resSig <- subset(res, padj<0.1)
			df.sig <- matrix("", nrow=nrow(df), ncol=ncol(df)); rownames(df.sig) <- rownames(df); colnames(df.sig) <- colnames(df)
			if (nrow(resSig)>0) {
				for (i in 1:length(resSig)) {
					df.sig[resSig[i,"metabolite"], resSig[i,"BiopsyVariable"]] <- "*"
				}
			}
			for (met in rownames(df)) {
				for (bvar in colnames(df)) {
					tmp <- sign(mean(res[which(met == res$metabolite & bvar == res$BiopsyVariable), "Estimate"]))
					tmp <- ifelse(is.na(tmp), 1, tmp)
					df[met, bvar] <- tmp*df[met, bvar]
				}
			}
			par(cex.main=0.5)
			print(heatmap.2(df, Rowv=T, Colv=F, dendrogram="none", trace="none", col=colorRampPalette(c("blue","white","red"))(150), margins=c(10,15), cexCol=0.6, cexRow=0.4, cellnote=df.sig, notecex=0.8, notecol="red", main=sprintf("Biopsy RF (%s, %s, imp cutoff=%.3g)", st, mlevel, cutoff)))
		}
	}
}
par(cex.main=1.2)


###############################################################################################################
### multi-omics (everything combined) - 01/14/19

## store data for MOFA
mlevel <- "BIOCHEMICAL"; level <- "Genus"
selected <- list()
for (st in c("fecal", "plasma")) {
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% rownames(data))
	selected[[st]] <- unique(mapping.sel$StudyID)
}
mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"))
selected[["oral16S"]] <- unique(subset(mapping.sel, SampleType2=="OralSwab")$StudyID)
selected[["rectal16S"]] <- unique(subset(mapping.sel, SampleType2=="OralSwab")$StudyID)
all_subjects <- Reduce(union, selected)

data.MOFA <- list()
for (st in c("fecal", "plasma")) {
	tmp <- df.metabolon[[st]][[mlevel]]
	colnames(tmp) <- sprintf("[%s] %s", st, colnames(tmp))
	# pad with NAs
	rn <- rownames(tmp); adds <- setdiff(all_subjects, rn); pad <- matrix(NA, nrow=length(adds), ncol=ncol(tmp))
	tmp <- rbind(tmp, pad); rownames(tmp) <- c(rn, adds)
	data.MOFA[[st]] <- t(tmp)[, all_subjects]
}
for (st in c("OralSwab", "RectalSwab")) {
	ps.sel <- subset_samples(ps.relative, SampleType2==st)
	otu.filt <- as.data.frame(otu_table(ps.sel))
  otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
  agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
  lvl <- agg[[level]]
  agg <- agg[,-1]
  rownames(agg) <- lvl
  agg <- agg[which(rowSums(agg > 0.01) >= 3),]
  agg <- agg[rownames(agg),]
  colnames(agg) <- unlist(lapply(colnames(agg), function(x) unlist(strsplit(x, "\\."))[1]))
  rownames(agg) <- sprintf("[%s] %s", st, rownames(agg))
  # pad with NAs
	cn <- colnames(agg); adds <- setdiff(all_subjects, cn); pad <- matrix(NA, ncol=length(adds), nrow=nrow(agg))
	agg <- cbind(agg, pad); colnames(agg) <- c(cn, adds)
  data.MOFA[[st]] <- agg[, all_subjects]
}
tmp <- subset(mapping, SampleType2=="OralSwab"); rownames(tmp) <- gsub("\\.OralSwab", "", rownames(tmp))
covars <- tmp[all_subjects, colnames(mapping.sel)]
# save data for MOFA
save(data.MOFA, covars, file=sprintf("%s/data_for_MOFA.031119.RData", out_dir))


## store data for mixOmics using list of subjects that had everything run (fecal metabolites, plasma metabolites, 16S oral, 16S gut)
mlevel <- "BIOCHEMICAL"; level <- "Genus"
selected <- list()
for (st in c("fecal", "plasma")) {
	data <- df.metabolon[[st]][[mlevel]]
	mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"), StudyID %in% rownames(data))
	selected[[st]] <- unique(mapping.sel$StudyID)
}
mapping.sel <- subset(as(sample_data(ps.relative), "data.frame"))
selected[["oral16S"]] <- unique(subset(mapping.sel, SampleType2=="OralSwab")$StudyID)
selected[["rectal16S"]] <- unique(subset(mapping.sel, SampleType2=="OralSwab")$StudyID)
selected_subjects <- Reduce(intersect, selected)

data <- {}
data.mixomics <- list()
for (st in c("fecal", "plasma")) {
	tmp <- df.metabolon[[st]][[mlevel]]
	tmp <- tmp[selected_subjects,]
	colnames(tmp) <- sprintf("[%s] %s", st, colnames(tmp))
	data <- cbind(data, tmp)
	data.mixomics[[st]] <- tmp
}
for (st in c("OralSwab", "RectalSwab")) {
	ps.sel <- subset_samples(ps.relative, SampleType2==st & StudyID %in% selected_subjects)
	otu.filt <- as.data.frame(otu_table(ps.sel))
  otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
  agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
  lvl <- agg[[level]]
  agg <- agg[,-1]
  rownames(agg) <- lvl
  agg <- agg[which(rowSums(agg > 0.01) >= 3),]
  agg <- agg[rownames(agg),]
  colnames(agg) <- unlist(lapply(colnames(agg), function(x) unlist(strsplit(x, "\\."))[1]))
  agg <- t(agg[, selected_subjects])
  colnames(agg) <- sprintf("[%s] %s", st, colnames(agg))
  data <- cbind(data, agg)
  data.mixomics[[st]] <- agg
}
response <- unlist(sample_data(ps)[match(selected_subjects, sample_data(ps)$StudyID), "Cohort"]); names(response) <- selected_subjects
# save data for mixOmics
save(data.mixomics, response, file=sprintf("%s/data_for_mixOmics.011819.RData", out_dir))

# add BMI, Age, Sex as covariates
data.sel <- as.data.frame(data)
data.sel$BMI <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "BMI"])
data.sel$Age <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Age"])
data.sel$Sex <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Sex"])
agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Sex")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "feature", "value")

set.seed(sum(dim(data.sel)))

# random forests classification (multiclass); MULTIOMICS
## after running for the first time, COMMENT OUT THIS BLOCK ##
num_iter <- 100
ncores <- 20
out <- mclapply(1:num_iter, function (dummy) {
		importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
	}, mc.cores=ncores )
collated.importance <- do.call(cbind, out)
out <- mclapply(1:num_iter, function (dummy) {
		rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
	}, mc.cores=ncores )
collated.cv <- do.call(cbind, out)

write.table(collated.importance, file=sprintf("%s/randomForest_MULTIOMICS.%s.importance.txt", out_dir, "multiclass"), quote=F, sep="\t", row.names=T, col.names=F)
write.table(collated.cv, file=sprintf("%s/randomForest_MULTIOMICS.%s.cv.txt", out_dir, "multiclass"), quote=F, sep="\t", row.names=T, col.names=F)
## END BLOCK TO COMMENT ##

collated.importance <- read.table(sprintf("%s/randomForest_MULTIOMICS.%s.importance.txt", out_dir, "multiclass"), header=F, as.is=T, sep="\t", row.names=1)
collated.cv <- read.table(sprintf("%s/randomForest_MULTIOMICS.%s.cv.txt", out_dir, "multiclass"), header=F, as.is=T, sep="\t", row.names=1)
importance.mean <- rowMeans(collated.importance)
importance.sd <- unlist(apply(collated.importance, 1, sd))
cv.mean <- rowMeans(collated.cv)
cv.sd <- unlist(apply(collated.cv, 1, sd))
inds <- order(importance.mean, decreasing=T)
write.table(melt(importance.mean[inds]), file=sprintf("%s/randomForest_MULTIOMICS.%s.features.txt", out_dir, "multiclass"), quote=F, sep="\t", row.names=T, col.names=F)
inds <- inds[1:as.numeric(names(cv.mean)[which.min(cv.mean)])] # cross-validation suggests using 16 features instead of 12 that meet the 0.001 cutoff
#inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.001)))] # edit as appropriate

# after running for the first time, COMMENT OUT THIS BLOCK ##
# using a sparse model with N predictors
sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
save(sparseRF, file=sprintf("%s/randomForest_MULTIOMICS.%s.model", out_dir, "multiclass"))
load(sprintf("%s/randomForest_MULTIOMICS.%s.model", out_dir, "multiclass"))
# accuracy of final sparseRF model
pred <- predict(sparseRF, type="prob")
pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
pred_df_out <- merge(pred_df, data.sel, by="row.names")
write.table(pred_df_out, file=sprintf("%s/randomForest_MULTIOMICS.%s.predictions.txt", out_dir, "multiclass"), quote=F, sep="\t", row.names=F, col.names=T)
confusion_matrix <- table(pred_df[, c("true", "predicted")])
class_errors <- unlist(lapply(levels(mapping.sel$Cohort), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$Cohort)
accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
mccvalue <- mcc(vec.pred, vec.true)
df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s) (accuracy = %.2f%%, MCC = %.4f)", "multiclass", "MULTIOMICS", accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)
write.table(confusion_matrix, file=sprintf("%s/randomForest_MULTIOMICS.%s.confusion_matrix.txt", out_dir, "multiclass"), quote=F, sep="\t", row.names=T, col.names=T)
## END BLOCK TO COMMENT ##

# plotting - per-group sparse model
df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s", "multiclass", "MULTIOMICS")))
# plotting - per-group variables
df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
df$metabolite_name <- as.character(df$OTU)
df$feature_class <- unlist(lapply(df$metabolite_name, function(x) unlist(strsplit(x, " "))[1]))
p <- ggplot(df, aes(x=OTU, y=importance, label=OTU)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory features", "multiclass", "MULTIOMICS")) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
print(p) 
# shading rectangles of importance values
df.rect <- df
df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s explanatory features", "multiclass", "MULTIOMICS")) + scale_fill_gradient(low="white", high="black") + scale_color_manual(values=cols.sig)
print(p)

# violin plots of relative abundance
agg.melt <- agg.melt.stored
agg.melt$Cohort <- mapping.sel[match(agg.melt$SampleID, mapping.sel$StudyID), "Cohort"]
agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$Cohort))
agg.melt <- subset(agg.melt, feature %in% rownames(df))
agg.melt$feature <- factor(agg.melt$feature, levels=rownames(df))
agg.melt$L099 <- ifelse(agg.melt$SampleID=="L099", "L099", "not L099")
p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~feature, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s)", "MULTIOMICS")) + coord_flip() + scale_color_manual(values=cols.cohort)
print(p)
p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=L099)) + geom_violin() + geom_point(size=1) + facet_wrap(~feature, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s)", "MULTIOMICS")) + coord_flip() + scale_color_manual(values=c("red", "grey"))
print(p)
# separate violin plots for each omics type
for (fc in unique(df$feature_class)) {
	agg.melt2 <- subset(agg.melt, feature %in% rownames(subset(df, feature_class==fc)))
	p <- ggplot(agg.melt2, aes(x=feature, y=value, color=Cohort)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s)", "MULTIOMICS", fc)) + coord_flip() + scale_color_manual(values=cols.cohort)
	print(p)
}
# combined violin plot of all non-BMI features
agg.melt2 <- subset(agg.melt, feature %in% rownames(subset(df, feature_class!="BMI")))
p <- ggplot(agg.melt2, aes(x=feature, y=value, color=Cohort)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (%s, %s)", "MULTIOMICS", "non-BMI")) + coord_flip() + scale_color_manual(values=cols.cohort)
print(p)
# violin plots stratified by prediction/truth
for (met in levels(agg.melt$feature)) {
	tmp <- subset(agg.melt, feature==met)
	p <- ggplot(tmp, aes(x=Cohort, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s)", met, "MULTIOMICS")) + scale_fill_manual(values=cols.cohort)
	print(p)	
	# color scheme - manual
	p <- ggplot(tmp, aes(x=Cohort, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s)", met, "MULTIOMICS")) + scale_fill_manual(values=cols.cohort)
	print(p)
	p <- ggplot(tmp, aes(x=Cohort, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s)", met, "MULTIOMICS")) + scale_color_manual(values=cols.cohort)
	print(p)
}

## random forests classification (one-vs-one setup, X vs ObeseNML baseline); MULTIOMICS
res.mean <- {}; res.sd <- {}; res.sig <- list()
for (mvar_level in setdiff(levels(mapping.sel[, "Cohort"]), c("NMLBMI", "ObeseNML"))) {
	sel <- intersect(rownames(data), unique(subset(sample_data(ps), Cohort %in% c(mvar_level, "ObeseNML"))$StudyID))
	data.sel <- as.data.frame(data[sel,])
	response <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Cohort"]); names(response) <- sel
	# add BMI, Age, Sex as covariates
	data.sel$BMI <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "BMI"])
	data.sel$Age <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Age"])
	data.sel$Sex <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Sex"])
	agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Sex")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "feature", "value")
	
	## after running for the first time, COMMENT OUT THIS BLOCK ##
	num_iter <- 100
	ncores <- 20
	out <- mclapply(1:num_iter, function (dummy) {
			importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
		}, mc.cores=ncores )
	collated.importance <- do.call(cbind, out)
	out <- mclapply(1:num_iter, function (dummy) {
			rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
		}, mc.cores=ncores )
	collated.cv <- do.call(cbind, out)

	write.table(collated.importance, file=sprintf("%s/randomForest_MULTIOMICS_ObeseNML.%s.importance.txt", out_dir, mvar_level), quote=F, sep="\t", row.names=T, col.names=F)
	write.table(collated.cv, file=sprintf("%s/randomForest_MULTIOMICS_ObeseNML.%s.cv.txt", out_dir, mvar_level), quote=F, sep="\t", row.names=T, col.names=F)
	## END BLOCK TO COMMENT ##

	collated.importance <- read.table(sprintf("%s/randomForest_MULTIOMICS_ObeseNML.%s.importance.txt", out_dir, mvar_level), header=F, as.is=T, sep="\t", row.names=1)
	collated.cv <- read.table(sprintf("%s/randomForest_MULTIOMICS_ObeseNML.%s.cv.txt", out_dir, mvar_level), header=F, as.is=T, sep="\t", row.names=1)
	importance.mean <- rowMeans(collated.importance)
	importance.sd <- unlist(apply(collated.importance, 1, sd))
	cv.mean <- rowMeans(collated.cv)
	cv.sd <- unlist(apply(collated.cv, 1, sd))
	res.mean <- rbind(res.mean, importance.mean)
	res.sd <- rbind(res.sd, importance.sd)
	inds <- order(importance.mean, decreasing=T)
	write.table(melt(importance.mean[inds]), file=sprintf("%s/randomForest_MULTIOMICS_ObeseNML.%s.features.txt", out_dir, mvar_level), quote=F, sep="\t", row.names=T, col.names=F)
	inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.001)))] # edit as appropriate
	res.sig[[length(res.sig)+1]] <- names(importance.mean[inds])

	## after running for the first time, COMMENT OUT THIS BLOCK ##
	# using a sparse model with N predictors
	sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, ntree=10000, importance=T)
	save(sparseRF, file=sprintf("%s/randomForest_MULTIOMICS_ObeseNML.%s.model", out_dir, mvar_level))
	load(sprintf("%s/randomForest_MULTIOMICS_ObeseNML.%s.model", out_dir, mvar_level))
	# ROC and AUC of final sparseRF model
	pred <- predict(sparseRF, type="prob")
	pred2 <- prediction(pred[,2], ordered(response))
	perf <- performance(pred2, "tpr", "fpr")
	perf.auc <- performance(pred2, "auc")
	pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
	confusion_matrix <- table(pred_df[, c("true", "predicted")])
	accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
	vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
	mccvalue <- mcc(vec.pred, vec.true)
	plot(perf, main=sprintf("ROC %s %s (sparseRF final model ObeseNML baseline)", mvar_level, "MULTIOMICS")) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g, accuracy=%.2f%%, MCC=%.4g", unlist(perf.auc@y.values), accuracy, mccvalue))			
	## END BLOCK TO COMMENT ##
	
	# plotting - per-group sparse model
	df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
	colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
	print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s", mvar_level, "MULTIOMICS")))
	# plotting - per-group variables
	df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
	df$metabolite_name <- as.character(df$OTU)
	df$feature_class <- unlist(lapply(df$metabolite_name, function(x) unlist(strsplit(x, " "))[1]))
	p <- ggplot(df, aes(x=OTU, y=importance, label=OTU)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory features", mvar_level, "MULTIOMICS")) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
	print(p) 
	# shading rectangles of importance values
	df.rect <- df
	df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
	df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
	p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s explanatory features", mvar_level, "MULTIOMICS")) + scale_fill_gradient(low="white", high="black") + scale_color_manual(values=cols.sig)
	print(p)
	# violin plots of relative abundance
	agg.melt <- agg.melt.stored
	agg.melt$Cohort <- mapping.sel[match(agg.melt$SampleID, mapping.sel$StudyID), "Cohort"]
	agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$Cohort))
	agg.melt <- subset(agg.melt, feature %in% rownames(df))
	agg.melt$feature <- factor(agg.melt$feature, levels=rownames(df))
	agg.melt$L099 <- ifelse(agg.melt$SampleID=="L099", "L099", "not L099")
	p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~feature, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (ObeseNML %s, %s)", "MULTIOMICS", mvar_level)) + coord_flip() + scale_color_manual(values=cols.cohort)
	print(p)
	p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=L099)) + geom_violin() + geom_point(size=1) + facet_wrap(~feature, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (ObeseNML %s, %s)", "MULTIOMICS", mvar_level)) + coord_flip() + scale_color_manual(values=c("red", "grey"))
	print(p)
	# separate violin plots for each omics type
	for (fc in unique(df$feature_class)) {
		agg.melt2 <- subset(agg.melt, feature %in% rownames(subset(df, feature_class==fc)))
		p <- ggplot(agg.melt2, aes(x=feature, y=value, color=Cohort)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (ObeseNML %s, %s, %s)", "MULTIOMICS", fc, mvar_level)) + coord_flip() + scale_color_manual(values=cols.cohort)
		print(p)
	}
	# combined violin plot of all non-BMI features
	agg.melt2 <- subset(agg.melt, feature %in% rownames(subset(df, feature_class!="BMI")))
	p <- ggplot(agg.melt2, aes(x=feature, y=value, color=Cohort)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (ObeseNML %s, %s, %s)", "MULTIOMICS", "non-BMI", mvar_level)) + coord_flip() + scale_color_manual(values=cols.cohort)
	print(p)
}
# plotting - all groups importance values
rownames(res.mean) <- setdiff(levels(mapping.sel[, "Cohort"]), c("NMLBMI", "ObeseNML")); rownames(res.sd) <- rownames(res.mean)
df <- res.mean[, unique(unlist(res.sig))]; df[which(df<0.001, arr.ind=T)] <- 0
df.sig <- matrix("", nrow=nrow(df), ncol=ncol(df)); rownames(df.sig) <- rownames(df); colnames(df.sig) <- colnames(df)
for (i in 1:length(res.sig)) {
	df.sig[i, res.sig[[i]]] <- "*"
}
heatmap.2(df, Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(10,10), cexCol=0.6, cexRow=0.6, cellnote=df.sig, notecex=0.8, notecol="red", main=sprintf("RF importance values (%s, %s)", "ObeseNML", "MULTIOMICS"))


## random forests classification (NASH vs NAFLD); MULTIOMICS
mvar_level <- "NASH"
sel <- intersect(rownames(data), unique(subset(sample_data(ps), Cohort %in% c(mvar_level, "NAFLD"))$StudyID))
data.sel <- as.data.frame(data[sel,])
response <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Cohort"]); names(response) <- sel
# add BMI, Age, Sex as covariates
data.sel$BMI <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "BMI"])
data.sel$Age <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Age"])
data.sel$Sex <- unlist(sample_data(ps)[match(sel, sample_data(ps)$StudyID), "Sex"])
agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Sex")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "feature", "value")

## after running for the first time, COMMENT OUT THIS BLOCK ##
num_iter <- 100
ncores <- 20
out <- mclapply(1:num_iter, function (dummy) {
		importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
	}, mc.cores=ncores )
collated.importance <- do.call(cbind, out)
out <- mclapply(1:num_iter, function (dummy) {
		rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
	}, mc.cores=ncores )
collated.cv <- do.call(cbind, out)

write.table(collated.importance, file=sprintf("%s/randomForest_MULTIOMICS_NAFLD.%s.importance.txt", out_dir, mvar_level), quote=F, sep="\t", row.names=T, col.names=F)
write.table(collated.cv, file=sprintf("%s/randomForest_MULTIOMICS_NAFLD.%s.cv.txt", out_dir, mvar_level), quote=F, sep="\t", row.names=T, col.names=F)
## END BLOCK TO COMMENT ##

collated.importance <- read.table(sprintf("%s/randomForest_MULTIOMICS_NAFLD.%s.importance.txt", out_dir, mvar_level), header=F, as.is=T, sep="\t", row.names=1)
collated.cv <- read.table(sprintf("%s/randomForest_MULTIOMICS_NAFLD.%s.cv.txt", out_dir, mvar_level), header=F, as.is=T, sep="\t", row.names=1)
importance.mean <- rowMeans(collated.importance)
importance.sd <- unlist(apply(collated.importance, 1, sd))
cv.mean <- rowMeans(collated.cv)
cv.sd <- unlist(apply(collated.cv, 1, sd))
res.mean <- rbind(res.mean, importance.mean)
res.sd <- rbind(res.sd, importance.sd)
inds <- order(importance.mean, decreasing=T)
write.table(melt(importance.mean[inds]), file=sprintf("%s/randomForest_MULTIOMICS_NAFLD.%s.features.txt", out_dir, mvar_level), quote=F, sep="\t", row.names=T, col.names=F)
inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.001)))] # edit as appropriate

## after running for the first time, COMMENT OUT THIS BLOCK ##
# using a sparse model with N predictors
sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, ntree=10000, importance=T)
save(sparseRF, file=sprintf("%s/randomForest_MULTIOMICS_NAFLD.%s.model", out_dir, mvar_level))
load(sprintf("%s/randomForest_MULTIOMICS_NAFLD.%s.model", out_dir, mvar_level))
# ROC and AUC of final sparseRF model
pred <- predict(sparseRF, type="prob")
pred2 <- prediction(pred[,2], ordered(response))
perf <- performance(pred2, "tpr", "fpr")
perf.auc <- performance(pred2, "auc")
pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
confusion_matrix <- table(pred_df[, c("true", "predicted")])
accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
mccvalue <- mcc(vec.pred, vec.true)
plot(perf, main=sprintf("ROC %s %s (sparseRF final model NAFLD baseline)", mvar_level, "MULTIOMICS")) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g, accuracy=%.2f%%, MCC=%.4g", unlist(perf.auc@y.values), accuracy, mccvalue))			
## END BLOCK TO COMMENT ##

# plotting - per-group sparse model
df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s", mvar_level, "MULTIOMICS")))
# plotting - per-group variables
df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
df$metabolite_name <- as.character(df$OTU)
df$feature_class <- unlist(lapply(df$metabolite_name, function(x) unlist(strsplit(x, " "))[1]))
p <- ggplot(df, aes(x=OTU, y=importance, label=OTU)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory features", mvar_level, "MULTIOMICS")) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
print(p) 
# shading rectangles of importance values
df.rect <- df
df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s explanatory features", mvar_level, "MULTIOMICS")) + scale_fill_gradient(low="white", high="black") + scale_color_manual(values=cols.sig)
print(p)
# violin plots of relative abundance
agg.melt <- agg.melt.stored
agg.melt$Cohort <- mapping.sel[match(agg.melt$SampleID, mapping.sel$StudyID), "Cohort"]
agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$Cohort))
agg.melt <- subset(agg.melt, feature %in% rownames(df))
agg.melt$feature <- factor(agg.melt$feature, levels=rownames(df))
agg.melt$L099 <- ifelse(agg.melt$SampleID=="L099", "L099", "not L099")
p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~feature, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (NAFLD %s, %s)", "MULTIOMICS", mvar_level)) + coord_flip() + scale_color_manual(values=cols.cohort)
print(p)
p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=L099)) + geom_violin() + geom_point(size=1) + facet_wrap(~feature, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (NAFLD %s, %s)", "MULTIOMICS", mvar_level)) + coord_flip() + scale_color_manual(values=c("red", "grey"))
print(p)
# separate violin plots for each omics type
for (fc in unique(df$feature_class)) {
	agg.melt2 <- subset(agg.melt, feature %in% rownames(subset(df, feature_class==fc)))
	p <- ggplot(agg.melt2, aes(x=feature, y=value, color=Cohort)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (NAFLD %s, %s, %s)", "MULTIOMICS", fc, mvar_level)) + coord_flip() + scale_color_manual(values=cols.cohort)
	print(p)
}
# combined violin plot of all non-BMI features
agg.melt2 <- subset(agg.melt, feature %in% rownames(subset(df, feature_class!="BMI")))
p <- ggplot(agg.melt2, aes(x=feature, y=value, color=Cohort)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF features (NAFLD %s, %s, %s)", "MULTIOMICS", "non-BMI", mvar_level)) + coord_flip() + scale_color_manual(values=cols.cohort)
print(p)




## random forests classification (biopsy data); MULTIOMICS
# get selected data
biopsy <- read.table(sprintf("%s/liver_biopsy.txt", data_dir), header=T, row.names=1, as.is=T, sep="\t")
colnames(biopsy) <- c("NAS_score", "Diagnosis", "Steatosis", "LobularInflammation", "PortalInflammation", "Ballooning", "Fibrosis", "Cirrhosis", "Metabolon")
biopsy$Diagnosis2 <- factor(ifelse(biopsy$NAS_score >= 5, "DefiniteNASH", "BorderlineNASH"))
biopsy$Fibrosis2 <- factor(ifelse(biopsy$Fibrosis >= 3, "Advanced", "Mild"), levels=c("Mild", "Advanced"))
biopsy$Steatosis2 <- factor(pmax(biopsy$Steatosis, 1)) # convert 0 -> 1
biopsy$LobularInflammation2 <- factor(pmax(biopsy$LobularInflammation, 1)) # convert 0 -> 1
biopsy$Ballooning2 <- factor(biopsy$Ballooning)
biopsy_variables <- c("Diagnosis2", "Steatosis2", "Fibrosis2", "LobularInflammation2", "Ballooning2")
sel <- intersect(rownames(biopsy), rownames(data))
biopsy <- biopsy[sel,]
data.sel <- as.data.frame(data[sel,])

# classification
set.seed(nrow(biopsy))
for (biopsy_variable in biopsy_variables) {
	response <- biopsy[sel, biopsy_variable]; names(response) <- sel
	# add BMI, Age, Sex as covariates
	data.sel$BMI <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "BMI"])
	data.sel$Age <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Age"])
	data.sel$Sex <- unlist(sample_data(ps)[match(rownames(data.sel), sample_data(ps)$StudyID), "Sex"])
	agg.melt.stored <- melt(as.matrix(data.sel[, setdiff(colnames(data.sel), "Sex")]), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "feature", "value")

	## after running for the first time, COMMENT OUT THIS BLOCK ##
	num_iter <- 100
	ncores <- 20
	out <- mclapply(1:num_iter, function (dummy) {
			importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
		}, mc.cores=ncores )
	collated.importance <- do.call(cbind, out)
	out <- mclapply(1:num_iter, function (dummy) {
			rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
		}, mc.cores=ncores )
	collated.cv <- do.call(cbind, out)

	write.table(collated.importance, file=sprintf("%s/randomForest_BIOPSY_MULTIOMICS.%s.importance.txt", out_dir, biopsy_variable), quote=F, sep="\t", row.names=T, col.names=F)
	write.table(collated.cv, file=sprintf("%s/randomForest_BIOPSY_MULTIOMICS.%s.cv.txt", out_dir, biopsy_variable), quote=F, sep="\t", row.names=T, col.names=F)
	## END BLOCK TO COMMENT ##

	collated.importance <- read.table(sprintf("%s/randomForest_BIOPSY_MULTIOMICS.%s.importance.txt", out_dir, biopsy_variable), header=F, as.is=T, sep="\t", row.names=1)
	collated.cv <- read.table(sprintf("%s/randomForest_BIOPSY_MULTIOMICS.%s.cv.txt", out_dir, biopsy_variable), header=F, as.is=T, sep="\t", row.names=1)
	importance.mean <- rowMeans(collated.importance)
	importance.sd <- unlist(apply(collated.importance, 1, sd))
	cv.mean <- rowMeans(collated.cv)
	cv.sd <- unlist(apply(collated.cv, 1, sd))
	inds <- order(importance.mean, decreasing=T)
	write.table(melt(importance.mean[inds]), file=sprintf("%s/randomForest_BIOPSY_MULTIOMICS.%s.features.txt", out_dir, biopsy_variable), quote=F, sep="\t", row.names=T, col.names=F)
	inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.001)))] # edit as appropriate
	

	## after running for the first time, COMMENT OUT THIS BLOCK ##
	# using a sparse model with N predictors
	sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, ntree=10000, importance=T)
	save(sparseRF, file=sprintf("%s/randomForest_BIOPSY_MULTIOMICS.%s.model", out_dir, biopsy_variable))
	load(sprintf("%s/randomForest_BIOPSY_MULTIOMICS.%s.model", out_dir, biopsy_variable))
	# accuracy of final sparseRF model
	pred <- predict(sparseRF, type="prob")
	pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
	pred_df_out <- merge(pred_df, data.sel, by="row.names")
	write.table(pred_df_out, file=sprintf("%s/randomForest_BIOPSY_MULTIOMICS.%s.predictions.txt", out_dir, biopsy_variable), quote=F, sep="\t", row.names=F, col.names=T)
	confusion_matrix <- table(pred_df[, c("true", "predicted")])
	class_errors <- unlist(lapply(levels(biopsy[,biopsy_variable]), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(biopsy[,biopsy_variable])
	accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
	vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
	mccvalue <- mcc(vec.pred, vec.true)
	df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s) (accuracy = %.2f%%, MCC = %.4f)", biopsy_variable, "MULTIOMICS", accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
	print(p)
	write.table(confusion_matrix, file=sprintf("%s/randomForest_BIOPSY_MULTIOMICS.%s.confusion_matrix.txt", out_dir, biopsy_variable), quote=F, sep="\t", row.names=T, col.names=T)
	# ROC analysis
	if (nlevels(response)==2) {
		pred <- predict(sparseRF, type="prob")
		pred2 <- prediction(pred[,2], ordered(response))
		perf <- performance(pred2, "tpr", "fpr")
		perf.auc <- performance(pred2, "auc")
		print(plot(perf, main=sprintf("ROC %s %s (sparseRF final model)", biopsy_variable, "MULTIOMICS")) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", unlist(perf.auc@y.values))))
	}
	## END BLOCK TO COMMENT ##

	# plotting - per-group sparse model
	df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
	colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
	print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s", biopsy_variable, "MULTIOMICS")))
	# plotting - per-group variables
	df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
	df$metabolite_name <- as.character(df$OTU)
	df$feature_class <- unlist(lapply(df$metabolite_name, function(x) unlist(strsplit(x, " "))[1]))
	p <- ggplot(df, aes(x=OTU, y=importance, label=OTU)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory features", biopsy_variable, "MULTIOMICS")) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
	print(p) 
	# shading rectangles of importance values
	df.rect <- df
	df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
	df.rect$importance <- pmin(df.rect$importance, sort(df.rect$importance, decreasing=T)[2]) # censor importance value to 2nd highest level (drops BMI from coloring)
	p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("%s - %s explanatory features", biopsy_variable, "MULTIOMICS")) + scale_fill_gradient(low="white", high="black") + scale_color_manual(values=cols.sig)
	print(p)
	# violin plots of relative abundance
	agg.melt <- agg.melt.stored
	agg.melt$Truth <- biopsy[agg.melt$SampleID, biopsy_variable]
	agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt$Truth))
	agg.melt <- subset(agg.melt, feature %in% rownames(df))
	agg.melt$feature <- factor(agg.melt$feature, levels=rownames(df))
	p <- ggplot(agg.melt, aes(x=Truth, y=value, color=Truth)) + geom_violin() + geom_point() + facet_wrap(~feature, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", biopsy_variable, "MULTIOMICS")) + coord_flip()
	print(p)
	# violin plots stratified by prediction/truth
	for (met in levels(agg.melt$feature)) {
		tmp <- subset(agg.melt, feature==met)
		# color scheme - Set1
		p <- ggplot(tmp, aes(x=Truth, y=value, fill=Prediction)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s)", met, "MULTIOMICS")) + scale_fill_brewer(palette="Set1")
		print(p)
		p <- ggplot(tmp, aes(x=Truth, y=value, fill=Prediction, color=Prediction)) + geom_point(position=position_jitter(width=0.1)) + theme_classic() + ggtitle(sprintf("Rel. abund. of %s (%s)", met, "MULTIOMICS")) + scale_color_brewer(palette="Set1")
		print(p)
	}
}


# heatmap of all biopsy classification models
for (cutoff in c(0.001, 0.002, 0.005)) {
	res <- {}
	for (biopsy_variable in biopsy_variables) {
		collated.importance <- read.table(sprintf("%s/randomForest_BIOPSY_MULTIOMICS.%s.importance.txt", out_dir, biopsy_variable), header=F, as.is=T, sep="\t", row.names=1)
		collated.cv <- read.table(sprintf("%s/randomForest_BIOPSY_MULTIOMICS.%s.cv.txt", out_dir, biopsy_variable), header=F, as.is=T, sep="\t", row.names=1)
		importance.mean <- rowMeans(collated.importance)
		importance.sd <- unlist(apply(collated.importance, 1, sd))
		cv.mean <- rowMeans(collated.cv)
		cv.sd <- unlist(apply(collated.cv, 1, sd))
		inds <- order(importance.mean, decreasing=T)
		inds <- inds[1:min(30, length(which(importance.mean[inds] > cutoff)))]
		tmp <- data.frame(feature=names(importance.mean[inds]), importance=importance.mean[inds], mvar=biopsy_variable)
		res <- rbind(res, tmp)
	}
	df <- dcast(res, feature ~ mvar, value.var="importance")
	rownames(df) <- df$feature; df <- as.matrix(df[,-1])
#			df <- log2(df)
	df[which(is.na(df), arr.ind=T)] <- 0 # NA -> 0
	par(cex.main=0.5)
	print(heatmap.2(df, Rowv=T, Colv=F, dendrogram="none", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(10,15), cexCol=0.6, cexRow=0.4, main=sprintf("Biopsy RF (%s, imp cutoff=%.3g)", "MULTIOMICS", cutoff)))
	print(heatmap.2(df, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(10,15), cexCol=0.6, cexRow=0.4, main=sprintf("Biopsy RF (%s, imp cutoff=%.3g)", "MULTIOMICS", cutoff)))
}
par(cex.main=1.2)



## mixOmics
# load data (in new R 3.5.2 environment because all other analysis is in R 3.4.3 but mixOmics is not available in older R)
library(mixOmics)
cols.cohort <- c("#888888", "orange", "red", "#a400a4"); names(cols.cohort) <- c("NMLBMI", "ObeseNML", "NAFLD", "NASH")
load(sprintf("%s/data_for_mixOmics.011819.RData", out_dir))
data <- data.mixomics
set.seed(length(data))

# set up design matrix
design <- matrix(0.1, ncol = length(data), nrow = length(data), dimnames = list(names(data), names(data)))
diag(design) <- 0

# tune the number of components
sgccda.res <- block.splsda(X = data, Y = response, ncomp = 5, design = design)
perf.diablo <- perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
plot(perf.diablo)
ncomp <- perf.diablo$choice.ncomp$WeightedPredict["Overall.BER", ]

# tune keepX (number of parameters to keep from each data type)
test.keepX <- list(fecal = c(5:9, seq(10, 18, 2), seq(20,30,5)), plasma = c(5:9, seq(10, 18, 2), seq(20,30,5)), OralSwab = c(5:9, seq(10, 18, 2), seq(20,30,5)), RectalSwab = c(5:9, seq(10, 18, 2), seq(20,30,5)))
tune <- tune.block.splsda(X = data, Y = response, ncomp = ncomp, test.keepX = test.keepX, design = design, validation = 'Mfold', folds = 10, nrepeat = 1, cpus = 16, dist = "mahalanobis.dist")

# final model
sgccda.res <- block.splsda(X = data, Y = response, ncomp = ncomp, keepX = tune$choice.keepX, design = design)

# diagnostic plots
plotDiablo(sgccda.res, ncomp = 1)
plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO', ellipse=TRUE)
plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE)
circosPlot(sgccda.res, cutoff = 0.7, line = TRUE, size.labels = 1.5, comp=1, color.Y=cols.cohort)
circosPlot(sgccda.res, cutoff = 0.8, line = TRUE, size.labels = 1.5, comp=1)
circosPlot(sgccda.res, cutoff = 0.7, line = TRUE, size.labels = 1.5, comp=2)
circosPlot(sgccda.res, cutoff = 0.8, line = TRUE, size.labels = 1.5, comp=2)

plotLoadings(sgccda.res, comp = 1, contrib = 'max', method = 'median')
plotLoadings(sgccda.res, comp = 2, contrib = 'max', method = 'median')
cimDiablo(sgccda.res, comp=1, margins=c(4,20), size.legend=0.8, transpose=T, legend.position="topright", color.Y=cols.cohort)
cimDiablo(sgccda.res, comp=2, margins=c(4,20), size.legend=0.8, transpose=T, legend.position="topright", color.Y=cols.cohort)

for (b in rownames(design)) {
	auc.splsda <- auroc(sgccda.res, roc.block = b, roc.comp = 1)
	auc.splsda <- auroc(sgccda.res, roc.block = b, roc.comp = 2)
}


## MOFA (https://github.com/bioFAM/MOFA)
library(MOFA)
library(MultiAssayExperiment)

# load data
cols.cohort <- c("#888888", "orange", "red", "#a400a4"); names(cols.cohort) <- c("NMLBMI", "ObeseNML", "NAFLD", "NASH")
load(sprintf("%s/data_for_MOFA.031119.RData", out_dir))

# set up MOFA object
mae <- MultiAssayExperiment(experiments = data.MOFA, colData = covars)
MOFAobject <- createMOFAobject(mae)
MOFAobject
DataOptions <- getDefaultDataOptions()
ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 25
TrainOptions <- getDefaultTrainOptions()
TrainOptions$maxiter <- 15000
TrainOptions$DropFactorThreshold <- 0.02
TrainOptions$tolerance <- 0.1
TrainOptions$seed <- length(data.MOFA)
MOFAobject <- prepareMOFA(MOFAobject, DataOptions=DataOptions, ModelOptions=ModelOptions, TrainOptions=TrainOptions)
# save featureNames to fix truncation during HDF5 save
featurenames <- MOFA::featureNames(MOFAobject)

# train MOFA model
#MOFAobject <- runMOFA(MOFAobject, outfile=sprintf("%s/MOFAobject.hdf5", out_dir))
MOFAobject <- loadModel(sprintf("%s/MOFAobject.hdf5", out_dir), MOFAobject)
MOFA::featureNames(MOFAobject) <- featurenames[names(MOFA::featureNames(MOFAobject))]

# secondary anaylsis
r2 <- calculateVarianceExplained(MOFAobject)
plotVarianceExplained(MOFAobject)
for (viewname in viewNames(MOFAobject)) {
	plotWeightsHeatmap(MOFAobject, view=viewname, factors="all", show_colnames=T, main=sprintf("%s", viewname), fontsize_row=8, fontsize_col=4)
}
for (viewname in viewNames(MOFAobject)) {
	for (lf in factorNames(MOFAobject)) {
		print(plotWeights(MOFAobject, view=viewname, factor=lf, nfeatures=10) + ggtitle(sprintf("%s (%s)", viewname, lf)) + theme(text=element_text(size=6)))
	}
}
for (viewname in viewNames(MOFAobject)) {
	for (lf in factorNames(MOFAobject)) {
		print(plotTopWeights(MOFAobject, view=viewname, factor=lf))
	}
}
for (viewname in viewNames(MOFAobject)) {
	for (lf in factorNames(MOFAobject)) {
		plotDataHeatmap(MOFAobject, view=viewname, factor=lf, features=20, show_rownames=T, main=sprintf("%s (%s)", viewname, lf), fontsize_row=4, fontsize_col=6)
	}
}
plotFactorScatter(MOFAobject, factors=1:2, color_by="Cohort")
plotFactorScatters(MOFAobject, factors=1:3, color_by="Cohort")


plotFactorBeeswarm(MOFAobject, factor=1, color_by="Cohort")
clusters <- clusterSamples(MOFAobject, k=2, factors=1)



##########################################################################################################################
### Shotgun analysis
library(phyloseq)
library(DESeq2)
library(randomForest)
library(vegan)
library(ggplot2)
library(reshape2)
library(ROCR)

data <- read.table(sprintf("%s/merged_abundance_table.kraken_raw.L7.txt", data_dir), header=T, as.is=T, sep="\t", comment.char="", row.names=1, quote="")
data <- data[, setdiff(1:ncol(data), grep(".Bacterial.kraken.1", colnames(data)))]
colnames(data) <- gsub("PCMP_", "", gsub(".Bacterial.kraken", "", colnames(data)))
rownames(data) <- gsub("\\|", ";", rownames(data))

mapping <- read.table(sprintf("%s/NAFLDomics_Mapping.combined.txt", data_dir), header=T, as.is=T, sep="\t", comment.char="", row.names=1)
metadata <- read.delim(sprintf("%s/metadata.052218.txt", data_dir), header=T, sep="\t", comment.char=""); colnames(metadata)[1] <- "StudyID"
# impute missing BMI values
metadata$BMI....kg.m.2. <- as.numeric(as.character(metadata$BMI....kg.m.2.))
metadata[which(is.na(metadata$BMI....kg.m.2.)), "BMI....kg.m.2."] <- mean(metadata$BMI....kg.m.2., na.rm=T)
mapping$Cohort <- as.character(metadata$Cohort[match(mapping$StudyID, metadata$StudyID)])
mapping$BMI <- metadata[match(mapping$StudyID, metadata$StudyID), "BMI....kg.m.2."]
mapping$Age <- metadata[match(mapping$StudyID, metadata$StudyID), "Age.at.Date.of.Visit"]
mapping$Sex <- metadata[match(mapping$StudyID, metadata$StudyID), "Sex"]

agg <- aggregate(cbind(Cohort, BMI, Age) ~ StudyID, mapping, unique); agg2 <- aggregate(Sex ~ StudyID, mapping, unique)
mapping <- merge(agg, agg2, by="StudyID"); rownames(mapping) <- mapping$StudyID
# set up metadata
metadata_variables <- read.table(sprintf("%s/metadata_types.txt", data_dir), header=T, as.is=T, sep="\t", row.names=1)
sel <- intersect(rownames(metadata_variables), colnames(mapping))
metadata_variables <- metadata_variables[sel,, drop=F]
mapping.sel <- subset(mapping, StudyID %in% colnames(data))
## fix column types
for (mvar in rownames(metadata_variables)) {
	if (metadata_variables[mvar, "type"] == "factor") {
		mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
		if (metadata_variables[mvar, "baseline"] != "") {
			mapping.sel[,mvar] <- relevel(mapping.sel[,mvar], metadata_variables[mvar, "baseline"])
		}
	} else if (metadata_variables[mvar, "type"] == "ordered") {
		if (metadata_variables[mvar, "baseline"] == "") {
			lvls <- unique(mapping.sel[,mvar])
		} else {
			lvls <- unlist(strsplit(metadata_variables[mvar, "baseline"], ","))
		}
		mapping.sel[,mvar] <- ordered(mapping.sel[,mvar], levels=lvls)
	} else if (metadata_variables[mvar, "type"] == "numeric") {
		mapping.sel[,mvar] <- as.numeric(as.character(mapping.sel[,mvar]))
	}
}
mapping.sel$Cohort <- droplevels(mapping.sel$Cohort)

ranklist <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxlist <- lapply(rownames(data), function(x) parse_taxonomy_qiime(x))
taxa <- matrix(NA, nrow=nrow(data), ncol=length(ranklist)); colnames(taxa) <- ranklist
for (i in 1:length(taxlist)) {
	taxa[i, names(taxlist[[i]])] <- taxlist[[i]]
}
tt <- tax_table(taxa); rownames(tt) <- rownames(data)
ps <- phyloseq(otu_table(data, taxa_are_rows=T), tt, sample_data(mapping.sel))
data <- as.data.frame(otu_table(ps))
ps.relative <- transform_sample_counts(ps, function(x) x / sum(x) )
ps.rarefied <- rarefy_even_depth(ps, sample.size=40224, rngseed=nsamples(ps))

## write files for FishTaco
otu.filt <- as.data.frame(otu_table(ps.relative))
otu.filt <- otu.filt[which(rowSums(otu.filt > 0.01) >= ceiling(ncol(otu.filt)/10)),]
write.table(otu.filt, file=sprintf("%s/FishTaco/taxa_abundance.L7.txt", out_dir), quote=F, sep="\t", row.names=T, col.names=T)
labels <- mapping.sel[, c("StudyID", "Cohort"), drop=F]; colnames(labels) <- c("Sample", "Label")
write.table(labels, file=sprintf("%s/FishTaco/sample_labels.txt", out_dir), quote=F, sep="\t", row.names=F, col.names=T)

## PCoA + PERMANOVA
set.seed(112717)
ps.sel <- ps.relative
df <- as(sample_data(ps.sel), "data.frame")
for (distance_metric in distance_metrics) {
  ordi <- ordinate(ps.sel, method = "PCoA", distance = distance_metric)
  p <- plot_ordination(ps.sel, ordi, "samples", color = "Cohort") + geom_point(size = 2) + scale_fill_brewer(type = "qual", palette = "Set1") + scale_colour_brewer(type = "qual", palette = "Set1") + theme_classic() + ggtitle(sprintf("%s - %s", "Cohort", distance_metric))
  print(p)
  p <- plot_ordination(ps.sel, ordi, "samples", color = "BMI") + geom_point(size = 2) + theme_classic() + ggtitle(sprintf("%s - %s", "BMI", distance_metric))
  print(p)
}
sink(sprintf("%s/phyloseq_output.%s.%s.txt", out_dir, "NAFLDomics_Shotgun", format(Sys.Date(), "%m%d%y")))
dm <- list()
for (distance_metric in distance_metrics) {
  dm[[length(dm)+1]] <- as.matrix(phyloseq::distance(ps.sel, method=distance_metric))
}
names(dm) <- distance_metrics
print("PERMANOVA manually selected variables")
for (distance_metric in distance_metrics) {
  print(distance_metric)
  form <- as.formula("as.dist(dm[[distance_metric]]) ~ Cohort + BMI + Age + Sex")
  res <- adonis(form , data=df, permutations=999)
  print(res)
}
sink()

## taxa barplots
ordi <- ordinate(ps.sel, method = "PCoA", distance = "bray")
ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))
otu.filt <- normalizeByCols(as.data.frame(otu_table(ps.sel)))
otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level="Genus")
agg <- aggregate(. ~ Genus, otu.filt, sum)
genera <- agg$Genus
agg <- agg[,-1]
agg <- normalizeByCols(agg)
inds_to_grey <- which(rowMeans(agg)<0.01)
genera[inds_to_grey] <- "Other"
agg$Genus <- genera
df <- melt(agg, variable.name="SampleID")
df2 <- aggregate(as.formula("value ~ Genus + SampleID"), df, sum)
df2$SampleID <- as.character(df2$SampleID)
df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
df.SampleIDstr <- unique(df2[,c("SampleID", "SampleIDfactor")])
df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, "Cohort"])
df2$Cohort <- mapping.sel[df2$SampleID, "Cohort"]
p <- ggplot(df2, aes_string(x="SampleIDfactor", y="value", fill="Genus", order="Genus")) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring) + ggtitle(sprintf("Taxa barplots overall")) + guides(col = guide_legend(ncol = 3)) + ylim(c(-.1, 1.01)) + annotate("rect", xmin = as.numeric(df.SampleIDstr$SampleIDfactor)-0.5, xmax = as.numeric(df.SampleIDstr$SampleIDfactor)+0.5, ymin = -0.04, ymax = -0.02, fill=cols.cohort[df.SampleIDstr$Group])
print(p)
p <- ggplot(df2, aes_string(x="SampleID", y="value", fill="Genus", order="Genus")) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring) + ggtitle(sprintf("Taxa barplots overall")) + guides(col = guide_legend(ncol = 3)) + ylim(c(-.1, 1.01))
print(p)
p <- ggplot(df2, aes_string(x="SampleIDfactor", y="value", fill="Genus", order="Genus")) + geom_bar(stat="identity", position="stack") + facet_wrap(~Cohort, scales="free_x") + theme_classic() + theme(legend.position="right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + scale_fill_manual(values=coloring) + ggtitle(sprintf("Taxa barplots overall")) + guides(col = guide_legend(ncol = 3)) + ylim(c(-.1, 1.01))
print(p)

## alpha diversity
adiv <- estimate_richness(ps.rarefied)
adiv <- adiv[levels(df2$SampleIDfactor),]
adiv$SampleID <- rownames(adiv)
adiv$SampleID <- factor(adiv$SampleID, levels=levels(df2$SampleIDfactor))
adiv$Cohort <- mapping.sel[as.character(adiv$SampleID), "Cohort"]
adiv$BMI <- mapping.sel[as.character(adiv$SampleID), "BMI"]
adiv$Age <- mapping.sel[as.character(adiv$SampleID), "Age"]
adiv$Sex <- mapping.sel[as.character(adiv$SampleID), "Sex"]
plotlist <- list()
for (am in alpha_metrics){
  p <- ggplot(adiv, aes_string(x="SampleID", y=am)) + geom_bar(stat="identity", position="identity") + theme_classic() + ggtitle(sprintf("Alpha Diversity (%s)", am)) + theme( axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8))
  print(p)
  test <- wilcox.test(as.formula(sprintf("%s ~ Cohort", am)), adiv)
  sink(sprintf("%s/phyloseq_output.%s.%s.txt", out_dir, "NAFLDomics_Shotgun", format(Sys.Date(), "%m%d%y")), append=T)
  mod <- lm(as.formula(sprintf("%s ~ Cohort + BMI + Age + Sex", am)), adiv)
  print(summary(mod))
  sink()
  p <- ggplot(adiv, aes_string(x="Cohort", y=am)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Alpha Diversity by Cohort (%s, Wilcox p=%.4g lm(Cohort+BMI+Age+Sex) p=%.4g)", am, test$p.value, summary(mod)$coefficients["Cohort.L", "Pr(>|t|)"])) + theme( axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8))
  plotlist[[length(plotlist)+1]] <- p
}
multiplot(plotlist=plotlist, cols=2, rows=2)

## Bacteroidetes/Firmicutes	ratio
level <- "Phylum"
ps.sel <- ps.relative
otu.filt <- as.data.frame(otu_table(ps.sel))
otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.sel), level=level)
agg <- aggregate(as.formula(sprintf(". ~ %s", level)), otu.filt, sum)
lvl <- agg[[level]]
agg <- agg[,-1]
rownames(agg) <- lvl
agg <- agg[which(rowSums(agg > 0.01) >= ceiling(ncol(agg)/10)),]
agg <- normalizeByCols(agg)
df <- melt(agg["Bacteroidetes",] / agg["Firmicutes",]); colnames(df) <- c("SampleID", "BFRatio"); df$SampleID <- as.character(df$SampleID); df$BFLogRatio <- log10(df$BFRatio)
df$Cohort <- mapping[df$SampleID, "Cohort"]; df$Cohort <- factor(as.character(df$Cohort), levels=c("NMLBMI", "NASH"))
df$BMI <- mapping.sel[as.character(df$SampleID), "BMI"]
df$Age <- mapping.sel[as.character(df$SampleID), "Age"]
df$Sex <- mapping.sel[as.character(df$SampleID), "Sex"]
test <- wilcox.test(BFLogRatio ~ Cohort, df)
sink(sprintf("%s/phyloseq_output.%s.%s.txt", out_dir, "NAFLDomics_Shotgun", format(Sys.Date(), "%m%d%y")), append=T)
mod <- lm(as.formula(sprintf("%s ~ Cohort + BMI + Age + Sex", "BFLogRatio")), df)
summary(mod)
sink()
p <- ggplot(df, aes(x=Cohort, y=BFLogRatio)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("BFLogRatio by Cohort (Wilcoxon p=%.4g)", test$p.value))
print(p)


## DESeq2
ps.deseq2 <- ps; sample_data(ps.deseq2)$Cohort <- factor(as.character(sample_data(ps.deseq2)$Cohort), levels=c("NMLBMI", "NASH"))
dds <- phyloseq_to_deseq2(ps.deseq2, ~ Cohort + BMI + Age + Sex)
dds <- DESeq(dds)
res <- results(dds, name="Cohort_NASH_vs_NMLBMI")
resSig <- as.data.frame(subset(res, padj<0.1)); resSig <- resSig[order(resSig$log2FoldChange),]; resSig$taxa <- rownames(resSig); resSig$taxa <- factor(resSig$taxa, levels=resSig$taxa); resSig$hdir <- 1-ceiling(sign(resSig$log2FoldChange)/2)
p <- ggplot(resSig, aes(x=taxa, y=log2FoldChange)) + geom_bar(stat="identity", fill="#aaaaaa") + geom_text(aes(label=taxa), y=0, size=2, hjust=0.5) + coord_flip() + theme_classic() + ggtitle(sprintf("DESeq2 hits")) + theme(axis.text.y=element_blank())
print(p)
write.table(res, file=sprintf("%s/Shotgun_DESeq2.Species.txt", out_dir), quote=F, sep="\t", row.names=T, col.names=T)

## randomForest classification of Cohort (one-vs-one setup, X vs NMLBMI baseline)
level <- "Species"
otu.filt <- as.data.frame(otu_table(ps.relative))
otu.filt[[level]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.relative), level=level)
otu.filt <- otu.filt[which(rowSums(otu.filt > 0.01) >= ceiling(ncol(otu.filt)/10)),]
agg <- otu.filt[rownames(otu.filt),]; agg <- agg[, setdiff(colnames(agg), level)]
agg.melt.stored <- melt(as.matrix(agg), as.is=T); colnames(agg.melt.stored) <- c("taxa", "SampleID", "value")
res.mean <- {}; res.sd <- {}
set.seed(112817)
for (mvar_level in setdiff(levels(mapping.sel[, "Cohort"]), "NMLBMI")) {
	data.sel <- as.data.frame(t(agg))
	data.sel$BMI <- as(sample_data(ps.relative), "data.frame")[rownames(data.sel), "BMI"]
	data.sel$Age <- as(sample_data(ps.relative), "data.frame")[rownames(data.sel), "Age"]
	data.sel$Sex <- as(sample_data(ps.relative), "data.frame")[rownames(data.sel), "Sex"]
	response <- droplevels(mapping.sel[rownames(data.sel), "Cohort"]); names(response) <- rownames(data.sel)

	## after running for the first time, COMMENT OUT THIS BLOCK ##
	num_iter <- 100
	ncores <- 20
	out <- mclapply(1:num_iter, function (dummy) {
			importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
		}, mc.cores=ncores )
	collated.importance <- do.call(cbind, out)
	out <- mclapply(1:num_iter, function (dummy) {
			rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
		}, mc.cores=ncores )
	collated.cv <- do.call(cbind, out)

	write.table(collated.importance, file=sprintf("%s/randomForest_Shotgun.%s.%s.importance.txt", out_dir, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=F)
	write.table(collated.cv, file=sprintf("%s/randomForest_Shotgun.%s.%s.cv.txt", out_dir, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=F)
	## END BLOCK TO COMMENT ##
	
	collated.importance <- read.table(sprintf("%s/randomForest_Shotgun.%s.%s.importance.txt", out_dir, mvar_level, level), header=F, as.is=T, sep="\t", row.names=1)
	collated.cv <- read.table(sprintf("%s/randomForest_Shotgun.%s.%s.cv.txt", out_dir, mvar_level, level), header=F, as.is=T, sep="\t", row.names=1)
	importance.mean <- rowMeans(collated.importance)
	importance.sd <- unlist(apply(collated.importance, 1, sd))
	cv.mean <- rowMeans(collated.cv)
	cv.sd <- unlist(apply(collated.cv, 1, sd))
	res.mean <- rbind(res.mean, importance.mean)
	res.sd <- rbind(res.sd, importance.sd)
	inds <- order(importance.mean, decreasing=T)
	inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.001)))] # edit as appropriate
	
	## after running for the first time, COMMENT OUT THIS BLOCK ##
	# using a sparse model with N predictors
	sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
	save(sparseRF, file=sprintf("%s/randomForest_Shotgun.%s.%s.model", out_dir, mvar_level, level))
	load(sprintf("%s/randomForest_Shotgun.%s.%s.model", out_dir, mvar_level, level))
	# ROC and AUC of final sparseRF model
	pred <- predict(sparseRF, type="prob")
	pred2 <- prediction(pred[,2], ordered(response))
	perf <- performance(pred2, "tpr", "fpr")
	perf.auc <- performance(pred2, "auc")
	plot(perf, main=sprintf("ROC %s %s (sparseRF final model)", mvar_level, level)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", unlist(perf.auc@y.values)))	
	# 5x cross-validation ROC and AUC
	nperm <- 100
	out <- mclapply(1:nperm, function (dummy) {
		tmp <- data.frame(SampleID=rownames(data.sel), Cohort=response)
		tmp$SampleID <- as.character(tmp$SampleID)
		train.sel <- unlist(lapply(levels(tmp$Cohort), function(lvl) {
			tmp2 <- subset(tmp, Cohort==lvl)
			sample(tmp2$SampleID, size=round(nrow(tmp2)*0.8))
		}))
		test.sel <- setdiff(rownames(data.sel), train.sel)
		rf <- randomForest(x=data.sel[train.sel, names(importance.mean[inds])], y=response[train.sel], ntree=10000, importance=T)
		pred <- predict(rf, data.sel[test.sel, names(importance.mean[inds])], type="prob")
		pred[,2]
	}, mc.cores=ncores )
	saveRDS(out, file=sprintf("%s/randomForest_Shotgun.crossvalidation.%s.%s.rds", out_dir, mvar_level, level))
	out <- readRDS(sprintf("%s/randomForest_Shotgun.crossvalidation.%s.%s.rds", out_dir, mvar_level, level))
	pred <- prediction(out, lapply(out, function(x) ordered(response[names(x)])) )
	perf <- performance(pred, "tpr", "fpr")
	perf.auc <- performance(pred, "auc"); auc.mean <- mean(unlist(perf.auc@y.values))
	plot(perf,col="grey82",lty=3, main=sprintf("ROC %s %s (5x cross-validation)", mvar_level, level))
	plot(perf, lwd=3, avg="vertical", spread.estimate="boxplot", add=TRUE) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", auc.mean))			
	## END BLOCK TO COMMENT ##
	
	load(sprintf("%s/randomForest_Shotgun.%s.%s.model", out_dir, mvar_level, level))
	# plotting - per-group sparse model
	df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
	colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
	print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s", mvar_level, level)))
		
	# plotting - per-group variables
	df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
	res <- read.table(sprintf("%s/Shotgun_DESeq2.Species.txt", out_dir), header=T, sep="\t", as.is=T, quote="", comment.char="")
	df$log2FoldChange <- res[as.character(df$OTU), "log2FoldChange"]; df$padj <- res[as.character(df$OTU), "padj"]
	df$OTU_string <- sprintf("%s (log2FC=%.4g, padj=%.4g)", rownames(df), df$log2FoldChange, df$padj)
	df$sig <- ifelse(df$padj < 0.1, "sig", "ns"); df$sig[which(is.na(df$sig))] <- "ns"
	p <- ggplot(df, aes(x=OTU, y=importance, label=OTU)) + geom_bar(position=position_dodge(), stat="identity", fill="#999999", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string, color=sig), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory", mvar_level, level)) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
	print(p)
	
	# violin plots of relative abundance
	agg.melt <- agg.melt.stored
	agg.melt$Cohort <- mapping.sel[match(agg.melt$SampleID, mapping.sel$StudyID), "Cohort"]
	agg.melt <- subset(agg.melt, taxa %in% rownames(df) & Cohort %in% c(mvar_level, "NMLBMI"))
	agg.melt$taxa <- factor(agg.melt$taxa, levels=rownames(df))
	p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~taxa, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF taxa (%s)", mvar_level)) + coord_flip()
	print(p)			
	
	
}
	

dev.off()


#########################################################################################################################
### Shotgun analysis (HUMAnN2)

data <- read.table(sprintf("%s/merged_pathabundance.txt", data_dir), header=T, as.is=T, sep="\t", comment.char="", row.names=1, quote="")
colnames(data) <- gsub("PCMP_", "", gsub(".cat.filtered2_Abundance", "", colnames(data)))
rownames(data) <- gsub("\\|", ";", rownames(data))

mapping <- read.table(sprintf("%s/NAFLDomics_Mapping.combined.txt", data_dir), header=T, as.is=T, sep="\t", comment.char="", row.names=1)
metadata <- read.delim(sprintf("%s/metadata.052218.txt", data_dir), header=T, sep="\t", comment.char=""); colnames(metadata)[1] <- "StudyID"
# impute missing BMI values
metadata$BMI....kg.m.2. <- as.numeric(as.character(metadata$BMI....kg.m.2.))
metadata[which(is.na(metadata$BMI....kg.m.2.)), "BMI....kg.m.2."] <- mean(metadata$BMI....kg.m.2., na.rm=T)
mapping$Cohort <- as.character(metadata$Cohort[match(mapping$StudyID, metadata$StudyID)])
mapping$BMI <- metadata[match(mapping$StudyID, metadata$StudyID), "BMI....kg.m.2."]
mapping$Age <- metadata[match(mapping$StudyID, metadata$StudyID), "Age.at.Date.of.Visit"]
mapping$Sex <- metadata[match(mapping$StudyID, metadata$StudyID), "Sex"]

agg <- aggregate(cbind(Cohort, BMI, Age) ~ StudyID, mapping, unique); agg2 <- aggregate(Sex ~ StudyID, mapping, unique)
mapping <- merge(agg, agg2, by="StudyID"); rownames(mapping) <- mapping$StudyID
# set up metadata
metadata_variables <- read.table(sprintf("%s/metadata_types.txt", data_dir), header=T, as.is=T, sep="\t", row.names=1)
sel <- intersect(rownames(metadata_variables), colnames(mapping))
metadata_variables <- metadata_variables[sel,, drop=F]
mapping.sel <- subset(mapping, StudyID %in% colnames(data))
## fix column types
for (mvar in rownames(metadata_variables)) {
	if (metadata_variables[mvar, "type"] == "factor") {
		mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
		if (metadata_variables[mvar, "baseline"] != "") {
			mapping.sel[,mvar] <- relevel(mapping.sel[,mvar], metadata_variables[mvar, "baseline"])
		}
	} else if (metadata_variables[mvar, "type"] == "ordered") {
		if (metadata_variables[mvar, "baseline"] == "") {
			lvls <- unique(mapping.sel[,mvar])
		} else {
			lvls <- unlist(strsplit(metadata_variables[mvar, "baseline"], ","))
		}
		mapping.sel[,mvar] <- ordered(mapping.sel[,mvar], levels=lvls)
	} else if (metadata_variables[mvar, "type"] == "numeric") {
		mapping.sel[,mvar] <- as.numeric(as.character(mapping.sel[,mvar]))
	}
}
mapping.sel$Cohort <- droplevels(mapping.sel$Cohort)
sel <- intersect(rownames(mapping.sel), colnames(data))
data <- data[,sel]
mapping.sel <- mapping.sel[sel,]

## remove pathways with <50 reads, detected in less than 10 samples, scale to relative abundance
inds_to_remove <- which(rowSums(data) < 50)
data <- data[setdiff(1:nrow(data), inds_to_remove),]
inds_to_remove <- which(rowSums(data>0) < 10)
data <- data[setdiff(1:nrow(data), inds_to_remove),]
data.norm <- normalizeByCols(data)

## write files for FishTaco
genefamilies <- read.table(sprintf("%s/FishTaco/merged_genefamilies_KEGG.txt", out_dir), header=T, as.is=T, sep="\t", row.names=1, comment.char="")
colnames(genefamilies) <- gsub("PCMP_", "", gsub(".cat.filtered2_Abundance.RPKs", "", colnames(genefamilies)))
genefamilies <- genefamilies[, colnames(data.norm)]
inds_to_remove <- which(rowSums(genefamilies) < 1000)
genefamilies <- genefamilies[setdiff(1:nrow(genefamilies), inds_to_remove),]
inds_to_remove <- which(rowSums(genefamilies>100) < 10)
genefamilies <- genefamilies[setdiff(1:nrow(genefamilies), inds_to_remove),]
genefamilies.norm <- normalizeByCols(genefamilies)
write.table(genefamilies.norm, file=sprintf("%s/FishTaco/function_abundance_KEGG.txt", out_dir), quote=F, sep="\t", row.names=T, col.names=T)

## linear regression to identify differential pathways
res <- {}
for (f in rownames(data.norm)) {
	df <- melt(data.norm[f,]); rownames(df) <- as.character(df$variable)
	df <- merge(df, mapping.sel, by="row.names")
	mod <- lm(value ~ Cohort + BMI + Age + Sex, df)
	coef <- summary(mod)$coefficients
	res <- rbind(res, c(f, coef["CohortNASH",]))
}
res <- as.data.frame(res)
colnames(res) <- c("Pathway", "Estimate", "SE", "t", "pvalue")
res$pvalue <- as.numeric(as.character(res$pvalue))
res$padj <- p.adjust(res$pvalue)
write.table(res, file=sprintf("%s/Pathabundance.lm.txt", out_dir), quote=F, sep="\t", row.names=T, col.names=T)

## DESeq2 to identify differential pathways
dds <- DESeqDataSetFromMatrix(data, mapping.sel, design = ~Age + BMI + Sex + Cohort)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind=T)
res <- results(dds)
res <- res[order(res$pvalue),]
resSig <- subset(res, padj<0.2)
write.table(res, file=sprintf("%s/Pathabundance.DESeq2.txt", out_dir), quote=F, sep="\t", row.names=T, col.names=T)

## violin plots of significant hits
agg.melt <- melt(assay(vsd)[rownames(resSig),]); colnames(agg.melt) <- c("Pathway", "SampleID", "vsd"); agg.melt$SampleID <- as.character(agg.melt$SampleID); agg.melt$Pathway <- factor(as.character(agg.melt$Pathway), levels=rownames(resSig))
agg.melt$Cohort <- mapping.sel[agg.melt$SampleID, "Cohort"]
p <- ggplot(agg.melt, aes(x=Cohort, y=vsd, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~Pathway, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("VSD of pathway abundances (sighits by DESeq2)")) + coord_flip()
print(p)

## randomForest classification of Cohort (one-vs-one setup, X vs NMLBMI baseline)
agg <- data.norm
agg.melt.stored <- melt(as.matrix(agg), as.is=T); colnames(agg.melt.stored) <- c("pathway", "SampleID", "value")
set.seed(020818)
data.sel <- as.data.frame(t(agg))
data.sel$BMI <- mapping.sel[rownames(data.sel), "BMI"]
data.sel$Age <- mapping.sel[rownames(data.sel), "Age"]
data.sel$Sex <- mapping.sel[rownames(data.sel), "Sex"]
response <- droplevels(mapping.sel[rownames(data.sel), "Cohort"]); names(response) <- rownames(data.sel)
mvar_level <- "NASH"
level <- "Pathway"

## after running for the first time, COMMENT OUT THIS BLOCK ##
num_iter <- 100
ncores <- 20
out <- mclapply(1:num_iter, function (dummy) {
		importance(randomForest(x=data.sel, y=response, ntree=10000, importance=T), type=1, scale=F)
	}, mc.cores=ncores )
collated.importance <- do.call(cbind, out)
out <- mclapply(1:num_iter, function (dummy) {
		rfcv(trainx=data.sel, trainy=response, cv.fold=10, step=0.75)$error.cv
	}, mc.cores=ncores )
collated.cv <- do.call(cbind, out)

write.table(collated.importance, file=sprintf("%s/randomForest_Shotgun.%s.%s.importance.txt", out_dir, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=F)
write.table(collated.cv, file=sprintf("%s/randomForest_Shotgun.%s.%s.cv.txt", out_dir, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=F)
## END BLOCK TO COMMENT ##

collated.importance <- read.table(sprintf("%s/randomForest_Shotgun.%s.%s.importance.txt", out_dir, mvar_level, level), header=F, as.is=T, sep="\t", row.names=1, quote="")
collated.cv <- read.table(sprintf("%s/randomForest_Shotgun.%s.%s.cv.txt", out_dir, mvar_level, level), header=F, as.is=T, sep="\t", row.names=1, quote="")
importance.mean <- rowMeans(collated.importance)
importance.sd <- unlist(apply(collated.importance, 1, sd))
cv.mean <- rowMeans(collated.cv)
cv.sd <- unlist(apply(collated.cv, 1, sd))
inds <- order(importance.mean, decreasing=T)
inds <- inds[1:min(30, length(which(importance.mean[inds] > 0.001)))] # edit as appropriate

## after running for the first time, COMMENT OUT THIS BLOCK ##
# using a sparse model with N predictors
sparseRF <- randomForest(x=data.sel[, names(importance.mean[inds])], y=response, ntree=10000, importance=T)
save(sparseRF, file=sprintf("%s/randomForest_Shotgun.%s.%s.model", out_dir, mvar_level, level))
load(sprintf("%s/randomForest_Shotgun.%s.%s.model", out_dir, mvar_level, level))
## ROC and AUC of final sparseRF model
pred <- predict(sparseRF, type="prob")
pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true=response[rownames(pred)], stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
confusion_matrix <- table(pred_df[, c("true", "predicted")])
class_errors <- unlist(lapply(levels(mapping.sel$Cohort), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.sel$Cohort)
accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
mccvalue <- mcc(vec.pred, vec.true)
df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s) (accuracy = %.2f%%, MCC = %.4f)", level, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)
write.table(confusion_matrix, file=sprintf("%s/randomForest_Shotgun.%s.%s.confusion_matrix.txt", out_dir, mvar_level, level), quote=F, sep="\t", row.names=T, col.names=T)
## ROC analysis
if (nlevels(response)==2) {
	pred <- predict(sparseRF, type="prob")
	pred2 <- prediction(pred[,2], ordered(response))
	perf <- performance(pred2, "tpr", "fpr")
	perf.auc <- performance(pred2, "auc")
	print(plot(perf, main=sprintf("ROC %s %s (sparseRF final model)", mvar_level, level)) + text(x=0.8, y=0.2, label=sprintf("mean AUC=%.4g", unlist(perf.auc@y.values))))
}

## plotting - per-group sparse model
df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s", mvar_level, level)))
	
## plotting - per-group variables
df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
res <- read.table(sprintf("%s/Pathabundance.DESeq2.txt", out_dir), header=T, sep="\t", as.is=T, quote="", comment.char="")
df$log2FoldChange <- res[as.character(df$OTU), "log2FoldChange"]; df$padj <- res[as.character(df$OTU), "padj"]
df$OTU_string <- sprintf("%s (log2FC=%.4g, padj=%.4g)", rownames(df), df$log2FoldChange, df$padj)
df$sig <- ifelse(df$padj < 0.2, "sig", "ns"); df$sig[which(is.na(df$sig))] <- "ns"
p <- ggplot(df, aes(x=OTU, y=importance, label=OTU)) + geom_bar(position=position_dodge(), stat="identity", fill="#999999", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string, color=sig), size=3, hjust=0) + ggtitle(sprintf("%s - %s explanatory", mvar_level, level)) + scale_color_manual(values=cols.sig) + theme(axis.text.y=element_blank())
print(p)

## violin plots of relative abundance
agg.melt <- agg.melt.stored
agg.melt$Cohort <- mapping.sel[match(agg.melt$SampleID, mapping.sel$StudyID), "Cohort"]
agg.melt <- subset(agg.melt, pathway %in% rownames(df) & Cohort %in% c(mvar_level, "NMLBMI"))
agg.melt$pathway <- factor(agg.melt$pathway, levels=rownames(df))
p <- ggplot(agg.melt, aes(x=Cohort, y=value, color=Cohort)) + geom_violin() + geom_point() + facet_wrap(~pathway, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF pathways (%s)", mvar_level)) + coord_flip()
print(p)



