multiplot <- function(..., plotlist=NULL, file, cols=1, rows=1) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
	
	i = 1
	while (i < numPlots) {
		numToPlot <- min(numPlots-i+1, cols*rows)
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(i, i+cols*rows-1), ncol = cols, nrow = rows, byrow=T)
		if (numToPlot==1) {
		  print(plots[[i]])
		} else {
		  # Set up the page
		  grid.newpage()
		  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		  # Make each plot, in the correct location
		  for (j in i:(i+numToPlot-1)) {
		    # Get the i,j matrix positions of the regions that contain this subplot
		    matchidx <- as.data.frame(which(layout == j, arr.ind = TRUE))
		    print(plots[[j]], vp = viewport(layout.pos.row = matchidx$row,
		                                    layout.pos.col = matchidx$col))
		  }
		}
		i <- i+numToPlot
  }
}

normalizeByRows <- function (df, rsum=1)
{
	while (any(abs((rowSums(df)-rsum))>1e-13)) {
		df <- rsum*(df / rowSums(df))
	}
	return(df)
}
normalizeByCols <- function (df, csum=1, level=NULL, delim="\\|")
{
	if (is.null(level)) {
		while (any(abs((colSums(df)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
			missing <- which(colSums(df)==0)
			df <- sweep(df, 2, colSums(df)/csum, "/")
			df[,missing] <- 0
		}
	} else {
	 tmp <- df
	 tmp$taxa <- rownames(tmp)
	 tmp$splitter <- factor(unlist(lapply(rownames(tmp), function(x) unlist(strsplit(x, delim))[level])))
	 names <- rownames(tmp)[order(tmp$splitter)]
	 tmp <- ddply(tmp, .(splitter), function(x) {
	 		x <- x[, setdiff(colnames(x), c("taxa", "splitter"))]
			while (any(abs((colSums(x)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
				x <- sweep(x, 2, colSums(x)/csum, "/")
			}
			x
		})
		rownames(tmp) <- names
		df <- tmp[, setdiff(colnames(tmp), "splitter")]
	}
	return(df)
}

renameLevelsWithCounts <- function(fvec, originalLevelsAsNames=FALSE) {
	tab <- table(fvec)
	retval <- sprintf("%s (n=%d)", fvec, tab[unlist(lapply(fvec, function(x) match(x, names(tab))))])
#	newlevels <- sprintf("%s (n=%d)", levels(fvec), tab[levels(fvec)])
	newlevels <- sprintf("%s (n=%d)", levels(fvec), tab[unlist(lapply(names(tab), function(x) which(levels(fvec)==x)))])
	retval <- factor(retval, levels=newlevels)
	if (originalLevelsAsNames) {
		names(retval) <- fvec
	}
	return(retval)
}

# get selected taxonomy
# DEFAULT: do not include entries that are not classified to the requested level (return NA instead)
getTaxonomy <- function(otus, tax_tab, level, na_str = c("unclassified", "unidentified", "NA", "", "uncultured"), includeUnclassified = TRUE) {
	# coerce taxonomyTable to data.frame (see https://github.com/joey711/phyloseq/issues/983), otherwise hangs as of phyloseq v1.26
	if (class(tax_tab) == "taxonomyTable") {
		tax_tab <- as.data.frame(tax_tab@.Data)
	}
	ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
	sel <- ranks[1:match(level, ranks)]
	inds <- apply(tax_tab[otus,sel], 1, function(x) max(which(!(x %in% na_str | is.na(x)))))
	if (includeUnclassified) {
		retval <- as.data.frame(tax_tab)[cbind(otus, ranks[inds])]
		inds.nmatch <- inds!=match(level, ranks)
		retval[inds.nmatch] <- paste(na_str[1], retval[inds.nmatch], sep=" ")
		if (level == "Species") {
			tmp <- as.data.frame(tax_tab)[cbind(otus, "Genus")]
			retval[!inds.nmatch] <- sprintf("%s %s", tmp[!inds.nmatch], retval[!inds.nmatch])
		}
	} else {
		retval <- as.data.frame(tax_tab)[cbind(otus, level)]
		if (level == "Species") {
			tmp <- as.data.frame(tax_tab)[cbind(otus, "Genus")]
			retval[!is.na(retval)] <- sprintf("%s %s", tmp[!is.na(retval)], retval[!is.na(retval)])
		}
	}
	retval <- gsub("\\[|\\]", "", retval)
	return(retval)
}


## modified rgcv function from 'spm' package to enable functionality of rfcv from 'randomForest'
## do cross-validation along feature selection axis in addition to sample-wise cross-validation
rgcv2 <- function (trainx, trainy, cv.fold = 10, scale = "log", step = 0.5, mtry = function(p) max(1, floor(sqrt(p))), num.trees = 500, min.node.size = NULL, num.threads = NULL, verbose = FALSE, recursive = FALSE, ...) {
	recursive.importance <- ifelse(recursive, "permutation", "none")
	classRF <- is.factor(trainy)
	n <- nrow(trainx)
	p <- ncol(trainx)
	if (scale == "log") {
		k <- floor(log(p, base = 1/step))
		n.var <- round(p * step^(0:(k - 1)))
		same <- diff(n.var) == 0
		if (any(same)) {
		  n.var <- n.var[-which(same)]
		}
		if (!1 %in% n.var) { 
		  n.var <- c(n.var, 1)
		}
	} else {
		n.var <- seq(from = p, to = 1, by = step)
	}
	k <- length(n.var)
	cv.pred <- vector(k, mode = "list")
	for (i in 1:k) cv.pred[[i]] <- trainy
	if (classRF) {
		f <- trainy
	} else {
		f <- cut(trainy, c(-Inf, stats::quantile(trainy, 1:4/5), Inf))
	}
	nlvl <- table(f)
	idx <- numeric(n)
	for (i in 1:length(nlvl)) {
		  idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold, 
		      length = nlvl[i]))
	}
	for (i in 1:cv.fold) {
		  data.dev <- trainx[idx != i, , drop = FALSE]
		  data.pred <- trainx[idx == i, , drop = FALSE]
		  response.dev <- trainy[idx != i]
		  all.rf <- ranger::ranger(x = data.dev, y=response.dev, mtry = mtry(p), 
		      num.trees = num.trees, min.node.size = min.node.size, 
		      num.threads = num.threads, verbose = verbose, importance="permutation")
		  cv.pred[[1]][idx == i] <- stats::predict(all.rf, data = data.pred)$predictions
		  impvar <- (1:p)[order(importance(all.rf, type = "permutation"), decreasing = TRUE)]
		  for (j in 2:k) {
		    imp.idx <- impvar[1:n.var[j]]
		    data.dev.sub <- data.dev[, imp.idx, drop=F]
		    sub.rf <- ranger::ranger(x = data.dev.sub, y=response.dev, mtry = mtry(n.var[j]), 
		      num.trees = num.trees, min.node.size = min.node.size, 
		      num.threads = num.threads, verbose = verbose, importance=recursive.importance)
		    cv.pred[[j]][idx == i] <- stats::predict(sub.rf, data = data.pred[, imp.idx, drop=F])$predictions
		    if (recursive) {
		      impvar <- (1:length(imp.idx))[order(importance(sub.rf, type = "permutation"), decreasing = TRUE)]
		    }
		  }
	}
	if (classRF) {
		error.cv <- sapply(cv.pred, function(x) mean(trainy != x))
	} else {
		error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
	}
	names(error.cv) <- names(cv.pred) <- n.var
	list(n.var = n.var, error.cv = error.cv, predicted = cv.pred)

}

