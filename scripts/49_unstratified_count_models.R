# Count-based follow-up of the complete unstratified correlation screen.
args <- commandArgs(trailingOnly = TRUE)
directory <- args[1]
B <- if (length(args) >= 3) as.integer(args[3]) else 1000L
output_directory <- if (length(args) >= 4) args[4] else directory
detection_fraction <- if (length(args) >= 5) as.numeric(args[5]) else .05
suppressPackageStartupMessages(library(Matrix))
source(args[2])
meta <- read.delim(file.path(directory, "metadata.tsv"))
genes <- readLines(file.path(directory, "genes.txt"))
peaks <- readLines(file.path(directory, "peaks.txt"))
score_genes <- readLines(file.path(directory, "programme_genes.txt"))
rna <- as(readMM(file.path(directory, "rna.mtx")), "CsparseMatrix")
atac <- as(readMM(file.path(directory, "atac.mtx")), "CsparseMatrix")
rownames(rna) <- genes
rownames(atac) <- peaks
candidates <- read.delim(file.path(output_directory, "count_candidate_pairs.tsv"))
models <- c("technical", "technical_coq", "technical_state", "technical_state_coq")

design <- function(m, model, gene, y, library_intercepts = FALSE) {
  terms <- cbind(log1p(m$total_rna_umi), log1p(m$total_open_peaks), m$percent_mito)
  if (grepl("coq", model)) terms <- cbind(terms, log1p(m$COQ8A_umi))
  if (grepl("state", model)) {
    score <- m$myogenesis_score
    if (gene %in% score_genes) {
      score <- (score * length(score_genes) - log1p(y / m$total_rna_umi * 10000)) /
        (length(score_genes) - 1)
    }
    terms <- cbind(terms, score)
  }
  terms <- scale(terms)
  terms[!is.finite(terms)] <- 0
  intercept <- if (library_intercepts) model.matrix(~factor(m$gsm)) else rep(1, nrow(m))
  cbind(intercept, terms)
}

robust_fit <- function(y, x) {
  fit <- glm.fit(x, y, family = poisson(), control = glm.control(maxit = 50))
  if (!fit$converged || any(!is.finite(fit$coefficients))) return(NULL)
  bread <- tryCatch(solve(crossprod(x, x * fit$fitted.values)), error = function(e) NULL)
  if (is.null(bread)) return(NULL)
  meat <- crossprod(x * (y - fit$fitted.values))
  se <- sqrt(diag(bread %*% meat %*% bread))[1]
  list(beta = fit$coefficients[1], se = se)
}

rows <- list()
k <- 0L
for (gsm in unique(meta$gsm)) {
  ids <- which(meta$gsm == gsm)
  m <- meta[ids, ]
  base_designs <- lapply(models, function(model) design(m, model, "", rep(0, nrow(m))))
  names(base_designs) <- models
  for (j in seq_len(nrow(candidates))) {
    gene <- candidates$gene[j]
    y <- as.numeric(rna[gene, ids])
    a <- as.numeric(atac[candidates$peak[j], ids] > 0)
    if (mean(y > 0) <= detection_fraction || mean(a > 0) <= detection_fraction || sum(y > 0) < 20 || sum(a == 0) < 20) next
    for (model in models) {
      cov <- if (gene %in% score_genes && grepl("state", model)) design(m, model, gene, y) else base_designs[[model]]
      fit <- robust_fit(y, cbind(a, cov))
      if (is.null(fit)) next
      k <- k + 1L
      rows[[k]] <- data.frame(gsm = gsm, peak = candidates$peak[j], gene = gene,
        model = model, n = length(ids), beta = fit$beta, se_robust = fit$se,
        p = 2 * pnorm(-abs(fit$beta / fit$se)))
    }
    if (j %% 200 == 0) {cat(gsm, j, "/", nrow(candidates), "pairs\n"); flush.console()}
  }
  cat(gsm, "completed\n"); flush.console()
}
per <- do.call(rbind, rows)
write.table(per, file.path(output_directory, "count_models_by_library.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
combined <- list()
k <- 0L
for (g in split(per, interaction(per$peak, per$gene, per$model, drop = TRUE))) {
  if (nrow(g) < 3) next
  weights <- 1 / g$se_robust^2
  beta <- weighted.mean(g$beta, weights)
  se <- sqrt(1 / sum(weights))
  loo <- sapply(seq_len(nrow(g)), function(i) weighted.mean(g$beta[-i], weights[-i]))
  k <- k + 1L
  combined[[k]] <- data.frame(peak = g$peak[1], gene = g$gene[1], model = g$model[1],
    n_libraries = nrow(g), beta = beta, RNA_ratio = exp(beta),
    CI_low = exp(beta - 1.96 * se), CI_high = exp(beta + 1.96 * se),
    p = 2 * pnorm(-abs(beta / se)), positive_libraries = sum(g$beta > 0),
    loo_RNA_ratio_min = exp(min(loo)), loo_RNA_ratio_max = exp(max(loo)))
}
result <- do.call(rbind, combined)
result$q_selected_followup <- ave(result$p, result$model, FUN = function(p) p.adjust(p, "BH"))
write.table(result, file.path(output_directory, "count_models_combined.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

target_file <- file.path(output_directory, "bootstrap_targets.tsv")
if (file.exists(target_file)) {
  targets <- read.delim(target_file)
  bootrows <- list()
  k <- 0L
  set.seed(208248)
  for (j in seq_len(nrow(targets))) {
    gene <- targets$gene[j]
    peak <- targets$peak[j]
    good <- unique(per$gsm[per$gene == gene & per$peak == peak])
    if (length(good) < 3) next
    ids <- which(meta$gsm %in% good)
    m <- meta[ids, ]
    y <- as.numeric(rna[gene, ids])
    a <- as.numeric(atac[peak, ids] > 0)
    for (model in c("technical", "technical_state_coq")) {
      x <- cbind(a, design(m, model, gene, y, library_intercepts = TRUE))
      fit <- glm.fit(x, y, family = poisson())
      groups <- split(seq_along(ids), m$gsm)
      values <- rep(NA_real_, B)
      for (b in seq_len(B)) {
        ix <- unlist(lapply(groups, function(v) sample(v, length(v), replace = TRUE)), use.names = FALSE)
        z <- glm.fit(x[ix, , drop = FALSE], y[ix], family = poisson(), start = fit$coefficients)
        if (z$converged) values[b] <- z$coefficients[1]
        if (b %% 250 == 0) {cat(gene, model, b, "/", B, "bootstraps\n"); flush.console()}
      }
      values <- values[is.finite(values)]
      ci <- quantile(values, c(.025, .975))
      k <- k + 1L
      bootrows[[k]] <- data.frame(gene = gene, peak = peak, model = model,
        n_nuclei = length(ids), n_libraries = length(good), RNA_ratio = exp(fit$coefficients[1]),
        CI_low = exp(ci[1]), CI_high = exp(ci[2]), p_boot = basic_p(fit$coefficients[1], values),
        B = length(values))
      write.table(do.call(rbind, bootrows), file.path(output_directory, "scent_style_bootstrap_followup.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
}
