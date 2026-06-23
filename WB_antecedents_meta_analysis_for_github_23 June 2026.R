####################################
# Meta-analysis of Pearson's r     #
# Author: Elena Triantafillopoulou #
# Date: 23 June 2026               #
####################################


### 0) Clean workspace, install (if needed) and load packages --------

rm(list = ls())
graphics.off()
cat("\014")

set.seed(2026)   # for reproducibility

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
} else {
  setwd(dirname(sys.frame(1)$ofile))  # fallback when not in RStudio
}
cat("Working directory set to:", getwd(), "\n")

required_packages <- c(
  "metafor",       # core meta-analysis engine
  "readxl",        # read Excel files
  "dplyr",         # data wrangling
  "stringr",       # string cleaning
  "ggplot2",       # plots
  "clubSandwich",  # cluster-robust (CR2) standard errors
  "readr"          # write_csv
)

new_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(new_packages) > 0) {
  message("Installing missing packages: ", paste(new_packages, collapse = ", "))
  install.packages(new_packages, dependencies = TRUE)
} else {
  message("All required packages already installed.")
}

library(metafor)
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(clubSandwich)
library(readr)


### 1) Load dataset --------------------------------------------------

data_file <- "~/Desktop/WB_antecedents_meta_dataset_for_github_23 June 2026.xlsx"

# The categories to analyse — each becomes its own mini meta-analysis with its own forest plot, funnel plot, sensitivity analysis, etc

target_categories <- c(
  "Gender",
  "Psychological distress",
  "Bullied in adolescence",
  "Bullied at work",
  "Negative work climate",
  "Positive leadership",
  "Positive relationships at work",
  "Job control",
  "Job demands",
  "Role ambiguity and conflict",
  "Agreeableness",
  "Openness",
  "Conscientiousness",
  "Extraversion",
  "Neuroticism"
)
# NOTE: dataset_tag is now set per category inside the loop below.

run_beta_sensitivity    <- TRUE
# Convert beta -> r with the Peterson & Brown (2005) adjustment
# (r ~= beta + .05*lambda, lambda = 1 if beta >= 0, else 0). Set FALSE to
# treat the beta directly as r with no adjustment.
sens_use_peterson_brown <- FALSE
# Where the sampling variance comes from:
#   "n"    = from the sample size, exactly like the main analyses
#            (Fisher's z variance 1/(n-3)); sensitivity_p is then only
#            used as a printed consistency check.
#   "p_ci" = derived from the sensitivity_p column (p-value or 95% CI of
#            the beta); falls back to "n" for rows where it can't be parsed.
sens_variance_source    <- "n"

dat <- read_excel(data_file)

# Untouched copy of the raw spreadsheet for the beta sensitivity analysis
# (Section 13). Needed because the cleaning below drops rows with missing r
# and rows not flagged Include? == "Yes" — which would remove exactly the
# beta-only studies the sensitivity analysis is about.
dat_raw <- as.data.frame(dat)

cat("\nLoaded dataset with", nrow(dat), "rows and", ncol(dat), "columns.\n")
print(names(dat))


### 2) Clean / standardise columns -----------------------------------

# Column mapping
dat <- dat %>%
  rename(
    study_ID     = study_ID,
    Citation     = Citation,
    risk_factor  = risk_factor,
    r            = r,
    p_value      = p_value,
    n            = n,
    female       = female,
    time_lag     = time_lag,
    WB_construct = WB_construct,
    WB_measure   = WB_measure
  )

# ---- 2A: Clean p-values (handles "<", ">", "n.s.", numeric, HTML) ----
has_html_lt <- str_detect(dat$p_value, fixed("&lt;"))
has_html_gt <- str_detect(dat$p_value, fixed("&gt;"))

dat <- dat %>%
  mutate(
    p_value_numeric = suppressWarnings(as.numeric(p_value)),
    p_value_clean = case_when(
      str_detect(p_value, "<") | has_html_lt ~ p_value,
      str_detect(p_value, ">") | has_html_gt ~ p_value,
      str_detect(str_to_lower(p_value), "n.s")   ~ "p > 0.05",
      !is.na(p_value_numeric) ~ paste0("p = ", formatC(p_value_numeric, format = "f", digits = 4)),
      TRUE ~ as.character(p_value)
    )
  )

# ---- 2B: Type coercions ----
dat <- dat %>%
  mutate(
    study_ID     = as.character(study_ID),
    Citation     = as.character(Citation),
    n            = as.integer(n),
    female       = suppressWarnings(as.numeric(female)),
    time_lag     = suppressWarnings(as.numeric(time_lag)),
    WB_construct = as.factor(WB_construct),
    WB_measure   = as.factor(WB_measure)
  )

# Within-study effect index for multilevel nesting
dat <- dat %>%
  group_by(study_ID) %>%
  mutate(effect_id = row_number()) %>%
  ungroup()


### 3) Transform r → Fisher's z --------------------------------------

# ---- Clean r and n columns ----
# Converting cells to numeric, turning any non-numeric values to NA.
dat <- dat %>%
  mutate(
    r = suppressWarnings(as.numeric(as.character(r))),
    n = suppressWarnings(as.integer(as.character(n)))
  )

n_na_r <- sum(is.na(dat$r))
n_na_n <- sum(is.na(dat$n))
if (n_na_r > 0) cat("  WARNING:", n_na_r, "rows have non-numeric r and will be dropped.\n")
if (n_na_n > 0) cat("  WARNING:", n_na_n, "rows have non-numeric n and will be dropped.\n")

# Drop rows where r or n are missing (can't compute Fisher's z without both)
dat <- dat %>% filter(!is.na(r), !is.na(n))
cat("  Rows remaining after cleaning:", nrow(dat), "\n")

# Keep only studies flagged for inclusion (Include? == "Yes")
.inc_col <- names(dat)[tolower(gsub("[ _.?]", "", names(dat))) == "include"]
if (length(.inc_col)) {
  .keep <- tolower(trimws(as.character(dat[[.inc_col[1]]]))) %in% c("yes", "y", "true", "1")
  cat("  Include? filter ('", .inc_col[1], "'): keeping ", sum(.keep),
      " of ", nrow(dat), " rows.\n", sep = "")
  dat <- dat[.keep, , drop = FALSE]
} else {
  cat("  NOTE: no 'Include?' column found — no inclusion filter applied.\n")
}

# Bound r away from ±1 to avoid infinite z
dat <- dat %>%
  mutate(r = ifelse(r >=  0.9999,  0.9999,
             ifelse(r <= -0.9999, -0.9999, r)))

# Print all column names so we can confirm exact names before snapshotting
cat("\nColumn names in dat before escalc():\n")
print(names(dat))

# Detect the Meta_analysis column (case-insensitive, handles spaces/dots)
ma_col  <- names(dat)[tolower(names(dat)) %in% c("meta_analysis", "meta.analysis", "meta analysis")]
cit_col <- names(dat)[tolower(names(dat)) %in% c("citation", "citations")]

if (length(ma_col)  == 0) stop("Could not find Meta_analysis column — check name in Excel vs names printed above.")
if (length(cit_col) == 0) stop("Could not find Citation column — check name in Excel vs names printed above.")

cat("  Using Meta_analysis column:", ma_col[1],  "\n")
cat("  Using Citation column:     ", cit_col[1], "\n")

# Snapshot BEFORE escalc() so lengths are guaranteed to match
meta_analysis_snap <- dat[[ ma_col[1]  ]]
citation_snap      <- dat[[ cit_col[1] ]]

escalc_dat <- escalc(measure = "ZCOR", ri = r, ni = n, data = dat)

# Reattach from snapshot (same length as escalc_dat, no row-mismatch risk)
escalc_dat$Meta_analysis <- as.character(meta_analysis_snap)
escalc_dat$Citation      <- as.character(citation_snap)
cat("\nMeta_analysis groups found:", paste(unique(escalc_dat$Meta_analysis), collapse = ", "), "\n")

# ---- Keep a master copy of the FULL effect-size dataset ----
# Each category below is filtered from this; the master is never modified.
escalc_full <- escalc_dat

# ---- Diagnostic: categories in the file vs the target list ----
found_cats <- sort(unique(trimws(as.character(escalc_full$Meta_analysis))))
cat("\nMeta_analysis categories present in the file:\n"); print(found_cats)
.miss  <- target_categories[!tolower(trimws(target_categories)) %in% tolower(found_cats)]
.extra <- found_cats[!tolower(found_cats) %in% tolower(trimws(target_categories))]
if (length(.miss))
  cat("\n  WARNING — no rows found for these target categories (check spelling in Excel):\n   ",
      paste(.miss, collapse = " | "), "\n")
if (length(.extra))
  cat("  NOTE — file has categories not in your target list (they will be skipped):\n   ",
      paste(.extra, collapse = " | "), "\n")


### ---- Repeated-sample clustering (Citation + Cohort) --------
# Effects are treated as non-independent (same cluster) if they come from the
# same Citation OR share cohort acronyms

make_cluster_id <- function(citation, cohort, ignore_tokens = character(0)) {
  n <- length(citation)
  if (n == 0) return(integer(0))
  parent <- seq_len(n)
  find <- function(x) {
    root <- x
    while (parent[root] != root) root <- parent[root]
    while (parent[x] != root) { nxt <- parent[x]; parent[x] <<- root; x <- nxt }
    root
  }
  unite <- function(a, b) { ra <- find(a); rb <- find(b); if (ra != rb) parent[ra] <<- rb }

  cit <- tolower(trimws(as.character(citation)))
  for (key in unique(cit[!is.na(cit) & cit != ""])) {
    idx <- which(cit == key)
    if (length(idx) > 1) for (j in 2:length(idx)) unite(idx[1], idx[j])
  }

  toks <- strsplit(tolower(trimws(as.character(cohort))), "\\s*,\\s*")
  tok_map <- list()
  for (i in seq_len(n)) for (t in toks[[i]]) {
    t <- trimws(t)
    if (is.na(t) || t == "" || t %in% c("na", "n/a", "none", "-")) next
    if (t %in% ignore_tokens) next   # registry tokens don't merge studies (policy-driven)
    tok_map[[t]] <- c(tok_map[[t]], i)
  }
  for (rows in tok_map) if (length(rows) > 1) for (j in 2:length(rows)) unite(rows[1], rows[j])
  as.integer(factor(vapply(seq_len(n), find, integer(1))))
}

# ---- Multilevel (RVE) Egger-type test for funnel asymmetry --------------
if (!requireNamespace("clubSandwich", quietly = TRUE)) install.packages("clubSandwich")

ml_egger <- function(dat, label = "") {
  dat <- dat[is.finite(dat$yi) & is.finite(dat$vi) & dat$vi > 0, ]
  k <- nrow(dat); ncl <- length(unique(dat$cluster_id))
  if (k < 10 || ncl < 5) {
    cat(sprintf("  [%s] Egger test skipped: k = %d effects / %d clusters — too few to be reliable\n",
                label, k, ncl))
    return(invisible(NULL))
  }
  dat$sei <- sqrt(dat$vi)
  m  <- rma.mv(yi, vi, mods = ~ sei,
               random = ~ 1 | cluster_id/effect_id,
               data = dat, method = "REML", test = "t")
  mr <- robust(m, cluster = dat$cluster_id, clubSandwich = TRUE)   # CR2 + Satterthwaite df
  i  <- which(rownames(mr$b) == "sei")
  cat(sprintf("  [%s] Multilevel Egger (CR2): slope = %.3f, 95%% CI [%.3f, %.3f], t(%.1f) = %.2f, p = %.3f\n",
              label, mr$b[i], mr$ci.lb[i], mr$ci.ub[i], mr$dfs[i], mr$zval[i], mr$pval[i]))
  data.frame(analysis = label, k_effects = k, n_clusters = ncl,
             slope = as.numeric(mr$b[i]), ci_lb = mr$ci.lb[i], ci_ub = mr$ci.ub[i],
             t = mr$zval[i], df = mr$dfs[i], p = mr$pval[i],
             intercept_z = as.numeric(mr$b[1]),          # PET-style adjusted effect (Fisher's z)
             intercept_r = tanh(as.numeric(mr$b[1])),    # back-transformed to r
             stringsAsFactors = FALSE)
}

# --- Registry rules (applies to ALL clustering: models, forest, total N) -----

registry_cohorts <- c("sn", "ss")   # case-insensitive
registry_policy  <- "separate"      # "separate" = don't merge on registry; "merge" = do
ignore_registry  <- if (identical(registry_policy, "separate")) registry_cohorts else character(0)

### ---- Accent-aware sort key (so e.g. "Ågotnes" sorts with the A's) -----
# Default collation pushed Å/Ø/Æ-initial names to the end
ascii_key <- function(x) {
  x <- as.character(x)
  x <- gsub("[ÅÄÆ]", "A", x); x <- gsub("[åäæ]", "a", x)
  x <- gsub("[ØÖ]", "O", x);  x <- gsub("[øö]", "o", x)
  x <- gsub("Ü", "U", x); x <- gsub("ü", "u", x)
  k <- iconv(x, to = "ASCII//TRANSLIT", sub = "")
  k[is.na(k) | k == ""] <- x[is.na(k) | k == ""]
  toupper(trimws(gsub("[^A-Za-z0-9 ]", "", k)))
}

if (!exists("col_diamond")) col_diamond <- "#C0392B"

# Map Measure-type codes to symbols 
measure_symbols <- c("SL"   = "\u25B3",   # △ self-labelling items
                     "Q"    = "\u25CF",   # ● validated questionnaire
                     "Both" = "\u25D1")   # ◑ both
recode_measure <- function(x) {
  x   <- trimws(as.character(x))
  out <- measure_symbols[x]
  out[is.na(out)] <- x[is.na(out)]        # keep anything not in the lookup
  unname(out)
}

### ---- Grouped (combined) forest figure -----------------------------
# Renders several categories as stacked subgroups in ONE condensed figure,
# matching the manuscript layout: study labels on the left, an "analysis
# specifications" block (Time lag | Measure | Antecedent | Sample notes) to the
# left of the bars, the forest in the middle, and the correlation (95% CI) on
# the right. Each subgroup gets its own pooled diamond + heterogeneity line.
make_grouped_forest <- function(fdat, subgroup_order, out_tag, title = NULL,
                                condensed = TRUE,
                                xlab_text = "Correlation (r)",
                                ci_lab    = "Correlation (95% CI)") {
  if (!requireNamespace("metafor", quietly = TRUE)) return(invisible(NULL))
  has_st <- requireNamespace("showtext", quietly = TRUE)
  if (has_st) suppressMessages(library(showtext))
  cal    <- isTRUE(get0("calibri_loaded", ifnotfound = FALSE))
  col_d  <- if (exists("col_diamond")) col_diamond else "#C0392B"
  fdat <- as.data.frame(fdat)

  mg   <- tolower(trimws(as.character(fdat$Meta_analysis)))
  subs <- subgroup_order[tolower(trimws(subgroup_order)) %in% unique(mg)]
  if (length(subs) == 0) {
    cat("  (grouped forest '", out_tag, "': no matching rows)\n", sep = ""); return(invisible(NULL))
  }

  # Helper columns
  cohort_vals <- if (length(cohort_col)) as.character(fdat[[cohort_col[1]]]) else rep("", nrow(fdat))
  fdat$cluster_id <- make_cluster_id(fdat$Citation, cohort_vals, ignore_registry)
  fdat$effect_id  <- seq_len(nrow(fdat))
  .sncol <- names(fdat)[tolower(gsub("[ _.]", "", names(fdat))) == "samplenotes"]
  fdat$.sn   <- if (length(.sncol)) as.character(fdat[[.sncol[1]]]) else ""
  sym_SL <- if (!is.na(measure_symbols["SL"])) unname(measure_symbols["SL"]) else "(SL)"
  sym_Q  <- if (!is.na(measure_symbols["Q"]))  unname(measure_symbols["Q"])  else "(Q)"
  fdat$.sn <- gsub("\\(\\s*SL\\s*\\)", sym_SL, fdat$.sn, ignore.case = TRUE)
  fdat$.sn <- gsub("\\(\\s*Q\\s*\\)",  sym_Q,  fdat$.sn, ignore.case = TRUE)
  .adc <- names(fdat)[tolower(gsub("[ _.]", "", names(fdat))) == "riskfactordisplay"]
  if (!length(.adc)) .adc <- names(fdat)[tolower(gsub("[ _.]", "", names(fdat))) == "riskfactor"]
  fdat$.ante <- if (length(.adc)) as.character(fdat[[.adc[1]]]) else ""
  fdat$.wb   <- recode_measure(fdat$WB_measure)
  fdat$.tl   <- as.character(fdat$time_lag)
  for (cc in c(".sn", ".ante", ".wb"))
    fdat[[cc]][is.na(fdat[[cc]]) | fdat[[cc]] == "NA"] <- ""

  # Fold the Antecedent into Sample notes (semicolon-joined); drop standalone column
  .a <- trimws(fdat$.ante); .s <- trimws(fdat$.sn)
  fdat$.sn <- ifelse(.a != "" & .s != "", paste0(.a, "; ", .s), paste0(.a, .s))

  fdat$yi_r <- tanh(fdat$yi)
  fdat$vi_r <- (1 - fdat$yi_r^2)^2 * fdat$vi

  # Order rows: by subgroup order, then accent-aware citation, then time lag
  fdat <- fdat[order(match(mg, tolower(trimws(subs))),
                     ascii_key(fdat$Citation),
                     suppressWarnings(as.numeric(fdat$time_lag))), ]
  mg <- tolower(trimws(as.character(fdat$Meta_analysis)))

  # Per-subgroup random-effects fits
  fits <- lapply(subs, function(g) {
    d <- fdat[mg == tolower(trimws(g)), ]
    tryCatch(
      if (any(duplicated(d$cluster_id)))
        metafor::rma.mv(yi, vi, random = ~ 1 | cluster_id/effect_id, data = d, method = "REML")
      else metafor::rma(yi, vi, data = d, method = "REML"),
      error = function(e) NULL)
  })
  names(fits) <- subs

  # Row layout (metafor rows increase upward; build from the top down).
  gap_after_studies <- 1     # blank row between the studies and the diamond
  het_rows          <- 3     # n/k/N line + heterogeneity printed on three lines
  gap_between       <- 2.5   # blank space between subgroups (tightened)
  bottom_pad        <- 3     # clearance above the x-axis (stops overlap)
  block_rows <- sapply(subs, function(g)
    1 + sum(mg == tolower(trimws(g))) + gap_after_studies + 1 + het_rows + gap_between)
  total_rows <- sum(block_rows) + bottom_pad

  cur <- total_rows
  study_rows <- rep(NA_real_, nrow(fdat))
  header_y <- diamond_y <- het_y <- list()
  for (g in subs) {
    header_y[[g]] <- cur; cur <- cur - 1
    idx <- which(mg == tolower(trimws(g))); k <- length(idx)
    study_rows[idx] <- seq(cur, cur - k + 1); cur <- cur - k - gap_after_studies
    diamond_y[[g]] <- cur; cur <- cur - 1
    het_y[[g]]     <- cur; cur <- cur - het_rows - gap_between
  }
  fdat$.row <- study_rows

  # Geometry
  r_lb <- tanh(fdat$yi - 1.96 * sqrt(fdat$vi)); r_ub <- tanh(fdat$yi + 1.96 * sqrt(fdat$vi))
  r_min <- min(r_lb, na.rm = TRUE); r_max <- max(r_ub, na.rm = TRUE)
  step  <- if ((r_max - r_min) > 1.4) 0.25 else if ((r_max - r_min) > 1.1) 0.20 else 0.10
  t_lo  <- min(0, floor(r_min / step) * step)   # include 0 but don't waste space below the data
  t_hi  <- ceiling(r_max / step) * step
  ticks <- seq(t_lo, t_hi, by = step)
  pad <- 0.04
  bar_left <- r_min - pad; bar_right <- r_max + pad
  # Spec columns (Antecedent folded into Sample notes); tightened to cut whitespace
  x_study <- bar_left - 1.30   # study names  (nudge if citations clip or overlap Time lag)
  x_tl    <- bar_left - 0.85   # Time lag
  x_meas  <- bar_left - 0.70   # Measure type
  x_sn    <- bar_left - 0.55   # Sample notes (now also holds the antecedent; needs room before bars)
  x_ci    <- bar_right + 0.05
  xlim <- c(x_study - 0.05, x_ci + 0.95)
  cex0  <- if (condensed) 0.66 else 0.90
  pdf_w <- 10
  pdf_h <- max(4, total_rows * (if (condensed) 0.11 else 0.17) + 1.4)

  draw <- function() {
    if (cal) par(family = "Calibri")
    op <- par(mar = c(3.2, 0.5, 1.6, 0.5), xpd = NA, mgp = c(2.1, 0.4, 0), tcl = -0.25)
    metafor::forest(x = fdat$yi_r, vi = fdat$vi_r, rows = fdat$.row,
      slab = fdat$Citation, ylim = c(0, total_rows + 4), xlim = xlim,
      at = ticks, xlab = "",
      ilab = data.frame(fdat$.tl, fdat$.wb, fdat$.sn),
      ilab.xpos = c(x_tl, x_meas, x_sn), ilab.pos = 4,
      psize = 1, cex = cex0, lwd = 0.6, header = FALSE, annotate = FALSE,
      efac = c(1, 0), mlab = "")

    mtext(xlab_text, side = 1, line = 1.6, cex = cex0 * 0.8)
    
    text(x_ci, fdat$.row,
         paste0(formatC(tanh(fdat$yi), 2, format = "f"), " [",
                formatC(r_lb, 2, format = "f"), ", ",
                formatC(r_ub, 2, format = "f"), "]"),
         pos = 4, cex = cex0 * 0.95)

    text(x_tl,   total_rows + 0.4, "Time lag & measure",    pos = 4, cex = cex0 * 0.95)
    text(x_sn,   total_rows + 0.4, "Antecedent notes",      pos = 4, cex = cex0 * 0.95)
    text(x_ci,   total_rows + 0.4, "Correlation (95% CI)", pos = 4, cex = cex0 * 0.95)
    if (!is.null(title)) text(xlim[1], total_rows + 3.5, title, pos = 4, cex = cex0 * 1.2)
    par(font = 1)

    for (g in subs) {
      par(font = 2); text(xlim[1], header_y[[g]], g, pos = 4, cex = cex0); par(font = 1)
      fit <- fits[[g]]; if (is.null(fit)) next
      # Contributing studies (n), effect sizes (k) and deduplicated sample size (N).
      # N sums one (largest) n per cluster_id, so repeated samples aren't counted twice.
      dg  <- fdat[mg == tolower(trimws(g)), ]
      .nv <- suppressWarnings(as.numeric(dg$n))
      N_g <- sum(tapply(.nv, dg$cluster_id, function(z) max(z, na.rm = TRUE)), na.rm = TRUE)
      nkN_lab <- sprintf("n = %d, k = %d, N = %s",
                         length(unique(dg$Citation)), nrow(dg),
                         formatC(N_g, format = "d", big.mark = ","))
      est <- tanh(as.numeric(fit$b[1, 1])); lb <- tanh(fit$ci.lb); ub <- tanh(fit$ci.ub)
      metafor::addpoly(x = est, ci.lb = lb, ci.ub = ub, row = diamond_y[[g]],
        mlab = "", cex = cex0, efac = 0.8, col = col_d, border = col_d, annotate = FALSE)
      text(x_ci, diamond_y[[g]],
           paste0(formatC(est, 2, format = "f"), " [",
                  formatC(lb, 2, format = "f"), ", ", formatC(ub, 2, format = "f"), "]"),
           pos = 4, cex = cex0 * 0.95)
      if (inherits(fit, "rma.uni")) {
        l1 <- bquote(paste("Heterogeneity: ", tau^2, " = ", .(round(fit$tau2, 2)),
                           ",  ", I^2, " = ", .(paste0(round(fit$I2, 1), "%"))))
        l2 <- bquote(paste("Q(", .(fit$k - 1), ") = ", .(round(fit$QE, 2)),
                           ", p = ", .(formatC(fit$QEp, digits = 2, format = "f"))))
      } else {
        s2 <- sum(fit$sigma2)
        tv <- s2 + median(fdat$vi[mg == tolower(trimws(g))], na.rm = TRUE)
        l1 <- bquote(paste("Heterogeneity: ", tau^2, " between = ", .(round(fit$sigma2[1], 2)),
                           ", within = ", .(round(fit$sigma2[2], 2))))
        l2 <- paste0("I\u00B2 total = ", round(100 * s2 / tv, 1), "%")
      }
      text(xlim[1], het_y[[g]],       nkN_lab, pos = 4, cex = cex0 * 0.82)
      text(xlim[1], het_y[[g]] - 0.9, l1,      pos = 4, cex = cex0 * 0.82)
      text(xlim[1], het_y[[g]] - 1.8, l2,      pos = 4, cex = cex0 * 0.82)
    }
    
    # Symbol key (bottom-left) — reuses your measure_symbols so it stays in sync
    text(xlim[1], 1.9, paste0(measure_symbols["Q"],  "  Validated questionnaire (Q)"),
         pos = 4, cex = cex0 * 0.8)
    text(xlim[1], 1.0, paste0(measure_symbols["SL"], "  Self-labelling (SL)"),
         pos = 4, cex = cex0 * 0.8)
    
    par(op)
  }

  pdf(paste0("forest_grouped_", out_tag, ".pdf"), width = pdf_w, height = pdf_h, pointsize = 10)
  if (has_st) { showtext_auto(); showtext_opts(dpi = 300) }
  draw()
  if (has_st) showtext_auto(FALSE)
  dev.off()
  png(paste0("forest_grouped_", out_tag, ".png"),
      width = pdf_w, height = pdf_h, units = "in", res = 300, pointsize = 10)
  if (has_st) { showtext_auto(); showtext_opts(dpi = 300) }
  draw()
  if (has_st) showtext_auto(FALSE)
  dev.off()
  cat("  Saved forest_grouped_", out_tag, ".png\n", sep = "")
#  cat("  Saved forest_grouped_", out_tag, ".pdf\n", sep = "")
  invisible(NULL)
}

# Which categories to combine into single grouped figures (edit freely).
# Each entry = one figure; the vector lists the subgroups it stacks, in order.
forest_groups <- list(
  "Prior bullying victimisation" = c("Bullied in adolescence", "Bullied at work"),
  "Individual factors"           = c("Gender", "Psychological distress"),
  "Personality"                  = c("Agreeableness", "Openness", "Conscientiousness",
                                     "Extraversion", "Neuroticism"),
  "Work social environment"      = c("Negative work climate", "Positive leadership",
                                     "Positive relationships at work"),
  "Job and role characteristics" = c("Job control", "Job demands", "Role ambiguity and conflict")
)

# ---- Detect optional columns by normalised name (robust to spaces/case) ----
.norm_names <- tolower(gsub("[ _.]", "", names(escalc_full)))
cohort_col  <- names(escalc_full)[.norm_names == "cohort"]
rob_col     <- names(escalc_full)[.norm_names == "riskofbias"]
ante_col    <- names(escalc_full)[.norm_names == "riskfactordisplay"]
if (length(cohort_col) == 0) cat("\n  NOTE: no 'Cohort' column found — clustering falls back to Citation only.\n")
if (length(rob_col)    == 0) cat("  NOTE: no 'risk_of_bias' column found — risk-of-bias moderation will be skipped.\n")
if (length(ante_col)   == 0) cat("  NOTE: no 'risk_factor_display' column found — forest 'Antecedent' column falls back to risk_factor.\n")

### ---- Total N contributing to the meta-analyses --------------------
# Uses escalc_full, i.e. exactly the rows that enter the meta-analyses
# (Include? == Yes, valid r and n). Each study (Citation) is counted once
# at its largest n; studies sharing a Cohort token are then merged and only
# the largest study N in that group counts. Two exceptions:
#  * registry tokens (SN/SS): under registry_policy "separate" (the default) a
#    shared registry does NOT merge studies, so each study's n is counted and the
#    total is an over-estimate ("approximately"); under "merge" they count once.
#  * "Only women"/"Only men" in Sample_notes marks disjoint subsamples,
#    which are SUMMED within that study instead of maxed.

tryCatch({
  .nN   <- suppressWarnings(as.numeric(escalc_full$n))
  .citN <- trimws(as.character(escalc_full$Citation))
  .cohN <- if (length(cohort_col)) as.character(escalc_full[[cohort_col[1]]])
  else rep("", nrow(escalc_full))
  .cohN[is.na(.cohN)] <- ""
  .snc  <- names(escalc_full)[tolower(gsub("[ _.]", "", names(escalc_full))) == "samplenotes"]
  .noteN <- if (length(.snc)) tolower(trimws(as.character(escalc_full[[.snc[1]]])))
  else rep("", nrow(escalc_full))
  .noteN[is.na(.noteN)] <- ""
  .gtag <- ifelse(grepl("only wom|only men", .noteN), .noteN, "")
  
  # Keep all cohort tokens here; registry policy is applied centrally by
  # make_cluster_id(..., ignore_registry) below, so the count matches the models.
  .coh_clean <- vapply(strsplit(tolower(.cohN), "\\s*,\\s*"), function(tk) {
    tk <- trimws(tk)
    tk <- tk[!is.na(tk) & tk != ""]
    paste(tk, collapse = ", ")
  }, character(1))
  
  .stu <- data.frame(cit = .citN, n = .nN, g = .gtag, coh = .coh_clean,
                     stringsAsFactors = FALSE)
  .stu <- .stu[!is.na(.stu$n) & .stu$cit != "", , drop = FALSE]
  
  # One N per study: largest n per Citation, but disjoint gender
  # subsamples are summed (largest n per gender, then added together).
  per_study <- do.call(rbind, lapply(split(.stu, .stu$cit), function(d) {
    n_full <- if (any(d$g == "")) max(d$n[d$g == ""]) else NA_real_
    n_gsum <- if (any(d$g != "")) sum(tapply(d$n[d$g != ""], d$g[d$g != ""], max)) else NA_real_
    data.frame(cit = d$cit[1],
               n_study = max(c(n_full, n_gsum), na.rm = TRUE),
               coh = paste(unique(d$coh[d$coh != ""]), collapse = ", "),
               stringsAsFactors = FALSE)
  }))
  
  per_study$grp <- make_cluster_id(per_study$cit, per_study$coh, ignore_registry)
  n_by_grp <- tapply(per_study$n_study, per_study$grp, max)
  total_N  <- sum(n_by_grp)
  
  cat("\n=== Total N contributing to the meta-analyses ===\n")
  cat("  Contributing studies:", nrow(per_study),
      "| independent sample groups:", length(n_by_grp), "\n")
  .merged <- split(per_study$cit, per_study$grp)
  .merged <- .merged[vapply(.merged, length, integer(1)) > 1]
  if (length(.merged)) {
    cat("  Counted once (shared cohort) — largest study N used:\n")
    for (m in .merged) cat("    -", paste(m, collapse = " + "), "\n")
  }
  .gstud <- unique(.stu$cit[.stu$g != ""])
  if (length(.gstud))
    cat("  Disjoint gender subsamples summed within:", paste(.gstud, collapse = ", "), "\n")
  cat("  Registry policy:", registry_policy,
      if (identical(registry_policy, "merge"))
        "(studies sharing SN/SS counted once -> lower bound)"
      else "(SN/SS studies counted separately -> over-estimate)", "\n")
  cat(if (identical(registry_policy, "merge")) "  TOTAL N (at least) ="
      else "  TOTAL N (estimate) =",
      formatC(total_N, format = "d", big.mark = ","), "\n")
}, error = function(e) cat("  (total-N calculation skipped:", conditionMessage(e), ")\n"))

# ---- Global accumulator for the combined cross-category moderation table ----
overall_table_rows <- list()

# ---- Global accumulator for each category's MAIN pooled result ----
# (used by Section 13 to build the main vs beta-sensitivity comparison table)
overall_pooled <- list()


### ====================================================================
### LOOP: run one full mini meta-analysis per Meta_analysis category.
### Everything from here to the end (model, forest plot, funnel plot,
### leave-one-out sensitivity, moderation analyses, exports) runs once
### per category, writing output files tagged with that category name.
### ====================================================================

egger_all <- list()

for (current_category in target_categories) {

  # Filter the master dataset to this category (case-insensitive, trimmed)
  escalc_dat <- escalc_full[
    tolower(trimws(as.character(escalc_full$Meta_analysis))) ==
      tolower(trimws(current_category)), , drop = FALSE]

  # Filename-safe tag used on every output file for this category
  dataset_tag <- gsub("(^_|_$)", "",
                      gsub("[^A-Za-z0-9]+", "_", trimws(current_category)))

  cat("\n\n##########################################################\n")
  cat("###  MINI META-ANALYSIS:", current_category,
      "(k =", nrow(escalc_dat), "effect sizes)\n")
  cat("##########################################################\n")

  if (nrow(escalc_dat) < 2) {
    cat("  Skipping '", current_category,
        "' — needs at least 2 effect sizes.\n", sep = "")
    next
  }

  # Reset the per-category moderation collector
  moderation_results_all <- list()

  # ---- Build repeated-sample clusters for THIS category ----
  .cohort_vals <- if (length(cohort_col))
    as.character(escalc_dat[[cohort_col[1]]]) else rep("", nrow(escalc_dat))
  escalc_dat$cluster_id <- make_cluster_id(escalc_dat$Citation, .cohort_vals, ignore_registry)
  # Unique within-category effect id (nests under cluster_id in the models)
  escalc_dat$effect_id  <- seq_len(nrow(escalc_dat))
  cat("  Repeated-sample clusters:", length(unique(escalc_dat$cluster_id)),
      "| distinct citations:", length(unique(escalc_dat$Citation)),
      "| effects:", nrow(escalc_dat), "\n")

  # Wrap each category so one failure does not stop the remaining ones
  tryCatch({

### 4) Primary meta-analysis (auto multilevel if duplicates) ----------

multi_es <- any(duplicated(escalc_dat$cluster_id))

if (multi_es) {
  cat("\n=== Multilevel random-effects model (Fisher's z) ===\n")
  res_main <- rma.mv(yi, vi,
                     random = ~ 1 | cluster_id/effect_id,
                     data   = escalc_dat,
                     method = "REML")
} else {
  cat("\n=== Single-level random-effects model (Fisher's z) ===\n")
  res_main <- rma(yi, vi, data = escalc_dat, method = "REML")
}

print(res_main)

# Back-transform pooled estimate
pooled_z    <- as.numeric(res_main$b[1, 1])
pooled_r    <- tanh(pooled_z)
pooled_r_ci <- tanh(c(res_main$ci.lb, res_main$ci.ub))

cat("\nPooled r =", round(pooled_r, 3),
    " [", round(pooled_r_ci[1], 3), ",", round(pooled_r_ci[2], 3), "]\n")

# Store this category's pooled result for the Section 13 comparison table
overall_pooled[[current_category]] <- data.frame(
  Category = current_category,
  k        = res_main$k,
  r        = pooled_r,
  lb       = pooled_r_ci[1],
  ub       = pooled_r_ci[2],
  stringsAsFactors = FALSE
)


### 4B) Heterogeneity reporting --------------------------------------

if (inherits(res_main, "rma.uni")) {
  cat("\nHeterogeneity (single-level):\n")
  cat("  Q  =", round(res_main$QE, 2), "  p =", formatC(res_main$QEp, digits = 3, format = "f"), "\n")
  cat("  τ² =", round(res_main$tau2, 5), "\n")
  cat("  I² =", round(res_main$I2, 1), "%\n")
} else {
  cat("\nHeterogeneity (multilevel):\n")
  cat("  Variance components (τ²):\n")
  print(setNames(round(res_main$sigma2, 5), c("between-study", "within-study")))

  i2_out <- tryCatch(metafor::i2(res_main), error = function(e) NULL)
  if (!is.null(i2_out)) {
    cat("\n  I² by level:\n")
    print(round(i2_out, 1))
  } else {
    cat("\n  (I² decomposition not available in your metafor version.)\n")
  }

  cat("\nCluster-robust (CR2) test for pooled effect:\n")
  print(clubSandwich::coef_test(res_main, vcov = "CR2",
                                cluster = escalc_dat$cluster_id))
}


### 5) Grouped forest plot — one panel per Meta_analysis subgroup -----
#
#  Matches the reference style:
#   * Bold subgroup headers flush left
#   * Navy blue squares (individual studies), red diamonds (pooled)
#   * Heterogeneity stats + Q-test printed on the diamond row
#   * Study label left, r [95% CI] right
#   * Single vertical reference line at zero
#   * Separate rma (or rma.mv) fit per subgroup
#
#  The plot is built by manually assigning row positions so every
#  subgroup sits in its own clearly separated block, with a blank
#  spacer row between groups.

# ---- Colours (matching reference image) ----
col_study   <- "#1F3F6B"   # navy blue for study squares
col_diamond <- "#2166AC"   # blue for pooled diamonds

# ---- Identify subgroups (in order of appearance) ----
# Clean the column first: trim whitespace and fix any capitalisation
# inconsistencies that could cause the same group to appear twice
escalc_dat$Meta_analysis <- trimws(as.character(escalc_dat$Meta_analysis))

subgroups <- unique(escalc_dat$Meta_analysis)
n_groups  <- length(subgroups)

cat("\nNumber of Meta_analysis groups found:", n_groups, "\n")
cat("Group names:", paste(subgroups, collapse = ", "), "\n")

# Each iteration analyses a single category, so exactly one group is expected.
if (n_groups < 1) {
  warning("No Meta_analysis group found for '", current_category, "' — skipping.")
}

# ---- Fit a separate model per subgroup ----
# Uses multilevel (rma.mv) if any study_ID is duplicated within the
# subgroup, otherwise single-level rma — mirrors the primary analysis.
fit_subgroup <- function(df) {
  dup <- any(duplicated(df$cluster_id))
  if (dup) {
    tryCatch(
      rma.mv(yi, vi, random = ~ 1 | cluster_id/effect_id,
             data = df, method = "REML"),
      error = function(e) NULL
    )
  } else {
    tryCatch(
      rma(yi, vi, data = df, method = "REML"),
      error = function(e) NULL
    )
  }
}

sub_fits   <- lapply(subgroups, function(g)
                fit_subgroup(escalc_dat[escalc_dat$Meta_analysis == g, ]))
names(sub_fits) <- subgroups

for (g in subgroups) {
  cat("\n=== Subgroup model:", g, "===\n")
  print(sub_fits[[g]])
}

# ---- Build row-position table ----
# metafor rows count from the BOTTOM upward (row 1 = bottom of plot).
# We build positions bottom-up: last group first, then work upward.
# Each group block (bottom to top):
#   [spacer rows]        <- gap between groups
#   [2 hetero text rows] <- heterogeneity stats
#   [1 diamond row]      <- pooled effect
#   [k_g study rows]     <- individual studies
#   [1 header gap]       <- bold group label (drawn as text)

spacer      <- 2
hetero_rows <- 2

rows_per_group <- sapply(subgroups, function(g)
  sum(escalc_dat$Meta_analysis == g))

row_positions <- list()
diamond_rows  <- list()
header_y      <- list()
hetero_y      <- list()

cur_row <- 1   # build bottom-up

for (g in rev(subgroups)) {
  k_g <- rows_per_group[g]

  hetero_y[[g]]      <- cur_row + 0.7
  cur_row            <- cur_row + hetero_rows
  diamond_rows[[g]]  <- cur_row + 0.5
  cur_row            <- cur_row + 1 + 1          # +1 extra row = gap above diamond
  row_positions[[g]] <- cur_row:(cur_row + k_g - 1)
  cur_row            <- cur_row + k_g
  header_y[[g]]      <- cur_row + 0.3
  cur_row            <- cur_row + 1 + spacer
}

total_rows <- cur_row + 1

cat("\nRow layout summary:\n")
for (g in subgroups) {
  cat(" ", g, "— study rows:", min(row_positions[[g]]),
      "to", max(row_positions[[g]]),
      "| diamond at:", diamond_rows[[g]],
      "| header at:", header_y[[g]], "\n")
}

# ---- Sort each subgroup alphabetically by Citation ----
# Ensures authors always appear A→Z within each subgroup block.
escalc_dat <- escalc_dat %>%
  group_by(Meta_analysis) %>%
  arrange(ascii_key(Citation), Citation, .by_group = TRUE) %>%
  ungroup()

# ---- Assign rows and labels to every observation in escalc_dat ----
# This is the key step: we pre-assign a row number and label to each
# row of escalc_dat so a SINGLE forest() call places everything correctly.
escalc_dat$plot_row  <- NA_real_
escalc_dat$plot_slab <- NA_character_

for (g in subgroups) {
  idx <- which(escalc_dat$Meta_analysis == g)
  # assign rows top-to-bottom within the group's block
  # row_positions[[g]] is already ordered bottom-to-top, so reverse it
  escalc_dat$plot_row[idx]  <- rev(row_positions[[g]])
  escalc_dat$plot_slab[idx] <- escalc_dat$Citation[idx]
}

# ---- Layout ----
pt_cex <- if (total_rows > 80) 1.00 else
          if (total_rows > 50) 1.08 else
          if (total_rows > 30) 1.15 else 1.22

# No upper cap on height — let the page grow to fit all rows without overlap.
# Increase the multiplier here if rows still feel tight.
pdf_h  <- max(9, total_rows * 0.22 + 3.5)
# pdf_w is calculated AFTER xlim_f is known (see below)

# ---- Back-transform yi/vi to r-space so the plot axis is LINEAR in r ----
# Working in z-space and using atransf=transf.ztor causes unequal tick spacing
# because tanh is nonlinear. Instead, we convert the effect sizes and their
# CIs to r before passing them to forest(), then plot on a plain linear r axis.

# Back-transform point estimates and CI bounds to r
escalc_dat$yi_r  <- tanh(escalc_dat$yi)

# Delta-method variance in r-space: var(r) ≈ (1 - r²)^2 * var(z)
escalc_dat$vi_r  <- (1 - escalc_dat$yi_r^2)^2 * escalc_dat$vi

# CI bounds in r-space (used for axis range)
r_lb_vals  <- tanh(escalc_dat$yi - 1.96 * sqrt(escalc_dat$vi))
r_ub_vals  <- tanh(escalc_dat$yi + 1.96 * sqrt(escalc_dat$vi))
r_data_min <- min(r_lb_vals, na.rm = TRUE)
r_data_max <- max(r_ub_vals, na.rm = TRUE)

# Choose tick step and build evenly-spaced r-value ticks
at_step    <- if ((r_data_max - r_data_min) > 1.2) 0.25 else
              if ((r_data_max - r_data_min) > 0.6)  0.20 else 0.10
r_min_tick <- floor(r_data_min   / at_step) * at_step
r_max_tick <- ceiling(r_data_max / at_step) * at_step
r_ticks    <- seq(r_min_tick, r_max_tick, by = at_step)
# Ticks are already in r-space — no atanh conversion needed
at_vals    <- r_ticks

# ---- Coordinate layout (all in r units now) ----
# IMPORTANT: xlim_left must be far enough left to give the study labels (drawn
# at xlim_left) their own clear space, with the bars starting to the RIGHT of
# that. label_gap is the reserved width (in r units) for the study name column;
# bars begin at bar_left = r_data_min - plot_pad, which must sit >= xlim_left +
# label_gap so labels and bars never overlap.
plot_pad   <- 0.05                        # padding beyond outermost CI bar
label_gap  <- 1.40                        # space reserved for study name column
bar_left   <- r_data_min - plot_pad       # leftmost extent of any CI bar
xlim_left  <- bar_left - label_gap        # labels live in [xlim_left, bar_left)

r_plot_max <- r_data_max + plot_pad
bar_right  <- r_plot_max
xpos_ci    <- bar_right  + 0.04   # CI column — directly after bars
xpos_ante  <- xpos_ci    + 0.95   # Antecedent (risk_factor_display)
xpos_tlag  <- xpos_ante  + 1.40   # TL — wide gap so Antecedent text has room
xpos_wb    <- xpos_tlag  + 0.20   # WB
xpos_sn    <- xpos_wb    + 0.20   # SN
xlim_right <- xpos_sn    + 0.50   # generous right margin so SN text is never clipped
xlim_f     <- c(xlim_left, xlim_right)

# pdf_w: fixed wide enough to comfortably show all columns.
# The page is intentionally wider than strictly necessary so the
# right-hand annotation columns (CI, TL, WB, SN) always have room.
pdf_w  <- 20

# ---- Open PDF and draw ----
# showtext lets us embed Calibri (or any system font) without needing Cairo.
# Works on all platforms including Apple Silicon Macs.
tryCatch({

if (!requireNamespace("showtext", quietly = TRUE)) install.packages("showtext")
if (!requireNamespace("sysfonts",  quietly = TRUE)) install.packages("sysfonts")
library(showtext)
library(sysfonts)

# Try to register Calibri from common install locations.
calibri_loaded <- tryCatch({
  font_add("Calibri",
           regular = "/Library/Fonts/Microsoft/Calibri.ttf",
           bold    = "/Library/Fonts/Microsoft/Calibri Bold.ttf",
           italic  = "/Library/Fonts/Microsoft/Calibri Italic.ttf")
  TRUE
}, error = function(e) {
  tryCatch({
    font_add("Calibri",
             regular = "Calibri.ttf",
             bold    = "Calibrib.ttf",
             italic  = "Calibrii.ttf")
    TRUE
  }, error = function(e2) {
    message("Calibri not found — using default sans-serif font.")
    FALSE
  })
})

# Robust lookup for Sample_notes column
sn_col <- names(escalc_dat)[tolower(gsub("[ _.]", "", names(escalc_dat))) == "samplenotes"]
if (length(sn_col) == 0) {
  warning("Could not find Sample_notes column — SN column will be blank. ",
          "Column names available: ", paste(names(escalc_dat), collapse = ", "))
  escalc_dat$.sn_display <- ""
} else {
  cat("  Using Sample_notes column:", sn_col[1], "\n")
  escalc_dat$.sn_display <- as.character(escalc_dat[[ sn_col[1] ]])
}

# Antecedent column for the forest plot (from risk_factor_display; falls back
# to risk_factor if that column is absent)
.adc <- names(escalc_dat)[tolower(gsub("[ _.]", "", names(escalc_dat))) == "riskfactordisplay"]
if (length(.adc) == 0)
  .adc <- names(escalc_dat)[tolower(gsub("[ _.]", "", names(escalc_dat))) == "riskfactor"]
escalc_dat$.ante_display <- if (length(.adc))
  as.character(escalc_dat[[ .adc[1] ]]) else ""

# ---- Drawing function — called once for PDF, once for PNG ----
# This avoids any dependency on pdftoppm or magick for PNG export.
draw_forest_plot <- function() {
  if (calibri_loaded) par(family = "Calibri")
  op <- par(mar = c(3.5, 0.2, 1.0, 0.2), xpd = NA)

  forest(x         = escalc_dat$yi_r,
         vi        = escalc_dat$vi_r,
         ylim      = c(0, total_rows + 2),
         xlim      = xlim_f,
         at        = at_vals,
         xlab      = "Correlation (r)",
         rows      = escalc_dat$plot_row,
         slab      = escalc_dat$plot_slab,
         ilab      = data.frame(escalc_dat$.ante_display,
                                escalc_dat$time_lag,
                                as.character(escalc_dat$WB_measure),
                                escalc_dat$.sn_display),
         ilab.xpos = c(xpos_ante, xpos_tlag, xpos_wb, xpos_sn),
         ilab.pos  = 4,
         psize     = 1,
         cex       = pt_cex,
         header    = FALSE,
         anno      = FALSE,
         efac      = c(1, 0),
         mlab      = "")

  ci_labels <- formatC(tanh(escalc_dat$yi), format = "f", digits = 2)
  ci_lb     <- formatC(tanh(escalc_dat$yi - 1.96 * sqrt(escalc_dat$vi)), format = "f", digits = 2)
  ci_ub     <- formatC(tanh(escalc_dat$yi + 1.96 * sqrt(escalc_dat$vi)), format = "f", digits = 2)
  ci_text   <- paste0(ci_labels, " [", ci_lb, ", ", ci_ub, "]")
  text(xpos_ci, escalc_dat$plot_row, ci_text, pos = 4, cex = pt_cex * 0.8)

  par(font = 2)
  text(xlim_f[1], total_rows + 1.3, "Study", pos = 4, cex = pt_cex)
  text(xpos_ci,   total_rows + 1.3, "CI",    pos = 4, cex = pt_cex)
  text(xpos_ante, total_rows + 1.3, "Antecedent", pos = 4, cex = pt_cex)
  text(xpos_tlag, total_rows + 1.3, "TL",    pos = 4, cex = pt_cex)
  text(xpos_wb,   total_rows + 1.3, "WB",    pos = 4, cex = pt_cex)
  text(xpos_sn,   total_rows + 1.3, "SN",    pos = 4, cex = pt_cex)
  par(font = 1)

  for (g in subgroups) {
    fit_g <- sub_fits[[g]]
    df_g  <- escalc_dat[escalc_dat$Meta_analysis == g, ]

    par(font = 2)
    text(xlim_f[1], header_y[[g]], g, pos = 4, cex = pt_cex)
    par(font = 1)

    if (!is.null(fit_g)) {
      fit_g_r_est  <- tanh(as.numeric(fit_g$b[1, 1]))
      fit_g_r_cilb <- tanh(fit_g$ci.lb)
      fit_g_r_ciub <- tanh(fit_g$ci.ub)
      addpoly(x        = fit_g_r_est,
              ci.lb    = fit_g_r_cilb,
              ci.ub    = fit_g_r_ciub,
              row      = diamond_rows[[g]],
              mlab     = "",
              cex      = pt_cex,
              efac     = 0.8,
              col      = col_diamond,
              border   = col_diamond,
              annotate = FALSE)

      d_est <- formatC(fit_g_r_est,  format = "f", digits = 2)
      d_lb  <- formatC(fit_g_r_cilb, format = "f", digits = 2)
      d_ub  <- formatC(fit_g_r_ciub, format = "f", digits = 2)
      text(xpos_ci, diamond_rows[[g]],
           paste0(d_est, " [", d_lb, ", ", d_ub, "]"),
           pos = 4, cex = pt_cex * 0.9)

      if (inherits(fit_g, "rma.uni")) {
        i2_str   <- paste0(round(fit_g$I2, 1), "%")
        tau2_str <- round(fit_g$tau2, 2)
        H2_str   <- round(fit_g$H2, 2)
        line1 <- bquote(paste(
          "Heterogeneity: ",
          tau^2, " = ", .(tau2_str), ",  ",
          I^2,   " = ", .(i2_str),   ",  ",
          H^2,   " = ", .(H2_str)
        ))
        line2 <- bquote(paste(
          "Test of ", theta[i], " = ", theta[j],
          ": Q(", .(fit_g$k - 1), ") = ", .(round(fit_g$QE, 2)),
          ", p = ", .(formatC(fit_g$QEp, digits = 2, format = "f"))
        ))
      } else {
        typical_vi <- median(df_g$vi, na.rm = TRUE)
        sigma2_tot <- sum(fit_g$sigma2)
        total_var  <- sigma2_tot + typical_vi
        i2_total   <- round(100 * sigma2_tot / total_var, 1)
        i2_between <- round(100 * fit_g$sigma2[1] / total_var, 1)
        i2_within  <- round(100 * fit_g$sigma2[2] / total_var, 1)
        line1 <- bquote(paste(
          "Heterogeneity: ",
          tau^2, " between = ", .(round(fit_g$sigma2[1], 2)),
          ", within = ", .(round(fit_g$sigma2[2], 2))
        ))
        line2 <- paste0("I\u00B2 total = ", i2_total,
                        "% (between = ", i2_between,
                        "%, within = ", i2_within, "%)")
      }

      text(xlim_f[1], hetero_y[[g]],        line1, pos = 4, cex = pt_cex * 0.88)
      text(xlim_f[1], hetero_y[[g]] - 0.85, line2, pos = 4, cex = pt_cex * 0.88)
    }
  }
  par(op)
}

pdf_out <- paste0("forest_plot_", dataset_tag, ".pdf")

# ---- Render to PDF (per-category forest only exported for these two) ----
if (current_category %in% c("Gender", "Psychological distress")) {
  showtext_auto()
  showtext_opts(dpi = 300)
  pdf(pdf_out, width = pdf_w, height = pdf_h, pointsize = 11)
  draw_forest_plot()
  dev.off()
  showtext_auto(FALSE)
  cat("\nGrouped forest plot saved to", pdf_out, "\n")
} else {
  cat("\n(forest_plot_", dataset_tag, ".pdf not exported — not in keep-list)\n", sep = "")
}

# ---- PNG export disabled (PDFs only) ----


}, error = function(e) {
  # Always close the PDF device to avoid a corrupt/empty file
  if (dev.cur() > 1) dev.off()
  tryCatch(showtext_auto(FALSE), error = function(e2) NULL)
  cat("\n*** Forest plot error — details below. The rest of the script will continue. ***\n")
  cat("Error message:", conditionMessage(e), "\n")
})

# Reset graphics state so subsequent plots are not affected by any
# par() settings left over from inside the forest tryCatch block
graphics.off()
par(mar = c(5.1, 4.1, 4.1, 2.1))   # restore R defaults


### 6) Funnel plot + Egger test --------------------------------------

sei     <- sqrt(escalc_dat$vi)
k_f     <- if (inherits(res_main, "rma.mv")) nrow(escalc_dat) else res_main$k
pt_cex2 <- if (k_f > 100) 0.45 else if (k_f > 70) 0.55 else if (k_f > 40) 0.65 else 0.75

# ---- 6A: Base metafor funnel ----
tryCatch({
  png(paste0("funnel_plot_", dataset_tag, ".png"), width = 6.2, height = 5.4, units = "in", res = 300, pointsize = 9)
  op2 <- par(mar = c(3.6, 4.0, 2.6, 1.0), mgp = c(2.1, 0.60, 0), tcl = -0.25, cex = 0.9)
  funnel(res_main,
         yaxis   = "sei",
         xlab    = "Fisher's z (r-to-z)",
         ylab    = "Standard error",
         pch     = 16,
         cex     = pt_cex2,
         refline = 0,
         level   = c(85, 90, 95, 99),
         shade   = c("gray98","gray94","gray90","gray86"),
         back    = "white")
  title(main = current_category, adj = 0, font.main = 2, cex.main = 1.1, line = 1.0)
  abline(v = pooled_z, lty = 3, col = "gray40")
  par(op2)
  dev.off()
  cat("Funnel plot saved to", paste0("funnel_plot_", dataset_tag, ".pdf"), "\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("*** Funnel plot error:", conditionMessage(e), "\n")
})

# ---- 6B: Multilevel (RVE) Egger-type asymmetry test ----
egger_all[[current_category]] <- ml_egger(escalc_dat, label = current_category)

# ---- 6C: ggplot2 funnel (journal style, x-axis in r) ----
sei_seq <- seq(min(sei, na.rm = TRUE), max(sei, na.rm = TRUE), length.out = 300)

mk_contour_df <- function(p_two) {
  zc <- qnorm(1 - p_two / 2)
  data.frame(sei = sei_seq, r_lo = tanh(-zc * sei_seq), r_hi = tanh(zc * sei_seq), p = p_two)
}
df90 <- mk_contour_df(0.10)
df95 <- mk_contour_df(0.05)
df99 <- mk_contour_df(0.01)

df_pts   <- data.frame(r = tanh(escalc_dat$yi), sei = sei)
xlim_r   <- range(df_pts$r, na.rm = TRUE)
xlim_r   <- c(max(-0.999, xlim_r[1] - 0.03), min(0.999, xlim_r[2] + 0.03))
pt_size  <- if (k_f > 100) 0.7 else if (k_f > 70) 1.0 else 1.3

p_funnel <- ggplot(df_pts, aes(x = r, y = sei)) +
  geom_line(data = df90, aes(x = r_lo, y = sei), colour = "grey70", linewidth = 0.3, linetype = 2, inherit.aes = FALSE) +
  geom_line(data = df90, aes(x = r_hi, y = sei), colour = "grey70", linewidth = 0.3, linetype = 2, inherit.aes = FALSE) +
  geom_line(data = df95, aes(x = r_lo, y = sei), colour = "grey60", linewidth = 0.35, linetype = 3, inherit.aes = FALSE) +
  geom_line(data = df95, aes(x = r_hi, y = sei), colour = "grey60", linewidth = 0.35, linetype = 3, inherit.aes = FALSE) +
  geom_line(data = df99, aes(x = r_lo, y = sei), colour = "grey50", linewidth = 0.4, linetype = 4, inherit.aes = FALSE) +
  geom_line(data = df99, aes(x = r_hi, y = sei), colour = "grey50", linewidth = 0.4, linetype = 4, inherit.aes = FALSE) +
  geom_point(size = pt_size, shape = 16, alpha = 0.75) +
  geom_vline(xintercept = 0,        colour = "grey55", linewidth = 0.4) +
  geom_vline(xintercept = pooled_r, colour = "grey40", linewidth = 0.4, linetype = "dashed") +
  scale_y_reverse(expand = expansion(mult = c(0.02, 0.04))) +
  coord_cartesian(xlim = xlim_r) +
  labs(x = "Correlation (r)", y = "Standard error") +
  theme_minimal(base_size = 9) +
  theme(plot.margin = margin(6,6,6,6),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank())

# DISABLED (keep simple funnel only):
# ggsave(paste0("funnel_plot_ggplot_", dataset_tag, ".pdf"), p_funnel, width = 6.2, height = 5.4, units = "in")
# cat("ggplot funnel saved to", paste0("funnel_plot_ggplot_", dataset_tag, ".pdf"), "\n")


# ---- 6D: Egger test (safe fallback for rma.mv) ----
if (inherits(res_main, "rma.uni")) {
  cat("\nEgger test:\n")
  print(regtest(res_main, model = "rma"))
} else {
  cat("\nEgger test (approximated via sample/cluster-level aggregation for rma.mv):\n")
  cat("  Interpret with caution.\n")

  agg_uni <- escalc_dat %>%
    group_by(cluster_id) %>%
    summarise(
      yi = sum(yi / vi, na.rm = TRUE) / sum(1 / vi, na.rm = TRUE),
      vi = 1 / sum(1 / vi, na.rm = TRUE),
      .groups = "drop"
    )
  res_agg <- rma(yi, vi, data = agg_uni, method = "REML")
  print(regtest(res_agg, model = "rma"))
}


### 7) Leave-one-out sensitivity analysis ----------------------------
#
#  Iteratively removes one study at a time and refits the model.
#  Saves a table + a caterpillar plot of pooled r at each leave-out.

cat("\n=== Leave-one-out sensitivity analysis ===\n")

study_ids <- unique(escalc_dat$study_ID)
# Label each omitted study by its Citation rather than the numeric study_ID.
loo_labels <- as.character(escalc_dat$Citation[match(study_ids, escalc_dat$study_ID)])
if (any(duplicated(loo_labels)))            # guard if one Citation spans >1 study_ID
  loo_labels <- make.unique(loo_labels, sep = " #")
loo_rows  <- vector("list", length(study_ids))

for (i in seq_along(study_ids)) {
  tmp <- escalc_dat[escalc_dat$study_ID != study_ids[i], ]

  if (multi_es) {
    fit <- tryCatch(
      rma.mv(yi, vi, random = ~ 1 | cluster_id/effect_id,
             data = tmp, method = "REML"),
      error = function(e) NULL
    )
  } else {
    fit <- tryCatch(
      rma(yi, vi, data = tmp, method = "REML"),
      error = function(e) NULL
    )
  }

  if (is.null(fit)) {
    loo_rows[[i]] <- data.frame(
      study_omitted = loo_labels[i],
      r_pooled = NA, r_lb = NA, r_ub = NA
    )
  } else {
    z_b  <- as.numeric(fit$b[1, 1])
    loo_rows[[i]] <- data.frame(
      study_omitted = loo_labels[i],
      r_pooled      = tanh(z_b),
      r_lb          = tanh(fit$ci.lb),
      r_ub          = tanh(fit$ci.ub)
    )
  }
}

loo_df <- do.call(rbind, loo_rows)
cat("\nLeave-one-out pooled r (sorted):\n")
print(loo_df[order(loo_df$r_pooled), ])

# DISABLED (keep LOO PDF only):
# write_csv(loo_df, paste0("sensitivity_loo_", dataset_tag, ".csv"))
# cat("Saved leave-one-out results to", paste0("sensitivity_loo_", dataset_tag, ".csv"), "\n")

# Caterpillar plot
p_loo <- ggplot(loo_df, aes(x = r_pooled,
                             y = reorder(study_omitted, r_pooled))) +
  geom_errorbarh(aes(xmin = r_lb, xmax = r_ub),
                 height = 0.3, colour = "grey50") +
  geom_point(shape = 21, fill = "#4292C6", size = 2) +
  geom_vline(xintercept = pooled_r,
             linetype = "dashed", colour = "grey30") +
  labs(x = "Pooled r (omitting study)",
       y = "Study omitted",
       title = current_category,
       subtitle = paste0("Dashed line = full-model r = ", round(pooled_r, 3))) +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))

ggsave(paste0("sensitivity_loo_", dataset_tag, ".png"), p_loo,
       width = 7, height = max(5, 0.25 * length(study_ids) + 2),
       units = "in", dpi = 300)
cat("Leave-one-out plot saved to", paste0("sensitivity_loo_", dataset_tag, ".pdf"), "\n")


### 8) Subgroup / Meta-regression on time lag ------------------------

cat("\nTime lag NAs:", sum(is.na(escalc_dat$time_lag)), "\n")

escalc_dat$time_lag_cat <- factor(
  escalc_dat$time_lag,
  levels = sort(unique(escalc_dat$time_lag)),
  labels = paste0("T", sort(unique(escalc_dat$time_lag)))
)

# ---- 8A: Categorical subgroup ----
cat("\n--- Subgroup (categorical) by time_lag_cat ---\n")
if (multi_es) {
  res_sub_cat <- rma.mv(yi, vi,
                        mods   = ~ time_lag_cat - 1,
                        random = ~ 1 | cluster_id/effect_id,
                        data   = escalc_dat, method = "REML")
  print(res_sub_cat)
  cat("\nCR2 robust tests:\n")
  print(clubSandwich::coef_test(res_sub_cat, vcov = "CR2",
                                cluster = escalc_dat$cluster_id))
} else {
  res_sub_cat <- rma(yi, vi, mods = ~ time_lag_cat - 1,
                     data = escalc_dat, method = "REML")
  print(res_sub_cat)
}

# Back-transform subgroup estimates to r
est_cat <- data.frame(
  time_lag_cat = gsub("^time_lag_cat", "", rownames(res_sub_cat$b)),
  z_est        = as.numeric(res_sub_cat$b),
  z_ci_lb      = res_sub_cat$ci.lb,
  z_ci_ub      = res_sub_cat$ci.ub
) %>%
  mutate(r_est   = tanh(z_est),
         r_ci_lb = tanh(z_ci_lb),
         r_ci_ub = tanh(z_ci_ub))
cat("\nBack-transformed subgroup means (r):\n")
print(est_cat)


# ---- 8B: Continuous meta-regression with time_lag (per subgroup) ----
cat("\n--- Meta-regression (continuous): time_lag, per Meta_analysis subgroup ---\n")

run_metareg <- function(df, moderator, label, colour, file_tag) {
  # Fit model
  dup <- any(duplicated(df$cluster_id))
  mod_formula <- as.formula(paste0("~ ", moderator))
  fit <- tryCatch({
    if (dup)
      rma.mv(yi, vi, mods = mod_formula,
             random = ~ 1 | cluster_id/effect_id,
             data = df, method = "REML")
    else
      rma(yi, vi, mods = mod_formula, data = df, method = "REML", test = "knha")
  }, error = function(e) { cat("  Model failed:", e$message, "\n"); NULL })

  if (is.null(fit)) return(invisible(NULL))

  cat("\n  Model for", label, ":\n")
  print(fit)

  # Extract t-stat and p-value for the moderator coefficient.
  # For rma.mv use CR2 robust SE; for rma.uni use the model's own zval/pval.
  t_stat <- NA_real_; p_val <- NA_real_; mod_df <- NA_real_
  n_cl   <- length(unique(df$cluster_id))

  if (inherits(fit, "rma.mv")) {
    cr2 <- tryCatch(
      clubSandwich::coef_test(fit, vcov = "CR2", cluster = df$cluster_id),
      error = function(e) { cat("  CR2 failed:", e$message, "\n"); NULL }
    )
    cat("  CR2 robust test:\n")
    print(cr2)
    if (!is.null(cr2) && nrow(cr2) >= 2) {
      t_stat <- as.numeric(cr2$tstat[2])
      p_val  <- as.numeric(cr2$p_Satt[2])
      mod_df <- as.numeric(cr2$df_Satt[2])
    }
  } else {
    if (length(fit$zval) >= 2) {
      t_stat <- fit$zval[2]
      p_val  <- fit$pval[2]
      mod_df <- if (!is.null(fit$dfs)) as.numeric(fit$dfs[2]) else (fit$k - fit$p)
    }
  }

  # --- df diagnostic: expose what each df is built from ---
  used_k <- if (!is.null(fit$k)) fit$k else NA_integer_   # effects the model actually used
  cl_len <- length(df$cluster_id)                          # length of cluster vector sent to CR2
  cat(sprintf("  [df-check: %s] rows in df = %d | effects used by model = %s | clusters(full df) = %d | cluster-vector length = %d | reported df = %s%s\n",
              label, nrow(df), ifelse(is.na(used_k), "?", used_k), n_cl, cl_len,
              ifelse(is.na(mod_df), "NA", round(mod_df, 2)),
              ifelse(!is.na(used_k) && used_k != cl_len,
                     "  <-- MISMATCH: model dropped rows but full cluster vector was passed to CR2 (df suspect)",
                     "")))
  
  # df-floor guard: a meta-regression slope needs enough INDEPENDENT samples to
  # be interpretable. Below the floor we report nothing rather than an artefact
  # (a huge t on ~1 residual df is fit noise, not evidence).
  min_clusters_cont <- 3
  if (n_cl < min_clusters_cont || is.na(mod_df) || mod_df < 1) {
    cat(sprintf("  [%s] moderator suppressed: %d clusters, df = %s (below floor)\n",
                label, n_cl, ifelse(is.na(mod_df), "NA", round(mod_df, 1))))
    t_stat <- NA_real_; p_val <- NA_real_; mod_df <- NA_real_
  }

  # Prediction grid
  mod_range <- range(df[[moderator]], na.rm = TRUE)
  new_df    <- setNames(data.frame(seq(mod_range[1], mod_range[2], length.out = 100)),
                        moderator)
  preds <- tryCatch(
    predict(fit, newmods = as.matrix(new_df[[moderator]])),
    error = function(e) NULL
  )
  if (is.null(preds)) return(invisible(NULL))

  pred_df <- data.frame(
    x    = new_df[[moderator]],
    r    = tanh(preds$pred),
    r_lb = tanh(preds$ci.lb),
    r_ub = tanh(preds$ci.ub)
  )
  obs_df <- data.frame(x = df[[moderator]], r = tanh(df$yi))

  p <- ggplot(pred_df, aes(x = x, y = r)) +
    geom_ribbon(aes(ymin = r_lb, ymax = r_ub), alpha = 0.15) +
    geom_line(linewidth = 1.0, colour = colour) +
    geom_point(data = obs_df, aes(x = x, y = r),
               alpha = 0.55, size = 1.5, inherit.aes = FALSE) +
    labs(x = label, y = "Predicted r",
         title = paste0("Meta-regression: r ~ ", label),
         subtitle = file_tag) +
    theme_minimal(base_size = 10)

  # Return both the plot and the stats as a named list
  return(list(plot = p, subgroup = file_tag, moderator = label,
              t_stat = t_stat, p_val = p_val, df = mod_df, n_clusters = n_cl))
}

# Helper to save a multi-panel plot (one panel per subgroup) and collect stats
run_subgroup_metareg <- function(moderator, x_label, colour, file_prefix) {
  plots <- list()
  for (g in subgroups) {
    df_g <- escalc_dat[escalc_dat$Meta_analysis == g, ]
    n_obs <- sum(!is.na(df_g[[moderator]]))
    if (n_obs < 3) {
      cat("  Skipping", g, "— fewer than 3 non-NA values for", moderator, "\n")
      next
    }
    cat("\nSubgroup:", g, "\n")
    result <- run_metareg(df_g, moderator, x_label, colour, file_tag = g)
    if (!is.null(result)) {
      plots[[g]] <- result$plot
      # Accumulate stats into the global summary list
      moderation_results_all[[length(moderation_results_all) + 1]] <<- data.frame(
        Moderator = x_label,
        Subgroup  = g,
        t_stat    = result$t_stat,
        p_val     = result$p_val,
        df         = result$df,
        n_clusters = result$n_clusters,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(plots) == 0) return(invisible(NULL))

  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
  library(patchwork)
  combined <- Reduce(`/`, plots)
  out_file <- paste0(file_prefix, "_by_subgroup_", dataset_tag, ".pdf")
  # DISABLED (moderation plots not in keep-list):
  # ggsave(out_file, combined,
  #        width = 7, height = 4 * length(plots), units = "in", limitsize = FALSE)
  # cat("\nSaved subgroup meta-regression plot to", out_file, "\n")
}

cat("\n=== Time lag meta-regressions by subgroup ===\n")
# Initialise the global list that collects stats from all three continuous moderators
moderation_results_all <- list()
run_subgroup_metareg("time_lag", "Time lag (months)", "#238B45", "metareg_timelag")


# ---- 8C: Quadratic time-lag (overall) ----
escalc_dat$time_lag_c  <- escalc_dat$time_lag - mean(escalc_dat$time_lag, na.rm = TRUE)
escalc_dat$time_lag_c2 <- escalc_dat$time_lag_c^2

cat("\n--- Meta-regression with quadratic time_lag (overall) ---\n")
if (multi_es) {
  res_quad <- rma.mv(yi, vi, mods = ~ time_lag_c + time_lag_c2,
                     random = ~ 1 | cluster_id/effect_id,
                     data = escalc_dat, method = "REML")
  print(res_quad)
  cat("\nCR2 robust tests (quadratic):\n")
  print(clubSandwich::coef_test(res_quad, vcov = "CR2",
                                cluster = escalc_dat$cluster_id))
} else {
  res_quad <- rma(yi, vi, mods = ~ time_lag_c + time_lag_c2,
                  data = escalc_dat, method = "REML")
  print(res_quad)
}


### 8D) Prepare moderator variables: % female & mean age -------------

cat("\nModerator prep: mean age and % female\n")
print(names(escalc_dat))

# *** EDIT this if your column has a different name ***
age_col <- "M_age"

if (!age_col %in% names(escalc_dat)) {
  stop(paste("Column", age_col, "not found. Please update age_col."))
}

escalc_dat <- escalc_dat %>%
  mutate(
    mean_age  = suppressWarnings(as.numeric(.data[[age_col]])),
    female_pc = suppressWarnings(as.numeric(female))
  )

# Detect and fix % female recorded as proportion (0-1) instead of percentage
if (max(escalc_dat$female_pc, na.rm = TRUE) <= 1) {
  cat("  NOTE: female_pc appears to be a proportion (max =",
      max(escalc_dat$female_pc, na.rm = TRUE),
      ") — multiplying by 100 to convert to percentage.\n")
  escalc_dat$female_pc <- escalc_dat$female_pc * 100
}

cat("  Missing mean_age  :", sum(is.na(escalc_dat$mean_age)),  "\n")
cat("  Missing female_pc :", sum(is.na(escalc_dat$female_pc)), "\n")

# ---- Risk-of-bias score (continuous moderator) ----
if (length(rob_col)) {
  escalc_dat$rob_score <- suppressWarnings(as.numeric(as.character(escalc_dat[[rob_col[1]]])))
  cat("  Using risk-of-bias column:", rob_col[1],
      "| missing:", sum(is.na(escalc_dat$rob_score)), "\n")
} else {
  escalc_dat$rob_score <- NA_real_
}


### 8E) Moderation: % female — per subgroup --------------------------

cat("\n=== % Female meta-regressions by subgroup ===\n")
run_subgroup_metareg("female_pc", "% Female", "#8856A7", "metareg_female")


### 8F) Moderation: mean age — per subgroup --------------------------

cat("\n=== Mean age meta-regressions by subgroup ===\n")
run_subgroup_metareg("mean_age", "Mean age (years)", "#E6550D", "metareg_age")


### 8I) Moderation: risk of bias — per subgroup ----------------------

cat("\n=== Risk-of-bias meta-regressions by subgroup ===\n")
if (any(!is.na(escalc_dat$rob_score))) {
  run_subgroup_metareg("rob_score", "Risk of bias", "#756BB1", "metareg_riskofbias")
} else {
  cat("  Skipped — no numeric risk-of-bias values available for this category.\n")
}


### 8G) Optional interaction: % female × mean age --------------------
#
#  Uncomment if you want to test whether the two moderators interact.
#  Requires sufficient studies with data on both variables.
#
# cat("\n--- Interaction: % female × mean age ---\n")
# escalc_dat <- escalc_dat %>%
#   mutate(
#     female_c   = female_pc - mean(female_pc, na.rm = TRUE),
#     age_c      = mean_age  - mean(mean_age,  na.rm = TRUE),
#     fem_x_age  = female_c * age_c
#   )
# if (multi_es) {
#   res_inter <- rma.mv(yi, vi,
#                       mods   = ~ female_c + age_c + fem_x_age,
#                       random = ~ 1 | cluster_id/effect_id,
#                       data   = escalc_dat, method = "REML")
# } else {
#   res_inter <- rma(yi, vi, mods = ~ female_c + age_c + fem_x_age,
#                    data = escalc_dat, method = "REML")
# }
# print(res_inter)
# if (multi_es) {
#   print(clubSandwich::coef_test(res_inter, vcov = "CR2",
#                                 cluster = escalc_dat$cluster_id))
# }


### 8H) Moderation: WB_measure (categorical) — per Meta_analysis subgroup ----
#
#  Runs a separate categorical moderation model for each Meta_analysis subgroup,
#  then stacks the panels into a single PDF — matching the pattern of 8E/8F.

cat("\n=== WB_measure moderation by Meta_analysis subgroup ===\n")

# Clean WB_measure once for the whole dataset
escalc_dat$WB_measure <- trimws(as.character(escalc_dat$WB_measure))
escalc_dat$WB_measure <- ifelse(
  escalc_dat$WB_measure %in% c("", "NA", "na", "N/A", "n/a"),
  NA_character_,
  escalc_dat$WB_measure
)

# ============================================================
# Helper: Run WB_measure moderation for one subgroup
# ============================================================


run_wb_moderation <- function(df, subgroup_label) {
  
  # ============================================================
  # 1. Clean & prepare moderator
  # ============================================================
  
  df <- df[!is.na(df$WB_measure), ]
  
  # Drop unused levels and enforce factor structure
  df$WB_measure_f <- droplevels(factor(df$WB_measure))
  
  # Must have ≥2 categories in this subgroup
  if (nlevels(df$WB_measure_f) < 2) {
    cat("  Skipping", subgroup_label, "— fewer than 2 WB_measure categories.\n")
    return(invisible(NULL))
  }
  
  # Set SI as reference category (ALWAYS)
  if ("SL" %in% levels(df$WB_measure_f)) {
    df$WB_measure_f <- relevel(df$WB_measure_f, ref = "SL")
  }
  
  dup <- any(duplicated(df$cluster_id))
  
  # ============================================================
  # 2. Fit TWO models:
  #
  #   res_wb_display  — no-intercept (~WB_measure_f - 1)
  #     Each coefficient is the mean Fisher's z for that category,
  #     so we can back-transform and plot each group's absolute r.
  #
  #   res_wb_contrast — intercept (~WB_measure_f, SI as reference)
  #     The single contrast coefficient is Q vs SI.
  #     anova() on this model produces Q(1), i.e. the omnibus test
  #     of whether the two categories differ — which is what we want.
  #     Using the no-intercept model for anova() gives Q(2) because
  #     it tests whether *both* group means equal zero, not whether
  #     they differ from each other.
  # ============================================================
  
  fit_model <- function(formula_rhs) {
    tryCatch({
      if (dup) {
        rma.mv(
          yi, vi,
          mods   = formula_rhs,
          random = ~ 1 | cluster_id/effect_id,
          data   = df,
          method = "REML"
        )
      } else {
        rma(
          yi, vi,
          mods   = formula_rhs,
          data   = df,
          method = "REML",
          test   = "knha"
        )
      }
    }, error = function(e) {
      cat("  Model failed:", e$message, "\n")
      NULL
    })
  }
  
  # Display model: one mean estimate per category (no intercept)
  res_wb_display  <- fit_model(~ WB_measure_f - 1)
  # Contrast model: Q vs SI (intercept = SI; one free parameter)
  res_wb_contrast <- fit_model(~ WB_measure_f)
  
  # Both models must succeed to proceed
  if (is.null(res_wb_display) || is.null(res_wb_contrast)) return(invisible(NULL))
  
  # Use the display model for printing and coefficient extraction
  res_wb <- res_wb_display
  
  cat("\n--- Subgroup:", subgroup_label, "---\n")
  cat("  Display model (no intercept — one estimate per category):\n")
  print(res_wb_display)
  cat("  Contrast model (intercept = SI; coefficient = Q vs SI):\n")
  print(res_wb_contrast)
  
  
  # ============================================================
  # 3. CR2 robust tests (on the contrast model — one row = Q vs SI)
  # ============================================================
  wb_cr2 <- NULL
  if (inherits(res_wb_contrast, "rma.mv")) {
    cat("  CR2 robust tests (contrast model):\n")
    wb_cr2 <- tryCatch(
      clubSandwich::coef_test(res_wb_contrast, vcov = "CR2",
                              cluster = df$cluster_id),
      error = function(e) { cat("  CR2 failed:", e$message, "\n"); NULL }
    )
    print(wb_cr2)
  }
  
  # ============================================================
  # 4. Omnibus test — must use the CONTRAST model so df = 1.
  #
  #    For rma.mv: anova() by default tests ALL coefficients
  #    (intercept + contrast), giving Q(2). We must pass btt = 2
  #    to restrict the test to coefficient 2 only (the Q vs SI
  #    contrast), which gives Q(1).
  #
  #    For rma.uni: anova() already excludes the intercept and
  #    tests only the moderator coefficient(s), so Q(1) comes
  #    out automatically — no btt needed.
  # ============================================================
  omnibus <- tryCatch(
    if (inherits(res_wb_contrast, "rma.mv"))
      anova(res_wb_contrast, btt = 2)   # btt=2 → test coef 2 only (Q vs SI)
    else
      anova(res_wb_contrast),            # rma.uni already excludes intercept
    error = function(e) NULL
  )
  cat("  Omnibus test (SI vs Q, df should = 1):\n")
  print(omnibus)
  
  # Omnibus moderator test (SL vs Q) from the multilevel contrast model, with a
  # per-side cluster floor: we need enough INDEPENDENT samples in EACH measure
  # group, otherwise the contrast is confounded with study and we report nothing.
  # (The CR2 robust t printed above is kept for reference in the console log.)
  omni_Q <- NA_real_; omni_df <- NA_real_; omni_p <- NA_real_; omni_label <- ""
  per_side <- colSums(table(df$cluster_id, df$WB_measure) > 0)  # distinct clusters/measure
  min_clusters_side <- 2
  if (length(per_side) >= 2 && all(per_side >= min_clusters_side) && !is.null(omnibus)) {
    omni_Q  <- as.numeric(omnibus$QM)
    omni_df <- as.numeric(omnibus$QMdf[1])
    omni_p  <- as.numeric(omnibus$QMp)
    if (!is.na(omni_Q))
      omni_label <- paste0("SL vs Q omnibus: Q(", round(omni_df), ") = ",
                           formatC(omni_Q, 2, format = "f"), ", ",
                           ifelse(omni_p < .001, "p < .001",
                                  paste0("p = ", formatC(omni_p, 3, format = "f"))))
  } else {
    cat(sprintf("  SL vs Q suppressed: clusters per measure = %s (need >= %d each)\n",
                paste(per_side, collapse = "/"), min_clusters_side))
  }
  
  
  # ============================================================
  # 5. Extract coefficients robustly
  #
  #  rma.mv and rma.uni store coefficient tables differently.
  #  coef(summary()) works for both model types.
  # ============================================================

  coef_tab <- tryCatch(
    as.data.frame(coef(summary(res_wb))),
    error = function(e) {
      cat("  Coefficient extraction failed:", e$message, "\n")
      NULL
    }
  )

  if (is.null(coef_tab) || nrow(coef_tab) == 0) {
    cat("  Moderator dropped — skipping.\n")
    return(NULL)
  }

  coef_tab$WB_measure <- rownames(coef_tab)
  rownames(coef_tab) <- NULL

  # Keep only the category rows (WB_measure_fSI, WB_measure_fQ, etc.)
  coef_tab <- coef_tab[grep("^WB_measure_f", coef_tab$WB_measure), ]

  if (nrow(coef_tab) == 0) {
    cat("  Moderator dropped — skipping.\n")
    return(NULL)
  }
  
  # Add k per category
  coef_tab$k <- as.integer(table(df$WB_measure_f))[
    match(gsub("^WB_measure_f", "", coef_tab$WB_measure),
          levels(df$WB_measure_f))
  ]
  
  # Back-transform estimates to r
  coef_tab <- coef_tab %>%
    mutate(
      WB_measure = gsub("^WB_measure_f", "", WB_measure),
      r_est      = tanh(estimate),
      r_ci_lb    = tanh(ci.lb),
      r_ci_ub    = tanh(ci.ub),
      p_label    = ifelse(pval < .001, "p < .001",
                          paste0("p = ", formatC(pval, digits = 3, format = "f"))),
      point_label = paste0(
        "r = ", formatC(r_est, digits = 2, format = "f"),
        " [", formatC(r_ci_lb, digits = 2, format = "f"), ", ",
        formatC(r_ci_ub, digits = 2, format = "f"), "]   k = ", k,
        "   ", p_label
      )
    )
  
  
  # ============================================================
  # 6. Build plot panel
  # ============================================================
  
  p <- ggplot(coef_tab, aes(x = r_est, y = reorder(WB_measure, r_est))) +
    geom_errorbarh(aes(xmin = r_ci_lb, xmax = r_ci_ub),
                   height = 0.25, colour = "grey50") +
    geom_point(shape = 21, fill = "#2171B5", size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text(aes(x = r_ci_ub, label = point_label),
              hjust = -0.08, size = 2.8, colour = "grey20") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.8))) +
    coord_cartesian(clip = "off") +
    labs(x = "Pooled correlation (r)",
         y = NULL,
         title = subgroup_label,
         caption = omni_label) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 9),
      plot.caption = element_text(size = 8, hjust = 0,
                                  face = "italic", colour = "grey30"),
      plot.margin = margin(6, 6, 10, 6)
    )
  
  
  return(list(
    plot = p,
    data = coef_tab,
    subgroup = subgroup_label,
    omnibus = omni_label,
    omni_Q = omni_Q,
    omni_df = omni_df,
    omni_p = omni_p
  ))
}

# ============================================================
# Run moderation for each subgroup
# ============================================================

wb_panels   <- list()
wb_data_all <- list()
# Capture the omnibus (SI vs Q) for this category's combined-table row
wb_omni_Q <- NA_real_; wb_omni_df <- NA_real_; wb_omni_p <- NA_real_

for (g in subgroups) {
  df_g <- escalc_dat[escalc_dat$Meta_analysis == g, ]
  result <- run_wb_moderation(df_g, subgroup_label = g)
  
  if (!is.null(result)) {
    wb_panels[[g]] <- result$plot
    df_tagged <- result$data %>%
      mutate(Subgroup = g, Omnibus = result$omnibus)
    wb_data_all[[g]] <- df_tagged
    wb_omni_Q <- result$omni_Q; wb_omni_df <- result$omni_df; wb_omni_p <- result$omni_p
  }
}

# ============================================================
# Save combined plot
# ============================================================

if (length(wb_panels) > 0) {
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
  library(patchwork)
  
  combined_wb <- Reduce(`/`, wb_panels)
  out_file_wb <- paste0("moderation_WB_measure_by_subgroup_", dataset_tag, ".pdf")
  # DISABLED (moderation plots not in keep-list):
  # ggsave(out_file_wb, combined_wb,
  #        width = 9,
  #        height = max(4, 3.5 * length(wb_panels)),
  #        units = "in", limitsize = FALSE)
  # cat("\nSaved WB_measure subgroup moderation plot to", out_file_wb, "\n")
} else {
  cat("  No subgroups had sufficient WB_measure variation — no plot saved.\n")
}

# ============================================================
# Summary table (CSV + PDF)
# ============================================================

if (length(wb_data_all) > 0) {
  
  wb_summary <- do.call(rbind, wb_data_all) %>%
    mutate(
      Result = paste0(
        "r = ", formatC(r_est, digits = 2, format = "f"),
        ", 95% CI [", formatC(r_ci_lb, digits = 2, format = "f"),
        ", ", formatC(r_ci_ub, digits = 2, format = "f"),
        "], k = ", k
      )
    ) %>%
    select(Subgroup, `Measurement type` = WB_measure,
           Result, `p-value` = p_label, Omnibus)
  
  wb_csv_out <- paste0("WB_measure_summary_", dataset_tag, ".csv")
  # DISABLED (WB_measure summary not in keep-list):
  # write_csv(wb_summary, wb_csv_out)
  # cat("WB_measure summary CSV saved to", wb_csv_out, "\n")
  
  if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
  library(gridExtra)
  library(grid)
  
  if (FALSE) {   # DISABLED: WB_measure summary PDF not in keep-list
  wb_pdf_out <- paste0("WB_measure_summary_", dataset_tag, ".pdf")
  
  make_wb_section <- function(sg_name) {
    
    df_sec <- wb_summary[wb_summary$Subgroup == sg_name,
                         c("Measurement type", "Result", "p-value")]
    omni <- unique(wb_summary$Omnibus[wb_summary$Subgroup == sg_name])[1]
    
    tt <- ttheme_minimal(
      core = list(
        fg_params = list(hjust = 0, x = 0.04, fontsize = 9),
        bg_params = list(fill = c("white", "#f5f5f5"))
      ),
      colhead = list(
        fg_params = list(hjust = 0, x = 0.04, fontsize = 9, fontface = "bold"),
        bg_params = list(fill = "#2c2c5e", col = NA),
        fg_params2 = list(col = "white")
      )
    )
    
    tbl_grob <- tableGrob(df_sec, rows = NULL, theme = tt)
    title_grob <- textGrob(sg_name, gp = gpar(fontsize = 11, fontface = "bold"),
                           hjust = 0, x = 0)
    omni_grob <- textGrob(omni, gp = gpar(fontsize = 8, fontface = "italic",
                                          col = "grey40"),
                          hjust = 0, x = 0)
    
    arrangeGrob(title_grob, tbl_grob, omni_grob, ncol = 1,
                heights = unit(c(0.35, nrow(df_sec) * 0.32 + 0.4, 0.3), "inches"))
  }
  
  sg_names <- unique(wb_summary$Subgroup)
  sec_grobs <- lapply(sg_names, make_wb_section)
  
  main_title <- textGrob(
    paste0("Measurement Type Moderation — ", dataset_tag),
    gp = gpar(fontsize = 13, fontface = "bold")
  )
  
  spacer <- rectGrob(gp = gpar(col = NA, fill = NA),
                     height = unit(0.3, "inches"))
  
  all_grobs <- c(list(main_title),
                 do.call(c, lapply(sec_grobs, function(g) list(spacer, g))))
  
  total_h <- 1.2 + length(sg_names) *
    (0.7 + max(table(wb_summary$Subgroup)) * 0.35)
  
  pdf(wb_pdf_out, width = 9, height = max(5, total_h))
  grid.newpage()
  do.call(grid.arrange, c(all_grobs, list(ncol = 1, padding = unit(0.4, "inches"))))
  dev.off()
  
  cat("WB_measure summary PDF saved to", wb_pdf_out, "\n")
  }   # end if(FALSE) WB_measure summary PDF
}


### 9) Moderation summary table (PDF + CSV) ---------------------------
#
#  Combines results from all three continuous moderators (time lag, % female,
#  mean age) into a single formatted table matching the structure:
#    Moderator | Subgroup | t-stat | p-value | Significant?

cat("\n=== Saving moderation summary table ===\n")

if (length(moderation_results_all) > 0) {

  mod_summary <- do.call(rbind, moderation_results_all) %>%
    mutate(
      t_label   = ifelse(is.na(t_stat), "\u2014",
                         paste0("t(", round(df), ") = ", formatC(t_stat, digits = 2, format = "f"))),
      p_label   = ifelse(is.na(p_val),  "—",
                  ifelse(p_val < .001,   "p < .001",
                         paste0("p = ", formatC(p_val, digits = 3, format = "f")))),
      Significant = ifelse(is.na(p_val), "—",
                    ifelse(p_val < .05,  "Yes",
                    ifelse(p_val < .10,  "Marginal", "No")))
    ) %>%
    select(Moderator, Subgroup, `t-stat` = t_label, `p-value` = p_label, Significant)

  # ---- Save as CSV ----
  csv_out <- paste0("moderation_summary_", dataset_tag, ".csv")
  # DISABLED (moderation summary not in keep-list):
  # write_csv(mod_summary, csv_out)
  # cat("Moderation summary CSV saved to", csv_out, "\n")

  # ---- Save as PDF table ----
  if (FALSE) {   # DISABLED: moderation summary PDF not in keep-list
  if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
  if (!requireNamespace("grid",      quietly = TRUE)) install.packages("grid")
  library(gridExtra)
  library(grid)

  pdf_tbl_out <- paste0("moderation_summary_", dataset_tag, ".pdf")

  # Split by moderator for a section-per-moderator layout
  moderator_names <- unique(mod_summary$Moderator)

  # Build one grob per moderator section
  make_section_grob <- function(mod_name) {
    df_sec <- mod_summary[mod_summary$Moderator == mod_name,
                          c("Subgroup", "t-stat", "p-value", "Significant")]

    # Colour-code the Significant column
    sig_colours <- ifelse(df_sec$Significant == "Yes",      "#1a7a1a",
                   ifelse(df_sec$Significant == "Marginal", "#b36200", "#333333"))

    tt <- ttheme_minimal(
      core    = list(fg_params  = list(hjust = 0, x = 0.05, fontsize = 9),
                     bg_params  = list(fill = c("white", "#f5f5f5"))),
      colhead = list(fg_params  = list(hjust = 0, x = 0.05, fontsize = 9,
                                       fontface = "bold"),
                     bg_params  = list(fill = "#2c2c5e", col = NA),
                     fg_params2 = list(col = "white"))
    )

    tbl_grob <- tableGrob(df_sec, rows = NULL, theme = tt)

    # Colour the Significant column text
    sig_col_idx <- which(colnames(df_sec) == "Significant")
    for (i in seq_len(nrow(df_sec))) {
      cell_name <- paste0("core-fg-", i, "-", sig_col_idx)
      idx <- which(tbl_grob$layout$name == cell_name)
      if (length(idx) > 0)
        tbl_grob$grobs[[idx]]$gp$col <- sig_colours[i]
    }

    title_grob <- textGrob(mod_name, gp = gpar(fontsize = 11, fontface = "bold"),
                           hjust = 0, x = 0)
    arrangeGrob(title_grob, tbl_grob, ncol = 1,
                heights = unit(c(0.4, nrow(df_sec) * 0.35 + 0.4), "inches"))
  }

  section_grobs <- lapply(moderator_names, make_section_grob)

  # Stack all sections with spacers and a main title
  main_title <- textGrob(
    paste0("Moderation Analyses — ", dataset_tag),
    gp = gpar(fontsize = 13, fontface = "bold")
  )

  spacer <- rectGrob(gp = gpar(col = NA, fill = NA), height = unit(0.25, "inches"))

  all_grobs <- c(list(main_title),
                 do.call(c, lapply(section_grobs, function(g) list(spacer, g))))

  total_h <- 1 + length(moderator_names) *
             (0.5 + max(table(mod_summary$Moderator)) * 0.35 + 0.6)

  pdf(pdf_tbl_out, width = 9, height = max(6, total_h))
  grid.newpage()
  do.call(grid.arrange, c(all_grobs,
                          list(ncol   = 1,
                               padding = unit(0.5, "inches"))))
  dev.off()
  cat("Moderation summary PDF saved to", pdf_tbl_out, "\n")
  }   # end if(FALSE) moderation summary PDF

} else {
  cat("  No moderation results collected — summary table skipped.\n")
}


### 10) Export effect size dataset + weights --------------------------

# ---- Collect this category's row for the combined cross-category table ----
.mr <- if (length(moderation_results_all)) do.call(rbind, moderation_results_all) else
       data.frame(Moderator = character(), t_stat = numeric(), p_val = numeric())
.pick <- function(lbl) {
  hit <- .mr[.mr$Moderator == lbl, , drop = FALSE]
  if (nrow(hit) == 0) c(t = NA_real_, p = NA_real_, df = NA_real_) else
    c(t = hit$t_stat[1], p = hit$p_val[1], df = hit$df[1])
}
.age <- .pick("Mean age (years)"); .gen <- .pick("% Female")
.tl  <- .pick("Time lag (months)"); .rob <- .pick("Risk of bias")
overall_table_rows[[current_category]] <- data.frame(
  Category = current_category,
  age_t = .age["t"], age_p = .age["p"], age_df = .age["df"],
  gen_t = .gen["t"], gen_p = .gen["p"], gen_df = .gen["df"],
  tl_t  = .tl["t"],  tl_p  = .tl["p"],  tl_df  = .tl["df"],
  rob_t = .rob["t"], rob_p = .rob["p"], rob_df = .rob["df"],
  wb_Q  = wb_omni_Q, wb_df = wb_omni_df, wb_p = wb_omni_p,
  row.names = NULL, stringsAsFactors = FALSE
)

wts <- tryCatch({
  if (inherits(res_main, "rma.mv")) as.numeric(diag(res_main$W))
  else weights(res_main)
}, error = function(e) rep(NA_real_, nrow(escalc_dat)))

summary_table <- tibble(
  study_ID      = escalc_dat$study_ID,
  cluster_id    = escalc_dat$cluster_id,
  Citation      = escalc_dat$Citation,
  effect_id     = escalc_dat$effect_id,
  risk_factor   = escalc_dat$risk_factor,
  r             = escalc_dat$r,
  n             = escalc_dat$n,
  female        = escalc_dat$female,
  time_lag      = escalc_dat$time_lag,
  WB_construct  = escalc_dat$WB_construct,
  WB_measure    = escalc_dat$WB_measure,
  yi_z          = escalc_dat$yi,
  vi_z          = escalc_dat$vi,
  weight_row    = wts,
  p_value_clean = escalc_dat$p_value_clean
)

# DISABLED (effect-size CSVs not in keep-list):
# write_csv(summary_table, paste0("effect_sizes_", dataset_tag, ".csv"))
# cat("\nSaved effect size table to", paste0("effect_sizes_", dataset_tag, ".csv"), "\n")

  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("\n*** ERROR while analysing '", current_category,
        "' — skipping to the next category. ***\n", sep = "")
    cat("    Message:", conditionMessage(e), "\n")
  })

  # Reset graphics state before the next category
  graphics.off()
  par(mar = c(5.1, 4.1, 4.1, 2.1))

}  # end for (current_category)


### 10B) Funnel / Egger asymmetry summary table ----------------------
# One row per category that met the k>=10 / clusters>=5 threshold. Skipped
# categories simply don't appear (ml_egger returned NULL for them).
egger_table <- do.call(rbind, egger_all)
if (!is.null(egger_table) && nrow(egger_table) > 0) {
  rownames(egger_table) <- NULL
  egger_summary <- data.frame(
    `Antecedent group` = egger_table$analysis,
    `k (effects)`      = egger_table$k_effects,
    `Clusters`         = egger_table$n_clusters,
    `Egger slope [95% CI]` = sprintf("%.3f [%.3f, %.3f]",
                              egger_table$slope, egger_table$ci_lb, egger_table$ci_ub),
    `p`                = formatC(egger_table$p, format = "f", digits = 3),
    `Asymmetry?`       = ifelse(egger_table$ci_lb > 0 | egger_table$ci_ub < 0, "yes", "no"),
    `Adjusted r (PET)` = formatC(egger_table$intercept_r, format = "f", digits = 3),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  cat("\n\n=== Funnel / Egger asymmetry summary ===\n")
  print(egger_summary, row.names = FALSE)
  write_csv(egger_summary, "egger_asymmetry_summary.csv")
  cat("\nSaved egger_asymmetry_summary.csv\n")
} else {
  cat("\nNo categories met the k>=10 / clusters>=5 threshold — no Egger summary table.\n")
}


### 11) COMBINED moderation table across all categories --------------
#  One row per antecedent (category), grouped into families, with columns for
#  Age, Gender (woman), Time lag, and Measurement type (Omnibus SI vs Q) —
#  matching the manuscript table. Significant cells (p < .05) are shown bold
#  in the PDF. (Risk-of-bias t/p are also written to the CSV for reference.)

cat("\n=== Building combined moderation table across all categories ===\n")

if (length(overall_table_rows) > 0) {

  combined <- do.call(rbind, overall_table_rows)
  rownames(combined) <- NULL

  # Map each category -> "Family|Antecedent display label" (edit labels freely)
  family_map <- c(
    "Gender"                              = "Individual factors|Gender/sex",
    "Psychological distress"              = "Individual factors|Psychological distress",
    "Bullied in adolescence" = "Prior bullying victimisation|In adolescence",
    "Bullied at work"                     = "Prior bullying victimisation|At work",
    "Negative work climate"               = "Work social environment|Negative work climate",
    "Positive leadership"                 = "Work social environment|Positive leadership",
    "Positive relationships at work"      = "Work social environment|Positive relationships at work",
    "Job control"                         = "Job and role characteristics|Job control",
    "Job demands"                         = "Job and role characteristics|Job demands",
    "Role ambiguity and conflict"         = "Job and role characteristics|Role ambiguity and conflict",
    "Agreeableness"                       = "Personality|Agreeableness",
    "Openness"                            = "Personality|Openness",
    "Conscientiousness"                   = "Personality|Conscientiousness",
    "Extraversion"                        = "Personality|Extraversion",
    "Neuroticism"                         = "Personality|Neuroticism"
  )
  fam_lookup <- family_map[combined$Category]
  fam_lookup[is.na(fam_lookup)] <- paste0("Other|", combined$Category[is.na(fam_lookup)])
  combined$Family     <- sub("\\|.*$", "", fam_lookup)
  combined$Antecedent <- sub("^.*\\|", "", fam_lookup)

  # Keep the order implied by target_categories
  combined <- combined[order(match(combined$Category, target_categories)), ]

  fmt_tp <- function(t, p, df)
    ifelse(is.na(t) | is.na(p), "\u2014",
      paste0("t(", round(df), ") = ", formatC(t, format = "f", digits = 2), ", ",
             ifelse(p < .001, "p < .001",
                    paste0("p = ", formatC(p, format = "f", digits = 3)))))

  fmt_Q <- function(Q, df, p)
    ifelse(is.na(Q) | is.na(p), "\u2014",
      paste0("Q(", round(df), ") = ", formatC(Q, format = "f", digits = 2), ", ",
             ifelse(p < .001, "p < .001",
                    paste0("p = ", formatC(p, format = "f", digits = 3)))))

  tbl <- data.frame(
    Family     = combined$Family,
    Antecedent = combined$Antecedent,
    `Age`                          = fmt_tp(combined$age_t, combined$age_p, combined$age_df),
    `Gender (woman)`               = fmt_tp(combined$gen_t, combined$gen_p, combined$gen_df),
    `Time lag`                     = fmt_tp(combined$tl_t,  combined$tl_p,  combined$tl_df),
    `Risk of bias`                 = fmt_tp(combined$rob_t, combined$rob_p, combined$rob_df),
    `Measurement type (Omnibus SL vs Q)` = fmt_Q(combined$wb_Q, combined$wb_df, combined$wb_p),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  sig <- cbind(
    !is.na(combined$age_p) & combined$age_p < .05,
    !is.na(combined$gen_p) & combined$gen_p < .05,
    !is.na(combined$tl_p)  & combined$tl_p  < .05,
    !is.na(combined$rob_p) & combined$rob_p < .05,
    !is.na(combined$wb_p)  & combined$wb_p  < .05
  )

  # CSV: formatted cells plus the raw numbers (incl. risk-of-bias) for reference
  write_csv(
    cbind(tbl, combined[, c("age_t","age_p","age_df","gen_t","gen_p","gen_df",
                            "tl_t","tl_p","tl_df","rob_t","rob_p","rob_df",
                            "wb_Q","wb_df","wb_p")]),
    "combined_moderation_table.csv")
  cat("Saved combined_moderation_table.csv\n")

  # ---- Formatted PDF, grouped by family, bold significant cells ----
  tryCatch({
    if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
    library(gridExtra); library(grid)

    fams <- unique(tbl$Family)
    disp_rows <- list(); sig_rows <- list(); is_header <- logical(0)
    for (fm in fams) {
      disp_rows[[length(disp_rows) + 1]] <- c(fm, "", "", "", "", "")
      sig_rows[[length(sig_rows) + 1]]   <- rep(FALSE, 5)
      is_header <- c(is_header, TRUE)
      for (i in which(tbl$Family == fm)) {
        disp_rows[[length(disp_rows) + 1]] <- c(
          tbl$Antecedent[i], tbl$Age[i], tbl$`Gender (woman)`[i],
          tbl$`Time lag`[i], tbl$`Risk of bias`[i],
          tbl$`Measurement type (Omnibus SL vs Q)`[i])
        sig_rows[[length(sig_rows) + 1]] <- as.logical(sig[i, ])
        is_header <- c(is_header, FALSE)
      }
    }
    disp <- as.data.frame(do.call(rbind, disp_rows), stringsAsFactors = FALSE)
    colnames(disp) <- c("Antecedent", "Age", "Gender (woman)",
                        "Time lag", "Risk of bias", "Measurement type (Omnibus SL vs Q)")

    tt <- ttheme_minimal(
      core    = list(fg_params = list(hjust = 0, x = 0.03, fontsize = 8.5),
                     bg_params = list(fill = "white")),
      colhead = list(fg_params = list(fontsize = 9, fontface = "bold"),
                     bg_params = list(fill = "grey92", col = NA)))
    g <- tableGrob(disp, rows = NULL, theme = tt)

    for (r in seq_len(nrow(disp))) {
      if (is_header[r]) {
        for (cc in seq_len(ncol(disp))) {
          ii <- which(g$layout$name == paste0("core-fg-", r, "-", cc))
          if (length(ii)) g$grobs[[ii]]$gp <- gpar(fontface = "bold", fontsize = 9)
          jj <- which(g$layout$name == paste0("core-bg-", r, "-", cc))
          if (length(jj)) g$grobs[[jj]]$gp <- gpar(fill = "grey85", col = NA)
        }
      } else {
        for (cc in 2:ncol(disp)) if (isTRUE(sig_rows[[r]][cc - 1])) {
          ii <- which(g$layout$name == paste0("core-fg-", r, "-", cc))
          if (length(ii)) g$grobs[[ii]]$gp <- gpar(fontface = "bold", fontsize = 8.5)
        }
      }
    }

    ttl  <- textGrob("Moderation analyses across mini meta-analyses",
                     gp = gpar(fontsize = 12, fontface = "bold"), hjust = 0, x = 0.02)
    foot <- textGrob("SI = self-labelling items; Q = validated questionnaire.   Bold = p < .05.",
                     gp = gpar(fontsize = 7.5, fontface = "italic"), hjust = 0, x = 0.02)
    pdf("combined_moderation_table.pdf",
        width = 12.5, height = max(4, 0.42 * (nrow(disp) + 3)))
    grid.arrange(ttl, g, foot, ncol = 1,
                 heights = unit.c(unit(0.4, "in"),
                                  unit(0.30 * (nrow(disp) + 1), "in"),
                                  unit(0.35, "in")))
    dev.off()
    cat("Saved combined_moderation_table.pdf\n")
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("*** Combined-table PDF error:", conditionMessage(e),
        " (CSV was still written).\n")
  })

} else {
  cat("  No category results collected — combined table skipped.\n")
}

### 12) GROUPED (combined) forest figures ----------------------------
#  Stacks several categories into single condensed figures (see forest_groups
#  near the top). Each subgroup keeps its own diamond + heterogeneity line.

cat("\n=== Building grouped (combined) forest figures ===\n")
for (grp_name in names(forest_groups)) {
  # Keep only these four grouped figures; "Individual factors" is exported as the
  # per-category forest_plot_Gender / forest_plot_Psychological_distress PDFs instead.
  if (!grp_name %in% c("Individual factors", "Prior bullying victimisation", "Personality",
                       "Work social environment", "Job and role characteristics")) {
    cat("  (skipping grouped figure '", grp_name, "' — not in keep-list)\n", sep = ""); next
  }
  members <- forest_groups[[grp_name]]
  sub <- escalc_full[tolower(trimws(as.character(escalc_full$Meta_analysis))) %in%
                       tolower(trimws(members)), , drop = FALSE]
  if (nrow(sub) < 2) {
    cat("  Skipping group '", grp_name, "' — not enough rows.\n", sep = ""); next
  }
  tag <- gsub("(^_|_$)", "", gsub("[^A-Za-z0-9]+", "_", grp_name))
  tryCatch(
    make_grouped_forest(sub, subgroup_order = members, out_tag = tag, title = grp_name),
    error = function(e) {
      if (dev.cur() > 1) dev.off()
      cat("  *** grouped forest error for '", grp_name, "': ",
          conditionMessage(e), "\n", sep = "")
    })
}


cat("\n=== All", length(target_categories), "category analyses complete ===\n")


### ====================================================================
### 13) SENSITIVITY ANALYSIS — standardised betas ----------------------
### ====================================================================
#  Pools ONLY rows flagged Sensitivity == "yes" (the Include? column is
#  ignored for these rows — the flag overrides it), one mini
#  meta-analysis per Meta_analysis category, mirroring the main models:
#  rma.mv with Citation/Cohort clustering when clusters repeat,
#  otherwise rma; REML; pooled on Fisher's z and back-transformed to r.
#
#  Effect sizes come from sensitivity_effect (a standardised beta) and
#  are converted to r — by default with the Peterson & Brown (2005)
#  adjustment (see the OPTIONS block near the top of the script) — and
#  then to Fisher's z, so the pooled results are directly comparable
#  with the main analyses.
#
#  sensitivity_p may contain either a p-value ("0.012", "<.001", "n.s.")
#  or a 95% CI for the beta ("0.05, 0.33"). Depending on
#  sens_variance_source it is either the variance source ("p_ci") or a
#  printed consistency check against the n-based SE ("n", the default).
#
#  Outputs (all prefixed beta_sensitivity_ / SENSITIVITY_BETA so nothing
#  from the main analyses is overwritten):
#    beta_sensitivity_effect_sizes.csv     parsed + converted data, SE check
#    beta_sensitivity_results.csv          pooled result per category
#    beta_sensitivity_comparison.csv/.pdf  main vs sensitivity per category
#    beta_sensitivity_loo_<category>.csv   leave-one-out per category
#    forest_grouped_SENSITIVITY_BETA.pdf   stacked forest, all categories

if (isTRUE(run_beta_sensitivity)) {

cat("\n\n############################################################\n")
cat("###  SECTION 13: SENSITIVITY ANALYSIS (standardised betas)\n")
sens_es_label <- if (isTRUE(sens_use_peterson_brown)) "r" else "\u03b2"
cat("############################################################\n")

tryCatch({

  dS <- dat_raw
  .nnm <- tolower(gsub("[ _.?]", "", names(dS)))
  .colS <- function(key) { h <- names(dS)[.nnm == key]; if (length(h)) h[1] else NA_character_ }

  C_flag <- .colS("sensitivity")
  C_eff  <- .colS("sensitivityeffect")
  C_sp   <- .colS("sensitivityp")
  C_ma   <- .colS("metaanalysis")
  C_cit  <- if (!is.na(.colS("citation"))) .colS("citation") else .colS("citations")
  C_n    <- .colS("n")
  C_sid  <- .colS("studyid")
  C_wb   <- .colS("wbmeasure")
  C_tl   <- .colS("timelag")

  if (is.na(C_flag) || is.na(C_eff))
    stop("Could not find the 'Sensitivity' and/or 'sensitivity_effect' columns. ",
         "Columns present: ", paste(names(dS), collapse = ", "))
  if (is.na(C_ma) || is.na(C_cit))
    stop("Could not find the Meta_analysis and/or Citation column in the raw data.")
  if (is.na(C_sp))
    cat("  NOTE: no 'sensitivity_p' column found — p/CI parsing skipped.\n")

  # ---- Select flagged rows (Sensitivity == yes OVERRIDES Include?) ----
  .keepS <- tolower(trimws(as.character(dS[[C_flag]]))) %in% c("yes", "y", "true", "1")
  cat("  Rows flagged Sensitivity == yes:", sum(.keepS), "of", nrow(dS), "\n")
  dS <- dS[.keepS, , drop = FALSE]
  if (nrow(dS) == 0) stop("No rows are flagged 'yes' in the Sensitivity column.")

  # ---- Parse beta and n ----
  betaS <- suppressWarnings(as.numeric(as.character(dS[[C_eff]])))
  nS    <- if (!is.na(C_n)) suppressWarnings(as.integer(as.character(dS[[C_n]])))
           else rep(NA_integer_, nrow(dS))

  # ---- Parse sensitivity_p: one number -> two-sided p; two numbers -> 95% CI ----
  # ("<.001" uses the bound as the p-value, noted; "n.s." yields no SE.)
  parse_sens_p <- function(s, b) {
    out <- list(kind = "none", se_r = NA_real_, note = "")
    s <- trimws(as.character(s))
    if (is.na(s) || s == "" || tolower(s) %in% c("na", "n/a", "-")) return(out)
    if (grepl("n\\.?\\s?s", tolower(s))) { out$note <- "n.s. — no SE derivable"; return(out) }
    nums <- suppressWarnings(as.numeric(unlist(
      regmatches(s, gregexpr("-?[0-9]*\\.?[0-9]+", s)))))
    nums <- nums[!is.na(nums)]
    if (length(nums) >= 2) {                              # 95% CI bounds for beta
      lb <- min(nums[1:2]); ub <- max(nums[1:2])
      out$kind <- "ci"; out$se_r <- (ub - lb) / (2 * qnorm(0.975))
      if (!is.na(b) && (b < lb || b > ub))
        out$note <- "WARNING: beta lies outside its stated CI — check this cell"
    } else if (length(nums) == 1) {                       # two-sided p-value
      p <- nums[1]
      if (!is.na(b) && p > 0 && p < 1) {
        out$kind <- "p"; out$se_r <- abs(b) / qnorm(1 - p / 2)
        if (grepl("<", s)) out$note <- "p given as a bound ('<'); SE derived from that bound"
      } else out$note <- "p-value outside (0,1) — ignored"
    }
    out
  }
  pinfoS <- lapply(seq_len(nrow(dS)), function(i)
    if (!is.na(C_sp)) parse_sens_p(dS[[C_sp]][i], betaS[i])
    else list(kind = "none", se_r = NA_real_, note = ""))

  # ---- Convert beta -> r -> Fisher's z ----
  n_big_beta <- sum(abs(betaS) > 0.5, na.rm = TRUE)
  if (isTRUE(sens_use_peterson_brown)) {
    rS <- betaS + 0.05 * ifelse(betaS >= 0, 1, 0)   # Peterson & Brown (2005)
    cat("  Conversion: Peterson-Brown adjustment (r = beta + .05*lambda).\n")
  } else {
    rS <- betaS
    cat("  Conversion: beta treated directly as r (no adjustment).\n")
  }
  if (n_big_beta > 0)
    cat("  NOTE:", n_big_beta, "beta value(s) outside [-0.5, 0.5] — the beta->r",
        "approximation is less accurate there; interpret with caution.\n")
  rS <- pmin(pmax(rS, -0.9999), 0.9999)

  # ---- Sampling variances (Fisher-z scale) ----
  vi_nS   <- ifelse(!is.na(nS) & nS > 3, 1 / (nS - 3), NA_real_)
  se_rpci <- vapply(pinfoS, function(x) x$se_r, numeric(1))
  # delta method: var(z) ~= var(r) / (1 - r^2)^2
  vi_pciS <- ifelse(!is.na(se_rpci), (se_rpci^2) / (1 - rS^2)^2, NA_real_)

  if (identical(sens_variance_source, "p_ci")) {
    viS <- ifelse(!is.na(vi_pciS), vi_pciS, vi_nS)
    vsrc <- ifelse(!is.na(vi_pciS), "p/CI", "n (fallback)")
    .nfb <- sum(is.na(vi_pciS) & !is.na(vi_nS))
    if (.nfb > 0) cat("  NOTE:", .nfb, "row(s) had no usable p/CI — variance fell back to n.\n")
  } else {
    viS <- vi_nS
    vsrc <- rep("n", nrow(dS))
  }

  # ---- Assemble working dataset, then drop unusable rows ----
  sens_full <- data.frame(
    study_ID      = if (!is.na(C_sid)) as.character(dS[[C_sid]]) else as.character(dS[[C_cit]]),
    Citation      = as.character(dS[[C_cit]]),
    Meta_analysis = trimws(as.character(dS[[C_ma]])),
    n             = nS,
    beta          = betaS,
    r_converted   = rS,
    yi            = atanh(rS),
    vi            = viS,
    se_z_n        = sqrt(vi_nS),
    se_z_pci      = sqrt(vi_pciS),
    variance_src  = vsrc,
    pci_kind      = vapply(pinfoS, function(x) x$kind, character(1)),
    pci_note      = vapply(pinfoS, function(x) x$note, character(1)),
    WB_measure    = if (!is.na(C_wb)) as.character(dS[[C_wb]]) else "",
    time_lag      = if (!is.na(C_tl)) suppressWarnings(as.numeric(as.character(dS[[C_tl]]))) else NA_real_,
    stringsAsFactors = FALSE
  )
  # Carry the optional display columns through for the forest plot
  for (.opt in c("samplenotes", "riskfactordisplay", "riskfactor", "cohort")) {
    .h <- .colS(.opt)
    if (!is.na(.h)) sens_full[[.h]] <- dS[[.h]]
  }

  .dropS <- is.na(sens_full$yi) | is.na(sens_full$vi) | sens_full$Meta_analysis == ""
  if (any(.dropS)) {
    cat("  WARNING — dropping", sum(.dropS), "flagged row(s) with missing beta,",
        "variance, or Meta_analysis label:\n   ",
        paste(unique(sens_full$Citation[.dropS]), collapse = " | "), "\n")
    sens_full <- sens_full[!.dropS, , drop = FALSE]
  }
  if (nrow(sens_full) == 0) stop("No usable sensitivity rows after cleaning.")

  # ---- Consistency check: n-based vs p/CI-implied SE (Fisher-z scale) ----
  .chk <- !is.na(sens_full$se_z_n) & !is.na(sens_full$se_z_pci)
  if (any(.chk)) {
    cat("\n  Consistency check — SE(z) from n vs SE(z) implied by sensitivity_p:\n")
    .chk_df <- data.frame(Citation = sens_full$Citation[.chk],
                          Category = sens_full$Meta_analysis[.chk],
                          se_n     = round(sens_full$se_z_n[.chk], 4),
                          se_p_ci  = round(sens_full$se_z_pci[.chk], 4),
                          ratio    = round(sens_full$se_z_pci[.chk] / sens_full$se_z_n[.chk], 2),
                          source   = sens_full$pci_kind[.chk])
    print(.chk_df, row.names = FALSE)
    if (any(.chk_df$ratio > 2 | .chk_df$ratio < 0.5, na.rm = TRUE))
      cat("  NOTE: ratios far from 1 usually mean the reported p/CI reflects a\n",
          " covariate-adjusted SE (expected for betas) or a transcription issue.\n")
  }
  .notes <- sens_full$pci_note != ""
  if (any(.notes))
    for (i in which(.notes))
      cat("  [", sens_full$Citation[i], "] ", sens_full$pci_note[i], "\n", sep = "")

  # ---- Repeated-sample clusters (same rule as the main analyses) ----
  .cohS <- if (length(cohort_col) && cohort_col[1] %in% names(sens_full))
    as.character(sens_full[[cohort_col[1]]]) else rep("", nrow(sens_full))
  sens_full$cluster_id <- make_cluster_id(sens_full$Citation, .cohS, ignore_registry)

  # ---- Save the parsed/converted effect sizes ----
  write_csv(sens_full, "beta_sensitivity_effect_sizes.csv")
  cat("\n  Saved beta_sensitivity_effect_sizes.csv\n")

  # ---- One mini meta-analysis per category present in the flagged rows ----
  sens_cats <- target_categories[tolower(trimws(target_categories)) %in%
                                   unique(tolower(sens_full$Meta_analysis))]
  .extraS <- setdiff(unique(tolower(sens_full$Meta_analysis)),
                     tolower(trimws(target_categories)))
  if (length(.extraS))
    cat("  NOTE — flagged rows have categories not in target_categories",
        "(they will be skipped):", paste(.extraS, collapse = " | "), "\n")

  sens_pooled <- list()

  for (sc in sens_cats) {
    dC <- sens_full[tolower(sens_full$Meta_analysis) == tolower(trimws(sc)), , drop = FALSE]
    dC$effect_id <- seq_len(nrow(dC))
    tagS <- gsub("(^_|_$)", "", gsub("[^A-Za-z0-9]+", "_", trimws(sc)))

    cat("\n--- Beta sensitivity meta-analysis:", sc, "(k =", nrow(dC), ") ---\n")

    if (nrow(dC) == 1) {
      .lb <- tanh(dC$yi - qnorm(0.975) * sqrt(dC$vi))
      .ub <- tanh(dC$yi + qnorm(0.975) * sqrt(dC$vi))
      cat("  Only one effect — no pooling. r =", round(tanh(dC$yi), 3),
          " [", round(.lb, 3), ",", round(.ub, 3), "]\n")
      sens_pooled[[sc]] <- data.frame(Category = sc, k = 1L, r = tanh(dC$yi),
                                      lb = .lb, ub = .ub, model = "single effect",
                                      stringsAsFactors = FALSE)
      next
    }

    .dupC <- any(duplicated(dC$cluster_id))
    fitS <- tryCatch({
      if (.dupC) rma.mv(yi, vi, random = ~ 1 | cluster_id/effect_id,
                        data = dC, method = "REML")
      else       rma(yi, vi, data = dC, method = "REML")
    }, error = function(e) { cat("  Model failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(fitS)) next

    print(fitS)
    if (inherits(fitS, "rma.mv")) {
      cat("  Cluster-robust (CR2) test for the pooled effect:\n")
      .cr <- tryCatch(clubSandwich::coef_test(fitS, vcov = "CR2", cluster = dC$cluster_id),
                      error = function(e) NULL)
      if (!is.null(.cr)) print(.cr) else cat("  (CR2 unavailable for this fit.)\n")
    }

    .zb <- as.numeric(fitS$b[1, 1])
    cat("  Pooled r =", round(tanh(.zb), 3),
        " [", round(tanh(fitS$ci.lb), 3), ",", round(tanh(fitS$ci.ub), 3), "]\n")
    sens_pooled[[sc]] <- data.frame(Category = sc, k = fitS$k, r = tanh(.zb),
                                    lb = tanh(fitS$ci.lb), ub = tanh(fitS$ci.ub),
                                    model = if (.dupC) "multilevel" else "random-effects",
                                    stringsAsFactors = FALSE)

    # ---- Leave-one-out (only meaningful with 3+ studies) ----
    .sids <- unique(dC$study_ID)
    if (length(.sids) > 2) {
      .loo <- do.call(rbind, lapply(.sids, function(sid) {
        tmp <- dC[dC$study_ID != sid, ]
        f <- tryCatch({
          if (any(duplicated(tmp$cluster_id)))
            rma.mv(yi, vi, random = ~ 1 | cluster_id/effect_id, data = tmp, method = "REML")
          else rma(yi, vi, data = tmp, method = "REML")
        }, error = function(e) NULL)
        if (is.null(f)) data.frame(study_omitted = sid, r_pooled = NA_real_,
                                   r_lb = NA_real_, r_ub = NA_real_)
        else data.frame(study_omitted = sid,
                        r_pooled = tanh(as.numeric(f$b[1, 1])),
                        r_lb = tanh(f$ci.lb), r_ub = tanh(f$ci.ub))
      }))
      write_csv(.loo, paste0("beta_sensitivity_loo_", tagS, ".csv"))
      cat("  Saved beta_sensitivity_loo_", tagS, ".csv\n", sep = "")
    }
  }

  # ---- Pooled results + comparison with the MAIN analyses ----
  if (length(sens_pooled)) {
    sensP <- do.call(rbind, sens_pooled); rownames(sensP) <- NULL
    write_csv(sensP, "beta_sensitivity_results.csv")
    cat("\n  Saved beta_sensitivity_results.csv\n")

    mainP <- if (length(overall_pooled)) do.call(rbind, overall_pooled) else NULL
    .fmtP <- function(r, lb, ub) sprintf("%.2f [%.2f, %.2f]", r, lb, ub)

    comp <- data.frame(
      Antecedent                   = sensP$Category,
      `Main k`                     = NA_integer_,
      `Main r [95% CI]`            = "\u2014",
      `Sensitivity k`              = sensP$k,
      `Sensitivity r [95% CI]`     = .fmtP(sensP$r, sensP$lb, sensP$ub),
      check.names = FALSE, stringsAsFactors = FALSE)
    names(comp)[names(comp) == "Sensitivity r [95% CI]"] <-
      paste0("Sensitivity ", sens_es_label, " [95% CI]")
    if (!is.null(mainP)) {
      .m  <- match(sensP$Category, mainP$Category)
      .ok <- !is.na(.m)
      comp$`Main k`[.ok]          <- mainP$k[.m[.ok]]
      comp$`Main r [95% CI]`[.ok] <- .fmtP(mainP$r[.m[.ok]], mainP$lb[.m[.ok]], mainP$ub[.m[.ok]])
    }
    write_csv(comp, "beta_sensitivity_comparison.csv")
    cat("  Saved beta_sensitivity_comparison.csv\n")

    # Formatted PDF of the comparison table (same style as the other tables)
    tryCatch({
      if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
      library(gridExtra); library(grid)
      .tt <- ttheme_minimal(
        core    = list(fg_params = list(hjust = 0, x = 0.04, fontsize = 9),
                       bg_params = list(fill = c("white", "#f5f5f5"))),
        colhead = list(fg_params = list(hjust = 0, x = 0.04, fontsize = 9, fontface = "bold"),
                       bg_params = list(fill = "grey92", col = NA)))
      .g    <- tableGrob(comp, rows = NULL, theme = .tt)
      .ttl  <- textGrob("Sensitivity analysis: standardised betas vs main meta-analyses",
                        gp = gpar(fontsize = 12, fontface = "bold"), hjust = 0, x = 0.02)
      .conv <- if (isTRUE(sens_use_peterson_brown))
        "Betas converted to r with the Peterson & Brown (2005) adjustment." else
        "Betas treated directly as r."
      .vtxt <- if (identical(sens_variance_source, "p_ci"))
        "Sampling variances derived from reported p-values/CIs." else
        "Sampling variances from sample size (Fisher's z), as in the main analyses."
      .foot <- textGrob(paste(.conv, .vtxt),
                        gp = gpar(fontsize = 7.5, fontface = "italic"), hjust = 0, x = 0.02)
      pdf("beta_sensitivity_comparison.pdf", width = 10,
          height = max(3, 0.42 * (nrow(comp) + 3)))
      grid.arrange(.ttl, .g, .foot, ncol = 1,
                   heights = unit.c(unit(0.4, "in"),
                                    unit(0.32 * (nrow(comp) + 1), "in"),
                                    unit(0.35, "in")))
      dev.off()
      cat("  Saved beta_sensitivity_comparison.pdf\n")
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      cat("  (comparison PDF skipped:", conditionMessage(e), "— CSV was written.)\n")
    })
  }

  # ---- Stacked grouped forest of all sensitivity categories ----
  tryCatch(
    make_grouped_forest(sens_full, subgroup_order = sens_cats,
                        out_tag = "SENSITIVITY_BETA",
                        title = if (sens_es_label == "r")
                          "Sensitivity analysis \u2014 standardised betas (converted to r)"
                        else "Sensitivity analysis \u2014 standardised betas",
                        xlab_text = if (sens_es_label == "r") "Correlation (r)"
                        else "Standardised beta (\u03b2)",
                        ci_lab    = paste0(sens_es_label, " (95% CI)")),
    error = function(e) {
      if (dev.cur() > 1) dev.off()
      cat("  *** sensitivity forest error:", conditionMessage(e), "\n")
    })

  cat("\n=== Beta sensitivity analysis complete (",
      length(sens_pooled), "categories ) ===\n")

}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("\n*** ERROR in the beta sensitivity analysis — main results are unaffected. ***\n")
  cat("    Message:", conditionMessage(e), "\n")
})

}  # end if (run_beta_sensitivity)
