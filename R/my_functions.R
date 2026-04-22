# ── load_deg_table ────────────────────────────────────────────────────────────
# This function reads the DEG table from either a local file or a URL.
# It handles both CSV and tab-delimited files, and tries UTF-8 encoding
# first before falling back to latin1 if there are special characters.

#' @title Load DEG Table
#' @description Reads the S3 DEG table from a local file or URL.
#'   The file must have a title row on line 1, a header row on line 2, then data.
#' @param path Full file path or URL to the S3 table (.txt or .csv)
#' @return Data frame with columns: Row, TranscriptID, gene_assignment,
#'   GeneSymbol, RefSeq, padj, logFC
#' @export
#' @examples
#' deg <- load_deg_table("https://raw.githubusercontent.com/hjleon19/collaborative-project/refs/heads/master/S3%20Table.txt")
load_deg_table <- function(path) {

  # Detect file separator based on extension (.csv uses comma, everything else tab)
  ext    <- tolower(tools::file_ext(path))
  sep    <- if (ext == "csv") "," else "\t"

  # Detect whether path is a URL or a local file
  is_url <- grepl("^https?://", path)

  # Internal helper: open connection with specified encoding and read table
  read_data <- function(encoding) {
    con <- if (is_url) url(path, encoding = encoding) else file(path, encoding = encoding)
    on.exit(try(close(con), silent = TRUE))
    utils::read.table(
      con,
      sep          = sep,
      header       = TRUE,
      skip         = 1,       # skip the title row ("S3 Table. Full list of DEGs...")
      quote        = "",
      fill         = TRUE,
      comment.char = ""
    )
  }

  # Try UTF-8 first, fall back to latin1 if it fails
  deg <- tryCatch(
    read_data("UTF-8"),
    error = function(e) read_data("latin1")
  )

  # Standardise column names regardless of what the file header says
  colnames(deg) <- c("Row", "TranscriptID", "gene_assignment",
                     "GeneSymbol", "RefSeq", "padj", "logFC")

  # Convert p-value and fold change to numeric and remove any rows with NAs
  deg$padj  <- as.numeric(deg$padj)
  deg$logFC <- as.numeric(deg$logFC)
  deg       <- deg[!is.na(deg$padj) & !is.na(deg$logFC), ]

  # Print a summary so the user knows the data loaded correctly
  cat("DEG table loaded:", nrow(deg), "genes\n")
  cat("logFC range:", round(range(deg$logFC), 3), "\n")
  cat("padj range:",  round(range(deg$padj),  5), "\n")

  deg
}


# ── volcano_plot ──────────────────────────────────────────────────────────────
# This function creates a volcano plot from the DEG table.
# Each point is a gene. X-axis = log2 fold change (how much expression changed).
# Y-axis = -log10(FDR p-value) so the most significant genes appear at the TOP.
# Red points = up-regulated in vEDS. Blue points = down-regulated in vEDS.

#' @title Volcano Plot
#' @description Creates a volcano plot from differential expression results.
#'   Significant genes are coloured red (up) or blue (down). The y-axis shows
#'   -log10(FDR-adjusted p-value) so the most significant genes appear at the top.
#' @param data Data frame containing logFC and p-value columns
#' @param logFC Column name for log fold change values
#' @param pval Column name for p-values
#' @param alpha FDR significance threshold (default = 0.05)
#' @param lfc_thresh Log fold change threshold (default = 1)
#' @param title Plot title
#' @param already_adjusted Logical. TRUE if p-values are already FDR adjusted (default = FALSE)
#' @return A ggplot2 object
#' @export
#' @examples
#' volcano_plot(deg, logFC = "logFC", pval = "padj", already_adjusted = TRUE)
volcano_plot <- function(data,
                         logFC,
                         pval,
                         alpha            = 0.05,
                         lfc_thresh       = 1,
                         title            = "Volcano Plot",
                         already_adjusted = FALSE) {

  # Convert to data frame and extract the relevant columns
  data <- as.data.frame(data)
  data$logFC_vals <- data[[logFC]]
  data$pval_vals  <- data[[pval]]

  # Remove rows where either value is missing
  data <- data[!is.na(data$logFC_vals) & !is.na(data$pval_vals), ]

  if (nrow(data) == 0) {
    stop("No valid rows after removing NAs from logFC and pval columns.")
  }

  # Apply FDR correction if not already done
  # already_adjusted = TRUE means p-values are already BH-corrected (as in our S3 table)
  data$padj <- if (already_adjusted) {
    data$pval_vals
  } else {
    stats::p.adjust(data$pval_vals, method = "BH")
  }

  # Classify each gene as up-regulated, down-regulated, or not significant (NC)
  data$direction <- "NC"
  data$direction[data$padj < alpha & data$logFC_vals >  lfc_thresh] <- "up"
  data$direction[data$padj < alpha & data$logFC_vals < -lfc_thresh] <- "down"

  # Count significant genes in each direction for the annotation labels
  n_up   <- sum(data$direction == "up")
  n_down <- sum(data$direction == "down")

  cat("Significant genes - up:", n_up, "| down:", n_down, "\n")

  # Keep only significant genes for plotting
  sig_data <- data[data$direction != "NC", ]

  if (nrow(sig_data) == 0) {
    stop("No significant genes found. Try lowering alpha or lfc_thresh.")
  }

  # Compute -log10(padj) for the y-axis so significant genes appear at top
  sig_data$y <- -log10(sig_data$padj)
  y_max <- max(sig_data$y, na.rm = TRUE)
  x_min <- min(sig_data$logFC_vals, na.rm = TRUE)
  x_max <- max(sig_data$logFC_vals, na.rm = TRUE)

  # Build the ggplot2 volcano plot
  ggplot2::ggplot(
    sig_data,
    ggplot2::aes(x = logFC_vals, y = y, color = direction)
  ) +
    # Plot each gene as a semi-transparent point
    ggplot2::geom_point(alpha = 0.7, size = 1.5) +

    # Red for up-regulated, blue for down-regulated
    ggplot2::scale_color_manual(values = c("up" = "red", "down" = "blue")) +

    # Vertical dashed lines marking the fold change threshold
    ggplot2::geom_vline(
      xintercept = c(-lfc_thresh, lfc_thresh),
      linetype   = "dashed",
      color      = "grey40"
    ) +

    # Horizontal dashed line marking the significance threshold
    ggplot2::geom_hline(
      yintercept = -log10(alpha),
      linetype   = "dashed",
      color      = "grey40"
    ) +

    # Annotation showing count of down-regulated genes
    ggplot2::annotate(
      "text",
      x = x_min * 0.7, y = y_max * 0.05,
      label = paste(n_down, "down"),
      color = "blue", size = 4, fontface = "bold"
    ) +

    # Annotation showing count of up-regulated genes
    ggplot2::annotate(
      "text",
      x = x_max * 1.2, y = y_max * 0.05,
      label = paste(n_up, "up"),
      color = "red", size = 4, fontface = "bold"
    ) +

    # Axis labels and title
    ggplot2::labs(
      title = title,
      x     = "Log2 Fold Change",
      y     = "-log10(FDR-adjusted p-value)"
    ) +
    ggplot2::theme_minimal()
}


# ── pathway_analysis ──────────────────────────────────────────────────────────
# This function does three things:
#   1. Prints the top N most up and down regulated genes by fold change
#   2. Maps gene symbols to Entrez IDs (required by KEGG database)
#   3. Runs KEGG enrichment to find which biological pathways are overrepresented
#
# It returns everything the vignette needs to then call pathview directly.
# pathview is called in the vignette rather than here because it requires
# its internal data objects to be loaded via library(pathview) which only
# works reliably when called at the top level, not from inside a package function.

#' @title DEG Summary and Multi-Pathway Analysis
#' @description Prints the top N up and down regulated genes, runs KEGG pathway
#'   enrichment on the full DEG list, and returns the results for visualization.
#' @param deg Data frame from load_deg_table (needs GeneSymbol, logFC, padj columns)
#' @param top_n Number of top up and down regulated genes to print (default = 10)
#' @param n_pathways Number of top enriched pathways to return (default = 3)
#' @param padj_cutoff FDR threshold for filtering DEGs (default = 0.01)
#' @param lfc_cutoff Minimum absolute logFC for filtering DEGs (default = 0.2)
#' @return A list with: top_pathways data frame, gene_data named vector for
#'   pathview, kegg_result full enrichment object
#' @export
#' @examples
#' results <- pathway_analysis(deg)
pathway_analysis <- function(deg,
                             top_n       = 10,
                             n_pathways  = 3,
                             padj_cutoff = 0.01,
                             lfc_cutoff  = 0.2) {

  # Check that all required Bioconductor packages are installed
  for (pkg in c("org.Hs.eg.db", "clusterProfiler", "AnnotationDbi")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Please install ", pkg,
                  " first:\n  BiocManager::install('", pkg, "')"))
    }
  }

  # ── Step 1: Filter DEGs ─────────────────────────────────────────────────────
  # Keep only genes meeting both FDR and fold change thresholds
  # padj < 0.01 = less than 1% chance this is a false positive
  # abs(logFC) > 0.2 = removes trivially small expression changes
  deg_filt <- deg[
    !is.na(deg$padj) &
      deg$padj < padj_cutoff &
      abs(deg$logFC) > lfc_cutoff,
  ]

  cat("Genes after filtering:", nrow(deg_filt), "\n")

  if (nrow(deg_filt) == 0) {
    stop("No genes passed the filter. Try increasing padj_cutoff or lowering lfc_cutoff.")
  }

  # ── Step 2: Print top up and down regulated genes ───────────────────────────
  # Sort descending by logFC for top up-regulated genes
  # Sort ascending by logFC for top down-regulated genes
  top_up   <- deg_filt[order(deg_filt$logFC, decreasing = TRUE),  ][1:min(top_n, nrow(deg_filt)), ]
  top_down <- deg_filt[order(deg_filt$logFC, decreasing = FALSE), ][1:min(top_n, nrow(deg_filt)), ]

  cat("\nTop", top_n, "Up-regulated Genes:\n")
  print(top_up[, c("GeneSymbol", "logFC", "padj")])

  cat("\nTop", top_n, "Down-regulated Genes:\n")
  print(top_down[, c("GeneSymbol", "logFC", "padj")])

  # ── Step 3: Map gene symbols to Entrez IDs ──────────────────────────────────
  # KEGG requires numeric Entrez IDs not gene symbols
  # org.Hs.eg.db is the human genome annotation database for this mapping
  # Non-coding RNAs (SNORD, RNU genes) typically fail to map — this is expected
  map <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys    = deg_filt$GeneSymbol,
    keytype = "SYMBOL",
    columns = "ENTREZID"
  )

  deg_mapped <- merge(deg_filt, map, by.x = "GeneSymbol", by.y = "SYMBOL")
  deg_mapped <- deg_mapped[!is.na(deg_mapped$ENTREZID), ]
  deg_mapped <- deg_mapped[!duplicated(deg_mapped$ENTREZID), ]

  cat("\nGenes mapped to Entrez IDs:", nrow(deg_mapped), "\n")

  if (nrow(deg_mapped) == 0) {
    stop("No genes could be mapped to Entrez IDs.")
  }

  # ── Step 4: Build named vector for pathview ─────────────────────────────────
  # pathview needs a named numeric vector:
  #   names  = Entrez IDs
  #   values = logFC (colours genes on the pathway diagram)
  gene_data        <- deg_mapped$logFC
  names(gene_data) <- deg_mapped$ENTREZID

  # ── Step 5: KEGG pathway enrichment ─────────────────────────────────────────
  # enrichKEGG tests each KEGG pathway to see if our DEGs appear more often
  # than expected by chance (Fisher exact test, BH-corrected for multiple testing)
  cat("\nRunning KEGG pathway enrichment...\n")

  kegg_result <- clusterProfiler::enrichKEGG(
    gene          = deg_mapped$ENTREZID,
    organism      = "hsa",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH"
  )

  if (is.null(kegg_result) || nrow(kegg_result) == 0) {
    stop("No enriched KEGG pathways found. Try relaxing padj_cutoff.")
  }

  # Extract top N pathways ranked by adjusted p-value
  top_pathways <- as.data.frame(kegg_result)[1:min(n_pathways, nrow(kegg_result)), ]

  # Print ranked list so user knows which pathway is #1 vs #3
  cat("\nTop enriched KEGG pathways:\n")
  for (i in seq_len(nrow(top_pathways))) {
    cat(sprintf("  #%d  %s  (%s)  p.adjust = %.4f\n",
                i,
                top_pathways$Description[i],
                top_pathways$ID[i],
                top_pathways$p.adjust[i]))
  }

  # Return everything the vignette needs to run pathview
  invisible(list(
    top_pathways = top_pathways,
    gene_data    = gene_data,
    kegg_result  = kegg_result
  ))
}
