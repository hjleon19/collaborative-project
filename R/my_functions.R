#' @title Load DEG Table

#' @description Reads the S3 DEG table from a tab-delimited or xlsx-exported file.

#'   The file must have a title row on line 1, a header row on line 2, then data.

#'   If the file is still in .xlsx format, open it in Excel and save as .csv or .txt first.

#' @param path Full file path to the S3 table (.csv or tab-delimited .txt)

#' @return Data frame with columns: Row, TranscriptID, gene_assignment,

#'   GeneSymbol, RefSeq, padj, logFC

#' @export

#' @examples

#' deg <- load_deg_table("/Users/nikhilakalapatapu/collaborativeproject/example.data/S3Table.txt")



load_deg_table <- function(path) {



  ext <- tolower(tools::file_ext(path))

  sep <- if (ext == "csv") "," else "\t"



  deg <- tryCatch({

    utils::read.table(

      path,

      sep          = sep,

      header       = TRUE,

      skip         = 1,

      quote        = "",

      fill         = TRUE,

      comment.char = "",

      fileEncoding = "UTF-8"

    )

  }, error = function(e) {

    utils::read.table(

      path,

      sep          = sep,

      header       = TRUE,

      skip         = 1,

      quote        = "",

      fill         = TRUE,

      comment.char = "",

      fileEncoding = "latin1"

    )

  })



  colnames(deg) <- c("Row", "TranscriptID", "gene_assignment",

                     "GeneSymbol", "RefSeq", "padj", "logFC")



  deg$padj  <- as.numeric(deg$padj)

  deg$logFC <- as.numeric(deg$logFC)

  deg <- deg[!is.na(deg$padj) & !is.na(deg$logFC), ]



  cat("DEG table loaded:", nrow(deg), "genes\n")

  cat("logFC range:", round(range(deg$logFC), 3), "\n")

  cat("padj range:",  round(range(deg$padj),  5), "\n")



  deg

}



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



  data <- as.data.frame(data)

  data$logFC_vals <- data[[logFC]]

  data$pval_vals  <- data[[pval]]

  data <- data[!is.na(data$logFC_vals) & !is.na(data$pval_vals), ]



  if (nrow(data) == 0) {

    stop("No valid rows after removing NAs from logFC and pval columns.")

  }



  # Apply FDR adjustment if not already done

  data$padj <- if (already_adjusted) {

    data$pval_vals

  } else {

    stats::p.adjust(data$pval_vals, method = "BH")

  }



  # Classify each gene

  data$direction <- "NC"

  data$direction[data$padj < alpha & data$logFC_vals >  lfc_thresh] <- "up"

  data$direction[data$padj < alpha & data$logFC_vals < -lfc_thresh] <- "down"



  n_up   <- sum(data$direction == "Upregulated")

  n_down <- sum(data$direction == "Downregulated")



  cat("Significant genes - up:", n_up, "| down:", n_down, "\n")



  sig_data <- data[data$direction != "NC", ]



  if (nrow(sig_data) == 0) {

    stop("No significant genes found. Try lowering alpha or lfc_thresh.")

  }



  # -log10 transform for y axis

  sig_data$y <- -log10(sig_data$padj)

  y_max <- max(sig_data$y, na.rm = TRUE)

  x_min <- min(sig_data$logFC_vals, na.rm = TRUE)

  x_max <- max(sig_data$logFC_vals, na.rm = TRUE)



  ggplot2::ggplot(

    sig_data,

    ggplot2::aes(x = logFC_vals, y = y, color = direction)

  ) +

    ggplot2::geom_point(alpha = 0.7, size = 1.5) +

    ggplot2::scale_color_manual(values = c("up" = "red", "down" = "blue")) +

    ggplot2::geom_vline(

      xintercept = c(-lfc_thresh, lfc_thresh),

      linetype   = "dashed",

      color      = "grey40"

    ) +

    ggplot2::geom_hline(

      yintercept = -log10(alpha),

      linetype   = "dashed",

      color      = "grey40"

    ) +

    ggplot2::annotate(

      "text",

      x     = x_min * 0.7,

      y     = y_max * 0.05,

      label = paste(n_down, "down"),

      color = "blue",

      size  = 4,

      fontface = "bold"

    ) +

    ggplot2::annotate(

      "text",

      x     = x_max * 1.2,

      y     = y_max * 0.05,

      label = paste(n_up, "up"),

      color = "red",

      size  = 4,

      fontface = "bold"

    ) +

    ggplot2::labs(

      title = title,

      x     = "Log2 Fold Change",

      y     = "-log10(FDR-adjusted p-value)"

    ) +

    ggplot2::theme_minimal()

}



#' @title Hierarchical Clustering Heatmap

#' @description Produces a heatmap with dendrograms directly from a DEG table.

#'   Uses logFC values to represent expression differences between groups.

#' @param deg Data frame output from load_deg_table (must have GeneSymbol and logFC columns)

#' @param dist_method Distance method passed to dist() (default = "euclidean")

#' @param linkage Linkage method passed to hclust() (default = "complete")

#' @param title Plot title

#' @return Invisibly returns NULL

#' @export

#' @examples

#' group <- c(rep("vEDS", 3), rep("control", 9))

#' heatmap_plot(deg)

heatmap_plot <- function(deg,

                         dist_method = "euclidean",

                         linkage     = "complete",

                         title       = "Hierarchical clustering") {



  # Build a simple expression matrix from logFC

  # vEDS columns get logFC value, control columns get 0

  n_veds   <- 3

  n_ctrl   <- 9

  n_samples <- n_veds + n_ctrl



  expr_matrix <- matrix(0, nrow = nrow(deg), ncol = n_samples)

  rownames(expr_matrix) <- deg$GeneSymbol

  colnames(expr_matrix) <- c(paste0("P", 1:n_veds), paste0("C", 1:n_ctrl))



  # vEDS columns reflect the logFC, controls stay at 0 (baseline)

  for (i in 1:n_veds) {

    expr_matrix[, i] <- deg$logFC

  }



  # Z-score per gene, clamp to +/-2.49

  gene_var <- apply(expr_matrix, 1, stats::var)

  expr_matrix <- expr_matrix[gene_var > 0, ]



  expr_scaled <- t(scale(t(expr_matrix)))

  expr_scaled[expr_scaled >  2.49] <-  2.49

  expr_scaled[expr_scaled < -2.49] <- -2.49



  # Sidebar colours

  ctrl_colors <- c(

    "#FFB3B3", "#FFD9B3", "#FFFAB3", "#B3FFB3",

    "#B3FFF0", "#B3D9FF", "#C9B3FF", "#FFB3F0", "#D4B3FF"

  )

  sidebar <- c(rep("#FF69B4", n_veds), ctrl_colors[1:n_ctrl])



  stats::heatmap(

    expr_scaled,

    distfun       = function(x) stats::dist(x, method = dist_method),

    hclustfun     = function(x) stats::hclust(x, method = linkage),

    col           = grDevices::colorRampPalette(c("blue", "white", "red"))(100),

    ColSideColors = sidebar,

    scale         = "none",

    main          = title,

    labRow        = NA

  )



  invisible(NULL)

}

#' @title Protein-Protein Interaction Network
#' @description Builds a PPI network from DEGs using the STRINGdb Bioconductor
#'   package which downloads and caches the database locally — no live API
#'   connection required after first use.
#' @param deg Data frame from load_deg_table (needs GeneSymbol and logFC columns)
#' @param score_threshold Minimum STRING score 0-1000 (default = 400)
#' @param title Plot title
#' @return Invisibly returns the mapped DEG data frame with STRING IDs
#' @export
#' @examples
#' ppi_network(deg)
ppi_network <- function(deg,
                        score_threshold = 400,
                        title           = "PPI Network of DEGs") {

  # STRINGdb must be installed via BiocManager
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("Please install STRINGdb first:\n  BiocManager::install('STRINGdb')")
  }

  # ── 1. Initialise STRINGdb (downloads DB to cache on first run) ───────────
  cat("Initialising STRINGdb (may download on first run)...\n")
  string_db <- STRINGdb::STRINGdb$new(
    version          = "11.5",
    species          = 9606,        # Homo sapiens
    score_threshold  = score_threshold,
    input_directory  = ""           # uses temp cache
  )

  # ── 2. Map gene symbols to STRING IDs ─────────────────────────────────────
  deg_mapped <- string_db$map(
    my_data_frame        = deg,
    my_data_frame_id_col = "GeneSymbol",
    removeUnmappedRows   = TRUE
  )

  cat("Genes mapped to STRING:", nrow(deg_mapped), "of", nrow(deg), "\n")

  if (nrow(deg_mapped) == 0) {
    stop("No genes could be mapped to STRING IDs.")
  }

  # ── 3. Get interactions ───────────────────────────────────────────────────
  interactions <- string_db$get_interactions(deg_mapped$STRING_id)
  cat("Interactions found:", nrow(interactions), "\n")

  if (nrow(interactions) == 0) {
    stop("No interactions found. Try lowering score_threshold to 200.")
  }

  # ── 4. Build node colour and size from logFC ──────────────────────────────
  deg_mapped$color <- ifelse(deg_mapped$logFC > 0, "#B71C1C", "#1A237E")
  deg_mapped$size  <- 1 + abs(deg_mapped$logFC) * 0.8

  # ── 5. Plot using STRINGdb's built-in network plot ────────────────────────
  string_db$plot_network(deg_mapped$STRING_id)
  graphics::title(main = title)

  # ── 6. Add legend manually ────────────────────────────────────────────────
  graphics::legend(
    "bottomright",
    legend = c("Up-regulated", "Down-regulated"),
    col    = c("#B71C1C", "#1A237E"),
    pch    = 16,
    pt.cex = 1.5,
    bty    = "n"
  )
  invisible(deg_mapped)
}




#' @title PPI Hub Network
#' @description Builds a protein-protein interaction network showing only the
#'   top N most connected hub proteins from DEGs. Less cluttered and more
#'   interpretable than the full network. Also prints a table of hub genes
#'   with their connection counts and logFC.
#' @param deg Data frame from load_deg_table (needs GeneSymbol and logFC columns)
#' @param score_threshold Minimum STRING score 0-1000 (default = 400)
#' @param top_n Number of top hub proteins to show (default = 30)
#' @param title Plot title
#' @param cache_dir Local folder to cache STRINGdb files so they are not
#'   re-downloaded each session
#' @return Invisibly returns a data frame of the top hub genes with
#'   their logFC and connection counts
#' @export
#' @examples
#' ppi_hub_network(deg)
#' ppi_hub_network(deg, top_n = 20, score_threshold = 700)
ppi_hub_network <- function(deg,
                            score_threshold = 400,
                            top_n           = 30,
                            title           = "PPI Hub Network of DEGs",
                            cache_dir       = file.path(tools::R_user_dir("collaborativeproject", which = "cache"), "stringdb")) {

  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("Please install STRINGdb first:\n  BiocManager::install('STRINGdb')")
  }

  # ── 1. Create cache directory if needed ───────────────────────────────────
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
    cat("Created STRINGdb cache at:", cache_dir, "\n")
  }

  # ── 2. Initialise STRINGdb ─────────────────────────────────────────────────
  cat("Initialising STRINGdb (loading from cache if available)...\n")
  string_db <- STRINGdb::STRINGdb$new(
    version         = "11.5",
    species         = 9606,
    score_threshold = score_threshold,
    input_directory = cache_dir
  )

  # ── 3. Map gene symbols to STRING IDs ─────────────────────────────────────
  deg_mapped <- string_db$map(
    my_data_frame        = deg,
    my_data_frame_id_col = "GeneSymbol",
    removeUnmappedRows   = TRUE
  )
  cat("Genes mapped:", nrow(deg_mapped), "of", nrow(deg), "\n")

  if (nrow(deg_mapped) == 0) {
    stop("No genes could be mapped to STRING IDs.")
  }

  # ── 4. Get all interactions among mapped genes ────────────────────────────
  interactions <- string_db$get_interactions(deg_mapped$STRING_id)
  cat("Total interactions found:", nrow(interactions), "\n")

  if (nrow(interactions) == 0) {
    stop("No interactions found. Try lowering score_threshold.")
  }

  # ── 5. Count connections per protein to identify hubs ─────────────────────
  connection_counts <- table(c(interactions$from, interactions$to))
  connection_df <- data.frame(
    STRING_id   = names(connection_counts),
    connections = as.integer(connection_counts),
    stringsAsFactors = FALSE
  )

  # ── 6. Keep only top N hubs ────────────────────────────────────────────────
  top_nodes <- connection_df[order(connection_df$connections, decreasing = TRUE), ]
  top_nodes <- top_nodes[1:min(top_n, nrow(top_nodes)), ]

  # ── 7. Filter mapped genes to top nodes only ──────────────────────────────
  deg_top <- deg_mapped[deg_mapped$STRING_id %in% top_nodes$STRING_id, ]
  deg_top <- merge(deg_top, top_nodes, by = "STRING_id")
  deg_top <- deg_top[order(deg_top$connections, decreasing = TRUE), ]

  # ── 8. Print hub gene table ────────────────────────────────────────────────
  cat("\nTop", nrow(deg_top), "hub genes:\n")
  print(deg_top[, c("GeneSymbol", "logFC", "padj", "connections")])

  # ── 9. Plot subnetwork ─────────────────────────────────────────────────────
  string_db$plot_network(deg_top$STRING_id)
  graphics::title(main = title)

  invisible(deg_top)
}
