# ============================================================================
# Plotting helpers
# ----------------------------------------------------------------------------
# Tiny base-R plotting helpers shared across modules. Base R was chosen so
# the app shell stays dependency-free; modules can swap in ggplot2 / plotly
# later without changing this interface.
# ============================================================================

#' Scatter plot of an embedding, colored by a categorical metadata vector.
#'
#' @param emb     data.frame with `x`, `y` columns (from `get_embedding()`)
#' @param values  vector of length nrow(emb); coerced to factor for coloring
#' @param title   plot title
#' @param xlab,ylab  axis labels (defaults to "Dim 1" / "Dim 2")
plot_embedding_categorical <- function(emb, values, title = "",
                                       xlab = "Dim 1", ylab = "Dim 2") {
  stopifnot(!is.null(emb), nrow(emb) == length(values))
  f <- as.factor(values)
  pal <- grDevices::hcl.colors(max(length(levels(f)), 2L), palette = "Dark 3")
  graphics::par(mar = c(4, 4, 3, 8), xpd = NA)
  graphics::plot(emb$x, emb$y,
                 pch = 20, cex = 0.55,
                 col = pal[as.integer(f)],
                 xlab = xlab, ylab = ylab, main = title,
                 panel.first = graphics::grid(col = "gray90", lty = 1))
  graphics::legend("topright",
                   legend = levels(f),
                   col = pal[seq_along(levels(f))],
                   pch = 20, bty = "n",
                   cex = 0.8,
                   inset = c(-0.25, 0))
}

#' Scatter plot of an embedding, colored by a numeric expression vector
#' (FeaturePlot style).
plot_embedding_continuous <- function(emb, values, title = "",
                                      xlab = "Dim 1", ylab = "Dim 2",
                                      legend_title = "expr") {
  stopifnot(!is.null(emb), nrow(emb) == length(values))
  rng <- range(values, finite = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) rng <- rng + c(-0.5, 0.5)
  n_bins <- 64L
  pal <- grDevices::hcl.colors(n_bins, palette = "Viridis")
  bins <- pmax(1L, pmin(n_bins,
                        as.integer(cut(values, breaks = n_bins,
                                       include.lowest = TRUE))))
  # Draw low-expression cells first so high-expression cells sit on top.
  ord <- order(values, na.last = FALSE)
  graphics::par(mar = c(4, 4, 3, 8), xpd = NA)
  graphics::plot(emb$x[ord], emb$y[ord],
                 pch = 20, cex = 0.55,
                 col = pal[bins[ord]],
                 xlab = xlab, ylab = ylab, main = title,
                 panel.first = graphics::grid(col = "gray90", lty = 1))
  legend_color_strip(rng, pal, legend_title)
}

#' Volcano plot: -log10(p_val_adj) vs avg_log2FC.
#'
#' `de_df` columns required: `gene`, `avg_log2FC`, `p_val_adj`. Points are
#' coloured by significance under (`log2fc_thr`, `padj_thr`) -- red = up in
#' group_1, blue = up in group_2, grey = ns. Threshold guide lines are drawn
#' as dashed grey. The top `label_top` significant genes are text-labelled.
plot_volcano <- function(de_df, log2fc_thr = 0.5, padj_thr = 0.05,
                         label_top = 5) {
  if (is.null(de_df) || nrow(de_df) == 0L) {
    graphics::par(mar = c(4, 4, 3, 2))
    graphics::plot.new()
    graphics::title(main = "Volcano (no results)")
    return(invisible(NULL))
  }
  x <- de_df$avg_log2FC
  y <- -log10(pmax(de_df$p_val_adj, .Machine$double.eps))
  sig <- !is.na(de_df$p_val_adj) &
         de_df$p_val_adj <= padj_thr &
         abs(x) >= log2fc_thr
  col <- ifelse(sig, ifelse(x > 0, "#c62828", "#1565c0"), "#bbbbbb")
  graphics::par(mar = c(4, 4, 3, 2))
  graphics::plot(x, y, pch = 20, cex = 0.8, col = col,
                 xlab = "avg log2 fold-change  (group_1 vs group_2)",
                 ylab = "-log10(adj. p-value)",
                 main = "Volcano",
                 panel.first = graphics::grid(col = "gray90", lty = 1))
  graphics::abline(v = c(-log2fc_thr, log2fc_thr),
                   h = -log10(max(padj_thr, .Machine$double.eps)),
                   lty = 2, col = "gray60")
  if (label_top > 0L) {
    top_idx <- which(sig)
    if (length(top_idx) > label_top) {
      o <- order(de_df$p_val_adj[top_idx])[seq_len(label_top)]
      top_idx <- top_idx[o]
    }
    if (length(top_idx) > 0L) {
      graphics::text(x[top_idx], y[top_idx], labels = de_df$gene[top_idx],
                     pos = 3, cex = 0.75, offset = 0.4)
    }
  }
}

#' Horizontal bar plot of top enriched pathways.
#'
#' @param df       enrichment data.frame (output of `compute_enrichment()`)
#' @param top_n    integer; number of bars to show
#' @param metric   "padj" (default; uses -log10(p_val_adj)) or "odds_ratio"
#'
#' The plot uses one row per bar and integer y-positions, which makes it
#' easy for the module to map a click event to a row via `click$y`.
plot_pathway_enrichment <- function(df, top_n = 10L, metric = c("padj", "odds_ratio")) {
  metric <- match.arg(metric)
  if (is.null(df) || nrow(df) == 0L) {
    graphics::par(mar = c(4, 4, 3, 2)); graphics::plot.new()
    graphics::title(main = "Pathway enrichment (no results)")
    return(invisible(NULL))
  }
  df <- df[order(df$p_val_adj, df$p_val), , drop = FALSE]
  df <- utils::head(df, top_n)
  # Bar from bottom to top => most significant on top after we reverse.
  df <- df[rev(seq_len(nrow(df))), , drop = FALSE]
  x <- switch(metric,
              "padj"       = -log10(pmax(df$p_val_adj, .Machine$double.eps)),
              "odds_ratio" = ifelse(is.finite(df$odds_ratio), df$odds_ratio, 0))
  xlab <- switch(metric,
                 "padj"       = "-log10(adj. p-value)",
                 "odds_ratio" = "odds ratio")
  dir_pal <- c(up_in_g1 = "#c62828", up_in_g2 = "#1565c0", both = "#6a4c93")
  col <- dir_pal[as.character(df$direction)]
  col[is.na(col)] <- "#999999"
  graphics::par(mar = c(4, 18, 3, 2))
  bp <- graphics::barplot(x, names.arg = df$pathway, horiz = TRUE,
                          col = col, border = NA,
                          xlab = xlab, las = 1, cex.names = 0.85,
                          main = sprintf("Top %d pathways", nrow(df)))
  graphics::text(x = x, y = bp,
                 labels = sprintf("n=%d", df$n_overlap),
                 pos = 4, cex = 0.75, col = "#444444", offset = 0.2, xpd = NA)
  invisible(list(y_centers = as.numeric(bp), pathways = df$pathway))
}

#' Scatter of gene expression vs. pseudotime, with a binned trend line.
#'
#' @param pt           numeric pseudotime vector (any range; not rescaled)
#' @param expr         expression vector aligned with `pt`
#' @param gene_name    label for the y axis and title
#' @param n_bins       bins for the trend line (uses `bin_gene_by_pseudotime`)
#' @param point_color  point colour (raw scatter); the trend line is fixed
plot_gene_vs_pseudotime <- function(pt, expr, gene_name = "gene",
                                    n_bins = 20L,
                                    point_color = "#9aa0a6") {
  if (is.null(pt) || is.null(expr) ||
      length(pt) == 0L || length(pt) != length(expr)) {
    graphics::par(mar = c(4, 4, 3, 2)); graphics::plot.new()
    graphics::title(main = sprintf("%s vs pseudotime (no data)", gene_name))
    return(invisible(NULL))
  }
  ok <- !is.na(pt) & !is.na(expr)
  graphics::par(mar = c(4, 4, 3, 2))
  graphics::plot(pt[ok], expr[ok],
                 pch = 20, cex = 0.4,
                 col = grDevices::adjustcolor(point_color, alpha.f = 0.35),
                 xlab = "pseudotime", ylab = sprintf("%s expression", gene_name),
                 main = sprintf("%s vs pseudotime", gene_name),
                 panel.first = graphics::grid(col = "gray90", lty = 1))
  binned <- bin_gene_by_pseudotime(pt, expr, n_bins = n_bins)
  if (!is.null(binned)) {
    keep <- !is.na(binned$expr_mean)
    graphics::lines(binned$pt_mid[keep], binned$expr_mean[keep],
                    col = "#c62828", lwd = 2.5)
    graphics::points(binned$pt_mid[keep], binned$expr_mean[keep],
                     pch = 20, cex = 1.2, col = "#c62828")
  }
  invisible(binned)
}

#' Box + jittered-point plot of expression stratified by a group factor.
#'
#' Used as a lightweight stand-in for VlnPlot. No extra dependencies.
plot_expression_by_group <- function(expr, group, title = "",
                                     xlab = "", ylab = "expression") {
  stopifnot(length(expr) == length(group))
  keep <- !is.na(expr) & !is.na(group)
  expr <- expr[keep]; group <- group[keep]
  if (length(expr) == 0L) {
    graphics::plot.new(); graphics::title(main = paste(title, "(no data)"))
    return(invisible(NULL))
  }
  f <- as.factor(group)
  pal <- grDevices::hcl.colors(max(length(levels(f)), 2L), palette = "Dark 3")
  graphics::par(mar = c(5, 4, 3, 2))
  graphics::boxplot(expr ~ f, col = pal, outline = FALSE,
                    main = title, xlab = xlab, ylab = ylab,
                    las = if (length(levels(f)) > 4L) 2 else 1)
  for (i in seq_along(levels(f))) {
    sub <- expr[f == levels(f)[i]]
    if (length(sub) == 0L) next
    graphics::points(jitter(rep(i, length(sub)), amount = 0.15), sub,
                     pch = 20,
                     col = grDevices::adjustcolor(pal[i], alpha.f = 0.45),
                     cex = 0.5)
  }
}

# Tiny helper: draw a vertical color strip + min/max labels on the right
# margin. Used as a stand-in for a continuous color legend.
legend_color_strip <- function(rng, pal, label) {
  usr <- graphics::par("usr")
  x0 <- usr[2] + (usr[2] - usr[1]) * 0.02
  x1 <- usr[2] + (usr[2] - usr[1]) * 0.05
  y0 <- usr[3] + (usr[4] - usr[3]) * 0.1
  y1 <- usr[3] + (usr[4] - usr[3]) * 0.9
  n  <- length(pal)
  ys <- seq(y0, y1, length.out = n + 1L)
  for (i in seq_len(n)) {
    graphics::rect(x0, ys[i], x1, ys[i + 1L], col = pal[i], border = NA)
  }
  graphics::text(x1, y1, sprintf("%.2f", rng[2]), pos = 4, cex = 0.75)
  graphics::text(x1, y0, sprintf("%.2f", rng[1]), pos = 4, cex = 0.75)
  graphics::text((x0 + x1) / 2, y1 + (usr[4] - usr[3]) * 0.04,
                 label, cex = 0.8, font = 2)
}

#' Regulon heatmap: groups x regulons matrix of mean AUC.
#'
#' Base-R `image()` rendering. Rows = grouping levels (cluster ids),
#' columns = regulons. Color = mean AUC. Designed for the Regulons
#' module; could be reused for any score-by-group matrix.
#'
#' @param mat    numeric matrix with dimnames (rows = groups, cols = regulons)
#' @param title  plot title
#' @param low,high  endpoints for the viridis-like palette
plot_regulon_heatmap <- function(mat, title = "",
                                 low = "#440154", high = "#FDE725") {
  if (is.null(mat) || !is.matrix(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
    graphics::plot.new()
    graphics::title(main = title)
    graphics::text(0.5, 0.5, "(empty)", cex = 1.2, col = "gray60")
    return(invisible(NULL))
  }
  pal <- grDevices::colorRampPalette(c(low, "#3B528B", "#21908C",
                                       "#5DC863", high))(64L)
  rng <- range(mat, na.rm = TRUE, finite = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) rng <- rng + c(-0.5, 0.5)
  # image() draws cols-first; transpose so rows of `mat` (groups) end up
  # as rows of the heatmap visually.
  graphics::par(mar = c(7, 7, 3, 8), xpd = NA)
  graphics::image(seq_len(ncol(mat)), seq_len(nrow(mat)),
                  t(mat[nrow(mat):1L, , drop = FALSE]),
                  col = pal, xlab = "", ylab = "", axes = FALSE,
                  main = title, useRaster = TRUE,
                  zlim = rng)
  graphics::axis(1, at = seq_len(ncol(mat)),
                 labels = colnames(mat), las = 2, cex.axis = 0.85)
  graphics::axis(2, at = seq_len(nrow(mat)),
                 labels = rev(rownames(mat)), las = 2, cex.axis = 0.85)
  graphics::box()
  legend_color_strip(rng, pal, "AUC")
}
