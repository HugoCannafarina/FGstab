#' Plot stability paths from an ipss_result object
#'
#' Displays the selection frequency of each feature along the regularisation
#' path. Features that are ultimately selected are highlighted in colour;
#' all others are drawn in light grey. The top-k most stable features
#' (by maximum selection frequency) are labelled on the right margin.
#'
#' @param x An \code{\link{ipss_result}} object returned by
#'   \code{\link{FGstab}}.
#' @param n_highlight Integer. Number of top features to label on the right
#'   margin, ranked by maximum selection frequency. Default \code{5}.
#' @param selected_color Character. Colour used for selected feature paths.
#'   Default \code{"#E63946"} (red).
#' @param alpha_unselected Numeric in \code{[0, 1]}. Transparency of
#'   unselected feature paths. Default \code{0.4}.
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{ggplot} object (invisible). The plot is also printed.
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' res <- FGstab(X, y, n_alphas = 15, target_fdr = 0.2)
#' plot(res)
#' plot(res, n_highlight = 10)
#' }
#'
#' @export
plot.ipss_result <- function(x,
                             n_highlight      = 5L,
                             selected_color   = "#E63946",
                             alpha_unselected = 0.4,
                             ...) {

  paths    <- x$stability_paths   # matrix (n_alphas x p)
  n_alphas <- nrow(paths)
  p        <- ncol(paths)
  steps    <- seq_len(n_alphas)

  # Feature names
  feat_names <- if (!is.null(colnames(paths))) {
    colnames(paths)
  } else {
    paste0("V", seq_len(p))
  }

  # Max selection frequency per feature (computed from stability_paths)
  max_freq <- apply(paths, 2, max, na.rm = TRUE)

  # Selected and top-k features
  selected_idx <- x$selected_features
  top_idx      <- order(max_freq, decreasing = TRUE)[seq_len(min(n_highlight, p))]

  # Long data frame
  df <- data.frame(
    step     = rep(steps, times = p),
    prob     = as.vector(paths),
    feature  = rep(feat_names, each = n_alphas),
    feat_idx = rep(seq_len(p), each = n_alphas)
  )
  df$is_selected <- df$feat_idx %in% selected_idx
  df$is_top      <- df$feat_idx %in% top_idx

  # Labels at the last step for top features
  label_df <- df[df$is_top & df$step == n_alphas, ]

  subtitle <- if (length(selected_idx) > 0) {
    paste0(length(selected_idx), " feature(s) selected")
  } else {
    "No features selected (set target_fp or target_fdr)"
  }

  p_out <- ggplot(df, aes(x = step, y = prob, group = feature)) +

    geom_line(
      data  = df[!df$is_selected, ],
      color = "grey80",
      alpha = alpha_unselected,
      linewidth  = 0.35
    ) +

    geom_line(
      data  = df[df$is_selected, ],
      color = selected_color,
      linewidth  = 0.8
    ) +

    geom_text(
      data      = label_df,
      aes(label = feature),
      hjust     = -0.1,
      size      = 3,
      color     = ifelse(label_df$is_selected, selected_color, "grey50")
    ) +

    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +

    labs(
      title    = "Stability paths \u2014 IPSS",
      subtitle = subtitle,
      x        = "Regularisation step",
      y        = "Selection frequency"
    ) +

    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(color = "grey40")
    )

  print(p_out)
  invisible(p_out)
}
