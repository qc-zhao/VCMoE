.coef_data_frame <- function(object, coefficient_set = c("expert", "gating")) {
  coefficient_set <- match.arg(coefficient_set)
  values <- object$coefficients[[coefficient_set]]
  grid <- expand.grid(
    u = object$u_grid,
    component = dimnames(values)[[2L]],
    term = dimnames(values)[[3L]],
    KEEP.OUT.ATTRS = FALSE
  )
  grid$estimate <- as.vector(values)
  grid$coefficient_set <- coefficient_set
  grid
}

#' Plot fitted coefficient functions
#'
#' @param object A `vcmoe` object.
#' @param type `"expert"` or `"gating"`.
#' @return A `ggplot` object.
#' @export
plot_coefficients <- function(object, type = c("expert", "gating")) {
  type <- match.arg(type)
  df <- .coef_data_frame(object, type)
  ggplot2::ggplot(df, ggplot2::aes(x = u, y = estimate, color = component)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~term, scales = "free_y") +
    ggplot2::labs(x = "u", y = "coefficient", color = "component") +
    ggplot2::theme_minimal(base_size = 12)
}

#' Plot fitted posterior summaries over the coefficient grid
#'
#' @param object A `vcmoe` object.
#' @return A `ggplot` object.
#' @export
plot_posterior <- function(object) {
  posterior_mean <- t(apply(object$posterior, 3L, colMeans))
  df <- data.frame(
    u = rep(object$u_grid, times = object$k),
    component = rep(colnames(posterior_mean), each = length(object$u_grid)),
    posterior_mean = as.vector(posterior_mean)
  )
  ggplot2::ggplot(df, ggplot2::aes(x = u, y = posterior_mean, color = component)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(x = "u", y = "mean posterior probability", color = "component") +
    ggplot2::theme_minimal(base_size = 12)
}
