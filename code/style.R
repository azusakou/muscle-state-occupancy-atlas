palette_contract <- c(
  neutral_dark = "#272727",
  neutral_mid = "#767676",
  neutral_light = "#D8D8D8",
  baseline = "#4D4D4D",
  exercise = "#2C7FB8",
  overload = "#1B9E77",
  denervation = "#D95F02",
  injury = "#7570B3",
  disuse = "#A6761D",
  recovery = "#66A61E",
  gold = "#E6AB02"
)

theme_nature_contract <- function(base_size = 9.5, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.5, colour = "black"),
      legend.title = ggplot2::element_text(size = base_size - 0.3),
      legend.text = ggplot2::element_text(size = base_size - 0.7),
      strip.text = ggplot2::element_text(size = base_size - 0.3, face = "bold"),
      plot.title = ggplot2::element_text(size = base_size + 0.5, face = "bold"),
      plot.tag = ggplot2::element_text(size = 10, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}
