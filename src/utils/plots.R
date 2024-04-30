## settings

mesh_generation_plot_settings <- function() {
  standard_plot_settings_fiends <- theme_minimal() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = 14,
        hjust = 0.5,
        vjust = 1
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(angle = 45, hjust = 1)
    )
}