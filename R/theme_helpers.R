COL <- list(
  ink        = "#243041",
  muted      = "#4B5563",
  grid       = "#D9E6DB",
  bg         = "white",
  tint       = "#EEF6EE",
  tint_dark  = "#E4F0E4",
  point      = "#5E8A6A",
  fit        = "#2F5D46",
  zero       = "#CADBCB",
  ribbon     = "#5E8A6A",
  raw        = "#6B7280",
  pos        = "#C23B3B",
  neg        = "#2D63C8",
  border     = "#4B5563",
  qq_band    = "#CFE2D1",
  qq_line    = "#2F5D46",
  loess      = "#6B7280",
  overall    = "#8B5CF6",
  low        = "#0F766E",
  mid        = "#374151",
  high       = "#B45309"
)

theme_modern <- function() {

  ggplot2::theme_minimal(base_size = 12) +

    ggplot2::theme(

      plot.background =
        ggplot2::element_rect(fill = COL$bg, color = NA),

      panel.background =
        ggplot2::element_rect(fill = COL$tint, color = NA),

      panel.grid.minor =
        ggplot2::element_blank(),

      panel.grid.major =
        ggplot2::element_line(
          color = COL$grid,
          linewidth = 0.35
        ),

      axis.title =
        ggplot2::element_text(
          color = COL$ink,
          size = 11
        ),

      axis.text =
        ggplot2::element_text(
          color = COL$muted
        ),

      plot.title =
        ggplot2::element_text(
          color = COL$ink,
          face = "bold",
          size = 14
        ),

      plot.subtitle =
        ggplot2::element_text(
          color = COL$muted,
          size = 10
        ),

      plot.caption =
        ggplot2::element_text(
          color = COL$muted,
          size = 9,
          hjust = 0
        ),

      plot.caption.position = "plot",

      strip.text =
        ggplot2::element_text(
          color = COL$ink,
          face = "bold"
        ),

      legend.title =
        ggplot2::element_text(
          color = COL$ink
        ),

      legend.text =
        ggplot2::element_text(
          color = COL$muted
        )
    )
}
