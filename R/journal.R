library(ggplot2)
library(scales)

theme_journal <- function(...){
  theme_bw() +
    theme( legend.position = "inside", panel.grid = element_blank(),
         , legend.justification.inside = c(1,1)
         , legend.background = element_rect(fill = NA)
         , ...
    )
}

save_journal <- function(filename, plot=last_plot(),
                         width = 7, height = 5, scale = 0.6,
                         ...
){
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    scale = scale,
    ...
  )
}

x <- g3 |> ggplot(aes(weight, linetype=n))

geom_density_journal <- function(...){
    # geom_freqpoly(binwidth = 0.05) +
  list(
    stat_density(
      aes(y = after_stat(count)),
      geom="line",
      position="identity"
      , alpha = 0.85
      , ...
    ),
    scale_y_continuous(expand = expansion(mult = c(0,.1)))
  )
}
