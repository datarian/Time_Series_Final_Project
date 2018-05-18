library(ggplot2)

ggplot2::theme_set(ggplot2::theme_bw(base_size = 11))
ggplot2::theme_update(
    plot.background = ggplot2::element_rect(colour="white"),
    plot.title = element_text(hjust = 0.5),
    panel.background = ggplot2::element_rect(colour = "black", fill=NA),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(colour = "black", size=0.5),
    axis.line = ggplot2::element_line(colour = "black"),
    axis.ticks = ggplot2::element_line(colour = "black"),
    axis.title.x = ggplot2::element_text(face="bold", margin = ggplot2::margin(10,0,0,0)),
    axis.title.y = ggplot2::element_text(face="bold", margin = ggplot2::margin(0,10,0,0)),
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(3,0,0,0), colour="black"),
    axis.text.y = ggplot2::element_text(margin = ggplot2::margin(0,3,0,0), colour = "black"),
    strip.text.x = ggplot2::element_text(face="bold"),
    strip.background = ggplot2::element_rect(colour="black", fill="white"),
    legend.key = ggplot2::element_rect(colour="white", size = 0.3),
    legend.title = ggplot2::element_blank(),
    legend.position = "bottom",
    aspect.ratio = 0.618
)

ggplot2::scale_fill_manual(values=ggthemes::scale_fill_colorblind(8))
ggplot2::scale_color_manual(values=ggthemes::scale_colour_colorblind(8))
pal <-ggthemes::colorblind_pal()(8)
