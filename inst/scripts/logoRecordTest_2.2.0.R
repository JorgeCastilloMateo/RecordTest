# SCRIPT TO REPRODUCE RecordTest LOGO (R VERSION 4.2.2)
url <- "http://cran.r-project.org/src/contrib/Archive/RecordTest/RecordTest_1.0.1.tar.gz"
install.packages(url, repos = NULL, type = "source")

set.seed(23)
p <- RecordTest::records(
  stats::rnorm(23), 
  colour = c(
    "black", 
    grDevices::rgb(0, 176, 240, maxColorValue = 255), 
    grDevices::rgb(0, 176, 80, maxColorValue = 255), 
    "red"),
  alpha = c(1, 1, 1, 1))

p <- p + 
  ggplot2::labs(
    title    = ggplot2::element_blank(), 
    subtitle = ggplot2::element_blank()) +
  ggplot2::theme(
    axis.line        = ggplot2::element_blank(), 
    axis.text.x      = ggplot2::element_blank(),
    axis.text.y      = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.title.x     = ggplot2::element_blank(),
    axis.title.y     = ggplot2::element_blank(),
    legend.position  = "none",
    panel.background = ggplot2::element_blank(),
    panel.border     = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.background  = ggplot2::element_blank())

if (!require("hexSticker")) install.packages("hexSticker")
hexSticker::sticker(p, package = "RecordTest", 
  s_x = 1, s_y = 1.25, s_width = 1.7, s_height = 1.2,
  p_x = 1, p_y = 0.6, p_size = 23, h_size = 1.5,
  h_fill = grDevices::rgb(229, 184, 183, maxColorValue = 255), 
  h_color = grDevices::rgb(148, 54, 52, maxColorValue = 255),
  filename = "/inst/img/logoRecordTest_2.2.0.png")
