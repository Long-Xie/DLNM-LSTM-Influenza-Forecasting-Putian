###############figure3
library(ggplot2)
library(scales)      # for rescale
library(patchwork)   # for multi-panel layout

#====================
# Data preprocessing
#====================
wtcases$dirrus <- wtcases$winddir
wtcases$dirrus[which(wtcases$dirrus > 337.5)] <- 
  wtcases$winddir[which(wtcases$dirrus > 337.5)] - 360
wtcases$Date <- as.Date(wtcases$Date)

#====================
# Shared plotting function
#====================
plot_wt <- function(y1, y1min, y1max, brea){
  ggplot(wtcases, aes(x = Date)) +
    geom_line(aes(y = y1, color = "y1"), linewidth = line_size, show.legend = TRUE) + 
    geom_line(aes(y = rescale(typeA, c(y1min, y1max)), color = "typeA"), 
              linewidth = line_size, show.legend = FALSE) +
    geom_line(aes(y = rescale(typeB, c(y1min, y1max)), color = "typeB"), 
              linewidth = line_size, show.legend = FALSE) +
    scale_color_manual(
      name = NULL,
      values = c("y1"="#F7EBC6", "typeA"="#F07673", "typeB"="#7DAEE0"),
      breaks = c("y1","typeA","typeB"),
      labels = c("Weather indicators",
                 "Cases of influenza A virus",
                 "Cases of influenza B virus")
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) + 
    scale_x_date(
      breaks = c(as.Date("2018-06-01"), as.Date("2019-06-01"),
                 as.Date("2020-06-01"), as.Date("2021-06-01"),
                 as.Date("2022-06-01"), as.Date("2023-06-01")),
      date_labels = "%Y",
      expand = c(0,0),
      limits = c(as.Date("2018-01-01"), as.Date("2023-12-31"))
    ) + 
    scale_y_continuous(
      expand = c(0,0),
      limits = c(y1min, y1max),
      breaks = brea,
      sec.axis = sec_axis(~rescale(., c(0,60)), name = "Number of cases")
    ) +
    theme_bw() + 
    labs(x = "Date") +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(), 
      axis.line = element_line(color = "black", linewidth = line_size),
      axis.title = element_text(color = "black", size = text_size, family = family), 
      axis.text.x = element_text(color = "black", size = text_size, family = family), 
      axis.text.y = element_text(color = "black", angle = 0, size = text_size, family = family), 
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(color = "black", size = text_size, family = family),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt")
    )
}

#====================
# Subplots
#====================
p1 <- plot_wt(wtcases[,"temp"], 8, 33, c(8,10,15,20,25,30,33)) + 
  labs(y = "Average temperature (°C)") 
p2 <- plot_wt(wtcases[,"humidity"], 30, 100, c(30,50,70,90,100)) + 
  labs(y = "Humidity (%)")
p3 <- plot_wt(wtcases[,"precip"], 0, 55, c(0,10,20,30,40,50)) + 
  labs(y = "Precipitation (mm)")
p4 <- plot_wt(wtcases[,"solarradiation"], 10, 340, c(50,100,150,200,250,300)) + 
  labs(y = expression("Solar radiation (W/m"^2*")"))
p5 <- plot_wt(wtcases[,"tempdiff"], 1, 16, c(1,4,7,10,13,16)) + 
  labs(y = "Daily temperature ranges (°C)")
p6 <- plot_wt(wtcases[,"uvindex"], 0, 10, c(0,2,4,6,8,10)) + 
  labs(y = "Ultraviolet index")
p7 <- plot_wt(wtcases[,"windspeed"], 8, 32, c(8,12,16,20,24,28,32)) + 
  labs(y = "Wind speed (km/h)")

# Wind direction panel
p8 <- ggplot(wtcases, aes(x = Date)) +
  geom_line(aes(y = wtcases[,"dirrus"], color = "y1"), linewidth = line_size, show.legend = TRUE) + 
  geom_line(aes(y = rescale(typeA, c(-22.5,337.5)), color = "typeA"), 
            linewidth = line_size, show.legend = FALSE) +
  geom_line(aes(y = rescale(typeB, c(-22.5,337.5)), color = "typeB"), 
            linewidth = line_size, show.legend = FALSE) +
  scale_color_manual(
    name = NULL,
    values = c("y1"="#F7EBC6", "typeA"="#F07673", "typeB"="#7DAEE0"),
    breaks = c("y1","typeA","typeB"),
    labels = c("Weather indicators",
               "Cases of influenza A virus",
               "Cases of influenza B virus")
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) + 
  scale_x_date(
    breaks = c(as.Date("2018-06-01"), as.Date("2019-06-01"),
               as.Date("2020-06-01"), as.Date("2021-06-01"),
               as.Date("2022-06-01"), as.Date("2023-06-01")),
    date_labels = "%Y",
    expand = c(0,0),
    limits = c(as.Date("2018-01-01"), as.Date("2023-12-31"))
  ) + 
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(-22.5,337.5),
    breaks = c(0,45,90,135,180,225,270,315),
    labels = c("N","NE","E","SE","S","SW","W","NW"),
    sec.axis = sec_axis(~rescale(.,c(0,60)), name = "Number of cases")
  ) +
  theme_bw() + 
  labs(x = "Date", y = "Wind direction") +
  theme(
    panel.border = element_blank(), 
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black", linewidth = line_size), 
    axis.title = element_text(color = "black", size = text_size, family = family), 
    axis.text.x = element_text(color = "black", size = text_size, family = family), 
    axis.text.y = element_text(color = "black", angle = 0, size = text_size, family = family),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(color = "black", size = text_size, family = family),
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt")
  )

#====================
# Combine 9×1 layout
#====================
final_plot <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + guide_area() +
  plot_layout(ncol = 1, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(panel.spacing = unit(5, "mm"),
        legend.position = "bottom",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        plot.tag = element_text(size = text_size, colour = "black", family = family),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"))

ggsave(file.path(output_dir, "Fig3_meteorology_timeseries_9x1.tiff"),
       final_plot, height = 18, width = 12, units = "in", bg = NULL, dpi = dpi)





