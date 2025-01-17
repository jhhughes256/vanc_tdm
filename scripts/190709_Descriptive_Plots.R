# Experimentation with different descriptive plots
# -----------------------------------------------------------------------------
# For presentation of this work various plots are required. Results include:
#  - Concentration vs. Time
#  - Percentage of patients above target AUC 
#  - Average dose over time
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Prepare work environment
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
  # setwd("E:/Hughes/Git/vanc_tdm")

# Load package libraries
  library(magrittr)
  library(ggplot2)  # Graphical package

# Source functions utility
  source("scripts/functions_utility.R")

# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Set colourblind palette
  cbPalette <- data.frame(
		grey = "#999999",
		orange = "#E69F00",
		skyblue = "#56B4E9",
		green = "#009E73",
		yellow = "#F0E442",
		blue = "#0072B2",
		red = "#D55E00",
		pink = "#CC79A7",
		stringsAsFactors = F
	)

# Read in data and convert to same format
  sagov_tb <- readr::read_rds("output/vanc_regimen_sagov.rds") %>%
    dplyr::select(ID, bayes) %>% tidyr::unnest() %>%
    tibble::add_column(method = "sagov")
  
  auc24_tb <- readr::read_rds("output/vanc_regimen_auc24.rds") %>%
    dplyr::select(ID, bayes) %>% tidyr::unnest() %>%
    tibble::add_column(method = "auc24")
  
  bayes_tb <- readr::read_rds("output/vanc_regimen_bayes_fix.rds") %>%
    dplyr::select(ID, bayes) %>% tidyr::unnest() %>%
    tibble::add_column(method = "bayes")
  
  baynl_tb <- readr::read_rds("output/vanc_regimen_bayes_newload.rds") %>%
    dplyr::select(ID, bayes) %>% tidyr::unnest() %>%
    tibble::add_column(method = "baynl")

# Bind into plot dataset and calculate cumulative dose
  plot_tb <- dplyr::bind_rows(sagov_tb, auc24_tb, bayes_tb, baynl_tb) %>%
    dplyr::mutate(methodf = factor(method)) %>%
    dplyr::mutate(CRCLf = factor(CRCL))
  levels(plot_tb$methodf) <- c(
    "Proportional TDM",
    "Model-Based TDM",
    "Model-Based w/ Loading Dose",
    "SA Health Guidelines"
  )  # set levels
  levels(plot_tb$CRCLf) <- c(
    "GFR - 25 mL/min",
    "GFR - 50 mL/min",
    "GFR - 75 mL/min",
    "GFR - 100 mL/min"
  )
  
# Create AUC dataset
  auc_tb <- dplyr::filter(plot_tb, time %in% (0:7*24)) %>%
    dplyr::group_by(ID) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(data) {
      dplyr::mutate(data, dAUC = c(0, diff(AUC)))
    })) %>% tidyr::unnest() %>%
    dplyr::filter(time != 0) %>%
    dplyr::mutate(time = time/24)
  
  auc_target <- auc_tb %>%
    dplyr::group_by(time, method, methodf) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(data) {
      tibble::tibble(
        gt400 = sum(data$dAUC > 400)/dim(data)[1],
        gt700 = sum(data$dAUC > 700)/dim(data)[1],
        gt400lt700 = sum(data$dAUC >= 400 & data$dAUC <= 700)/dim(data)[1]
      )
    })) %>% tidyr::unnest()
  
  auc_target_bycrcl <- auc_tb %>%
    dplyr::group_by(time, method, methodf, CRCLf) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(data) {
      tibble::tibble(
        gt400 = sum(data$dAUC > 400)/dim(data)[1],
        gt700 = sum(data$dAUC > 700)/dim(data)[1],
        gt400lt700 = sum(data$dAUC >= 400 & data$dAUC <= 700)/dim(data)[1]
      )
    })) %>% tidyr::unnest()
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Evaluation Plots
  palette <- with(cbPalette, c(red, blue, pink, green))

# Trough Concentration-time plots
  p <- NULL
  p <- ggplot(data = dplyr::filter(plot_tb, method == "sagov"))
  p <- p + ggtitle("Concentration Time Profile - SA Health Guidelines")
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "line", fun.y = median,
    colour = cbPalette$green, size = 1)
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "ribbon",
    fun.ymin = CI90lo,  fun.ymax = CI90hi, fill = cbPalette$green, 
    size = 1, alpha = 0.25)
  p <- p + scale_x_continuous("\nTime (hours)", breaks = 0:7*24)
  p <- p + scale_y_continuous("Vancomycin Concentration (mg/L)\n")
  # p <- p + coord_cartesian(xlim = c(0, 168), ylim = NULL)
  p1_sagov <- p + facet_wrap(~CRCLf)
  p1_sagov
  
  p <- NULL
  p <- ggplot(data = dplyr::filter(plot_tb, method == "auc24"))
  p <- p + ggtitle("Concentration Time Profile - Proportional TDM")
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "line", fun.y = median,
    colour = cbPalette$red, size = 1)
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "ribbon",
    fun.ymin = CI90lo,  fun.ymax = CI90hi, fill = cbPalette$red, 
    size = 1, alpha = 0.25)
  p <- p + scale_x_continuous("\nTime (hours)", breaks = 0:7*24)
  p <- p + scale_y_continuous("Vancomycin Concentration (mg/L)\n")
  # p <- p + coord_cartesian(xlim = c(0, 168), ylim = NULL)
  p1_auc24 <- p + facet_wrap(~CRCLf)
  p1_auc24
  
  p <- NULL
  p <- ggplot(data = dplyr::filter(plot_tb, method == "bayes"))
  p <- p + ggtitle("Concentration Time Profile - Model-Based TDM")
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "line", fun.y = median,
    colour = cbPalette$blue, size = 1)
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "ribbon",
    fun.ymin = CI90lo,  fun.ymax = CI90hi, fill = cbPalette$blue, 
    size = 1, alpha = 0.25)
  p <- p + scale_x_continuous("\nTime (hours)", breaks = 0:7*24)
  p <- p + scale_y_continuous("Vancomycin Concentration (mg/L)\n")
  # p <- p + coord_cartesian(xlim = c(0, 168), ylim = NULL)
  p1_bayes <- p + facet_wrap(~CRCLf)
  p1_bayes
  
  p <- NULL
  p <- ggplot(data = dplyr::filter(plot_tb, method == "baynl"))
  p <- p + ggtitle("Concentration Time Profile - Model-Based w/ Loading Dose")
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "line", fun.y = median,
    colour = cbPalette$pink, size = 1)
  p <- p + stat_summary(aes(x = time, y = IPRE), geom = "ribbon",
    fun.ymin = CI90lo,  fun.ymax = CI90hi, fill = cbPalette$pink, 
    size = 1, alpha = 0.25)
  p <- p + scale_x_continuous("\nTime (hours)", breaks = 0:7*24)
  p <- p + scale_y_continuous("Vancomycin Concentration (mg/L)\n")
  # p <- p + coord_cartesian(xlim = c(0, 168), ylim = NULL)
  p1_baynl <- p + facet_wrap(~CRCLf)
  p1_baynl
  
  p <- NULL
  p <- ggplot(data = plot_tb)
  p <- p + ggtitle("Concentration Time Profile - TDM Method Comparison")
  p <- p + stat_summary(aes(x = time, y = IPRE, colour = methodf), 
    geom = "line", fun.y = median, size = 0.5, alpha = 0.8)
  p <- p + scale_x_continuous("\nTime (hours)", breaks = 0:7*24)
  p <- p + scale_y_continuous("Vancomycin Concentration (mg/L)\n")
  p <- p + scale_colour_manual("TDM Method", values = palette)
  p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1)))
  # p <- p + coord_cartesian(xlim = c(0, 168), ylim = NULL)
  p1_all <- p + facet_wrap(~CRCLf)
  p1_all
  
# Plot % patients with AUC > Target
  p <- NULL
  p <- ggplot(data = dplyr::filter(auc_target_bycrcl, method == "sagov"))
  p <- p + ggtitle("Probability of Target Attainment - SA Health Guidelines")
  p <- p + geom_line(aes(x = time, y = gt400*100),
    size = 1, colour = cbPalette$green)
  p <- p + geom_line(aes(x = time, y = gt700*100),
    size = 1, colour = cbPalette$red)
  p <- p + geom_line(aes(x = time, y = gt400lt700*100),
    size = 1, colour = cbPalette$blue)
  p <- p + scale_x_continuous("\nTime (days)", breaks = 0:7)
  p <- p + scale_y_continuous("Probability of Target Attainment (%)\n")
  p <- p + coord_cartesian(xlim = NULL, ylim = c(0, 100))
  p2_sagov <- p + facet_wrap(~CRCLf)
  p2_sagov
  
  p <- NULL
  p <- ggplot(data = dplyr::filter(auc_target_bycrcl, method == "auc24"))
  p <- p + ggtitle("Probability of Target Attainment - Proportional TDM")
  p <- p + geom_line(aes(x = time, y = gt400*100),
    size = 1, colour = cbPalette$green)
  p <- p + geom_line(aes(x = time, y = gt700*100),
    size = 1, colour = cbPalette$red)
  p <- p + geom_line(aes(x = time, y = gt400lt700*100),
    size = 1, colour = cbPalette$blue)
  p <- p + scale_x_continuous("\nTime (days)", breaks = 0:7)
  p <- p + scale_y_continuous("Probability of Target Attainment (%)\n")
  p <- p + coord_cartesian(xlim = NULL, ylim = c(0, 100))
  p2_auc24 <- p + facet_wrap(~CRCLf)
  p2_auc24
  
  p <- NULL
  p <- ggplot(data = dplyr::filter(auc_target_bycrcl, method == "bayes"))
  p <- p + ggtitle("Probability of Target Attainment - Model-Based TDM")
  p <- p + geom_line(aes(x = time, y = gt400*100),
    size = 1, colour = cbPalette$green)
  p <- p + geom_line(aes(x = time, y = gt700*100),
    size = 1, colour = cbPalette$red)
  p <- p + geom_line(aes(x = time, y = gt400lt700*100),
    size = 1, colour = cbPalette$blue)
  p <- p + scale_x_continuous("\nTime (days)", breaks = 0:7)
  p <- p + scale_y_continuous("Probability of Target Attainment (%)\n")
  p <- p + coord_cartesian(xlim = NULL, ylim = c(0, 100))
  p2_bayes <- p + facet_wrap(~CRCLf)
  p2_bayes
  
  p <- NULL
  p <- ggplot(data = dplyr::filter(auc_target_bycrcl, method == "baynl"))
  p <- p + ggtitle("Probability of Target Attainment - Model-Based w/ Loading Dose")
  p <- p + geom_line(aes(x = time, y = gt400*100),
    size = 1, colour = cbPalette$green)
  p <- p + geom_line(aes(x = time, y = gt700*100),
    size = 1, colour = cbPalette$red)
  p <- p + geom_line(aes(x = time, y = gt400lt700*100),
    size = 1, colour = cbPalette$blue)
  p <- p + scale_x_continuous("\nTime (days)", breaks = 0:7)
  p <- p + scale_y_continuous("Probability of Target Attainment (%)\n")
  p <- p + coord_cartesian(xlim = NULL, ylim = c(0, 100))
  p2_baynl <- p + facet_wrap(~CRCLf)
  p2_baynl
  
# Plot % patients with AUC within Target Range
  p <- NULL
  p <- ggplot(data = dplyr::filter(auc_target, method == "sagov"))
  p <- p + ggtitle("Probability of Target Attainment - SA Health Guidelines")
  p <- p + geom_line(aes(x = time, y = gt400*100),
    size = 1, colour = cbPalette$green)
  p <- p + geom_line(aes(x = time, y = gt700*100),
    size = 1, colour = cbPalette$red)
  p <- p + geom_line(aes(x = time, y = gt400lt700*100),
    size = 1, colour = cbPalette$blue)
  p <- p + scale_x_continuous("\nTime (days)", breaks = 0:7)
  p <- p + scale_y_continuous("Probability of Target Attainment (%)\n")
  p3_sagov <- p + coord_cartesian(xlim = NULL, ylim = c(0, 100))
  p3_sagov
  
  p <- NULL
  p <- ggplot(data = dplyr::filter(auc_target, method == "auc24"))
  p <- p + ggtitle("Probability of Target Attainment - Proportional TDM")
  p <- p + geom_line(aes(x = time, y = gt400*100),
    size = 1, colour = cbPalette$green)
  p <- p + geom_line(aes(x = time, y = gt700*100),
    size = 1, colour = cbPalette$red)
  p <- p + geom_line(aes(x = time, y = gt400lt700*100),
    size = 1, colour = cbPalette$blue)
  p <- p + scale_x_continuous("\nTime (days)", breaks = 0:7)
  p <- p + scale_y_continuous("Probability of Target Attainment (%)\n")
  p3_auc24 <- p + coord_cartesian(xlim = NULL, ylim = c(0, 100))
  p3_auc24
  
  p <- NULL
  p <- ggplot(data = dplyr::filter(auc_target, method == "bayes"))
  p <- p + ggtitle("Probability of Target Attainment - Model-Based TDM")
  p <- p + geom_line(aes(x = time, y = gt400*100),
    size = 1, colour = cbPalette$green)
  p <- p + geom_line(aes(x = time, y = gt700*100),
    size = 1, colour = cbPalette$red)
  p <- p + geom_line(aes(x = time, y = gt400lt700*100),
    size = 1, colour = cbPalette$blue)
  p <- p + scale_x_continuous("\nTime (days)", breaks = 0:7)
  p <- p + scale_y_continuous("Probability of Target Attainment (%)\n")
  p3_bayes <- p + coord_cartesian(xlim = NULL, ylim = c(0, 100))
  p3_bayes
  
  p <- NULL
  p <- ggplot(data = dplyr::filter(auc_target, method == "baynl"))
  p <- p + ggtitle("Probability of Target Attainment - Model-Based w/ Loading Dose")
  p <- p + geom_line(aes(x = time, y = gt400*100),
    size = 1, colour = cbPalette$green)
  p <- p + geom_line(aes(x = time, y = gt700*100),
    size = 1, colour = cbPalette$red)
  p <- p + geom_line(aes(x = time, y = gt400lt700*100),
    size = 1, colour = cbPalette$blue)
  p <- p + scale_x_continuous("\nTime (days)", breaks = 0:7)
  p <- p + scale_y_continuous("Probability of Target Attainment (%)\n")
  p3_baynl <- p + coord_cartesian(xlim = NULL, ylim = c(0, 100))
  p3_baynl
  
  auc_target$methodf <- relevel(auc_target$methodf, "Model-Based TDM")
  auc_target$methodf <- relevel(auc_target$methodf, "Proportional TDM")
  auc_target$methodf <- relevel(auc_target$methodf, "SA Health Guidelines")
  
  p <- NULL
  p <- ggplot(data = auc_target)
  p <- p + ggtitle("Probability of Target Attainment - TDM Method Comparison")
  p <- p + geom_line(aes(x = time, y = gt400*100),
    size = 1, colour = cbPalette$green)
  p <- p + geom_line(aes(x = time, y = gt700*100),
    size = 1, colour = cbPalette$red)
  p <- p + geom_line(aes(x = time, y = gt400lt700*100),
    size = 1, colour = cbPalette$blue)
  p <- p + scale_x_continuous("\nTime (days)", breaks = 0:7)
  p <- p + scale_y_continuous("Probability of Target Attainment (%)\n")
  p <- p + coord_cartesian(xlim = NULL, ylim = c(0, 100))
  p3_all <- p + facet_wrap(~ methodf, nrow = 1)
  p3_all_square <- p + facet_wrap(~ methodf, nrow = 2)
  p3_all
  p3_all_square
  
  # levels(auc_target$methodf)[4] <- "Model-Based w/ Load"
  
  p1_sagov
  ggsave("output/CvT_sagov.png", width = 8, height = 6)
  p1_auc24
  ggsave("output/CvT_auc24.png", width = 8, height = 6)
  p1_bayes
  ggsave("output/CvT_bayes.png", width = 8, height = 6)
  # ggsave("output/CvT_bayes_addsamp.png", width = 8, height = 6)
  p1_baynl
  ggsave("output/CvT_bayes_load.png", width = 8, height = 6)
  p1_all
  ggsave("output/CvT_all.png", width = 8, height = 6)
  # ggsave("output/CvT_all_addsamp.png", width = 8, height = 6)
  p2_sagov
  ggsave("output/AUC_bycrcl_sagov.png", width = 8, height = 6)
  p2_auc24
  ggsave("output/AUC_bycrcl_auc24.png", width = 8, height = 6)
  p2_bayes
  ggsave("output/AUC_bycrcl_bayes.png", width = 8, height = 6)
  # ggsave("output/AUC_bycrcl_bayes_addsamp.png", width = 8, height = 6)
  p2_baynl
  ggsave("output/AUC_bycrcl_bayes_load.png", width = 8, height = 6)
  p3_sagov
  ggsave("output/AUCtarget_sagov.png", width = 8, height = 6)
  p3_auc24
  ggsave("output/AUCtarget_auc24.png", width = 8, height = 6)
  p3_bayes
  ggsave("output/AUCtarget_bayes.png", width = 8, height = 6)
  # ggsave("output/AUCtarget_bayes_addsamp.png", width = 8, height = 6)
  p3_baynl
  ggsave("output/AUCtarget_bayes_load.png", width = 8, height = 6)
  p3_all
  ggsave("output/AUCtarget_all.png", width = 8, height = 6)
  # ggsave("output/AUCtarget_all_addsamp.png", width = 8, height = 6)
  p3_all_square
  ggsave("output/AUCtarget_all_square.png", width = 8, height = 6)
  