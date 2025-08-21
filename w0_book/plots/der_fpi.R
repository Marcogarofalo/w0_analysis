library(Rose)
library(ggplot2)
library(plotly)
library(knitr)

masses <- c(
  1,
  27.30,
  11.70 * 27.30
)


gg <- NULL
basenames <- c(
  "/home/garofalo/analysis/flow/data/fit_all_beta/der_fpi_diffplat_full_mu_mu2_mu3"
)
labels <- c(
  "P[0]/amu + P[1]/amu2 + P[3]/amu3 "
)
gg <- plot_fit(basenames[1], "amu",
  # noline = TRUE,
  data_type = c("B64"),
  gg = gg,
  id_x = 1,
  # single_name_for_fit = labels[j],
  noribbon = TRUE,
  noline = TRUE,
  stroke = 2,
  width = 2,
  filter_data = c(0),
)
colorlist <- c(
  # "#404040",
  "#4863A0",
  "#C04000",
  "#228B22",
  "#8B008B",
  "#00CCFF",
  "#996600", "#999999", "#FFCC33",
  "#FF6600", "#6633FF", "#9966FF",
  "#006666", "#FFCCFF", "#fc0303",
  "#03fc07", "#0335fc", "#fc03e3",
  "#d7fc03"
)
gg <- gg + scale_color_manual(values = colorlist)
gg <- gg + geom_vline(xintercept = masses[2] / masses[1], linetype = "dashed")
gg <- gg + geom_vline(xintercept = masses[3] / masses[1], linetype = "dashed")
gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$(\\mu/\\mu_\\ell)$", "$(\\mu_\\ell df_\\pi)/d\\mu$",
  to_print = FALSE,
  save_pdf = "der_fpi_B64",
  yrange = c(-0.5, 2.5),
  legend_position = c(0.7, 0.98)
)
############################################ à
#
#####################
gg <- NULL
gg <- plot_fit(basenames[1], "amu",
  # noline = TRUE,
  data_type = c("B64", "C80", "D96"),
  gg = gg,
  id_x = 1,
  # single_name_for_fit = labels[j],
  noribbon = TRUE,
  noline = TRUE,
  stroke = 2,
  width = 2,
  filter_data = c(0, 1, 2),
)
gg <- gg + scale_color_manual(values = colorlist)
gg <- gg + geom_vline(xintercept = masses[2] / masses[1], linetype = "dashed")
gg <- gg + geom_vline(xintercept = masses[3] / masses[1], linetype = "dashed")
gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$(\\mu/\\mu_\\ell)$", "$(\\mu_\\ell df_\\pi)/d\\mu$",
  to_print = FALSE,
  save_pdf = "der_fpi_all",
  yrange = c(-0.5, 2.5),
  legend_position = c(0.7, 0.98)
)
######################################
gg <- NULL
gg <- plot_fit(basenames[1], "amu",
               # noline = TRUE,
               data_type = c("B64", "C80", "D96"),
               gg = gg,
               id_x = 1,
               single_name_for_fit = "fit",
               noribbon = FALSE,
               noline = TRUE,
               stroke = 2,
               width = 2,
               filter_data = c(0, 1, 2),
)
gg <- gg + geom_vline(xintercept = masses[2] / masses[1], linetype = "dashed")
gg <- gg + geom_vline(xintercept = masses[3] / masses[1], linetype = "dashed")
gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$(\\mu/\\mu_\\ell)$", "$(\\mu_\\ell df_\\pi)/d\\mu$",
                to_print = FALSE,
                save_pdf = "der_fpi_fit",
                yrange = c(-0.5, 2.5),
                legend_position = c(0.7, 0.98)
)
##############################################################
##############################################################
##############################################################
#  chi2
##############################################################
##############################################################
basenames <- c(
  "/home/garofalo/analysis/flow/data/fit_all_beta/der_fpi_diffplat_full_mu_mu2_mu3_chi2d1"
)
######################################
gg <- NULL
gg <- plot_fit(basenames[1], "amu",
               # noline = TRUE,
               data_type = c("B64", "C80", "D96"),
               gg = gg,
               id_x = 1,
               single_name_for_fit = "fit",
               noribbon = FALSE,
               noline = TRUE,
               stroke = 2,
               width = 2,
               filter_data = c(0, 1, 2),
)
gg <- gg + geom_vline(xintercept = masses[2] / masses[1], linetype = "dashed")
gg <- gg + geom_vline(xintercept = masses[3] / masses[1], linetype = "dashed")
gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$(\\mu/\\mu_\\ell)$", "$(\\mu_\\ell df_\\pi)/d\\mu$",
                to_print = FALSE,
                save_pdf = "der_fpi_fit_chi2d1",
                yrange = c(-0.5, 2.5),
                legend_position = c(0.7, 0.98)
)

##################################
# Mpi
################################
basenames <- c(
  "/home/garofalo/analysis/flow/data/fit_all_beta/der_Mpi_diffplat_full_mu_mu2_mu3_chi2d1"
)

gg <- NULL
gg <- plot_fit(basenames[1], "amu",
               # noline = TRUE,
               data_type = c("B64", "C80", "D96"),
               gg = gg,
               id_x = 1,
               single_name_for_fit = "fit",
               noribbon = FALSE,
               noline = TRUE,
               stroke = 2,
               width = 2,
               filter_data = c(0, 1, 2),
)
gg <- gg + geom_vline(xintercept = masses[2] / masses[1], linetype = "dashed")
gg <- gg + geom_vline(xintercept = masses[3] / masses[1], linetype = "dashed")
gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$(\\mu/\\mu_\\ell)$", "$(\\mu_\\ell dM_\\pi)/d\\mu$",
                to_print = FALSE,
                save_pdf = "der_Mpi_fit_chi2d1",
                yrange = c(-0.8, 1),
                legend_position = c(0.7, 0.98)
)
