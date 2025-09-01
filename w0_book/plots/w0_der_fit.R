library(Rose)
library(ggplot2)
library(plotly)
library(knitr)
masses <- c(
  1,
  27.30,
  11.70*27.30
)
gg <- NULL
basenames <- c(
  # "/home/garofalo/analysis/flow/data/fit_all_beta/der_w0_full_mu_mu2_mu3_chi2d1",
  "/home/garofalo/analysis/flow/data/fit_all_beta/der_w0_full_mu0_mu1_chi2d1"
)
labels <- c(
  # "P[0]/amu + P[1]/amu2 + P[3]/amu3 ",
  "P[0] + P[1]*mu "
)


for (j in seq_along(basenames)) {
  print(basenames[j])
  fit <- Rose::read_fit_P_file(paste0(basenames[j], "_fit_P.dat"))
  make_table_fit_result(fit)
  gg <- plot_fit(basenames[j], "amu",
                 # noline = TRUE,
                 data_type = c("B64", "C80","D96"),
                 single_name_for_fit = "fit",
                 noribbon = FALSE,
                 noline = TRUE,
                 gg = gg,
                 id_x = 1,
                 alpha_ribbon = 0.2,
                 stroke = 2,
                 width = 2,
                 # single_name_for_fit = labels[j]
  )

}

gg <- gg + geom_vline(xintercept = masses[2] / masses[1], linetype = "dashed")
gg <- gg + geom_vline(xintercept = masses[3] / masses[1], linetype = "dashed")
gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$(\\mu/\\mu_\\ell)$", "$(\\mu_\\ell dw_0)/d\\mu[fm]$",
                to_print = FALSE,
                save_pdf = "der_w0_fit_chi2d1",
                #yrange = c(-1, 4),
                legend_position = c(0.7, 0.98)
)

