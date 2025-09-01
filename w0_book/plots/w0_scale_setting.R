library(Rose)
library(ggplot2)
library(plotly)
library(knitr)

masses <- c(
  1,
  27.30,
  11.70 * 27.30
)
dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all_rew/"
namefit <- paste0(dir, "w0_a_A12_noC20_rew")
file <- paste0(namefit, "_fit_P.dat")
fit <- read_fit_P_file(file)
cat("\n\n")
a_w0 <- fit$P[c(1:5),2]
da_w0 <- fit$P[c(1:5),3]

# dt <- make_table_fit_result(fit)

# print(fit$P)
dt <- data.frame(
  "par" =c("a(A)","a(B)","a(C)","a(D)","a(E)","P[5]"), "a_from_w0" = mapply(mean_print, fit$P[, 2], fit$P[, 3])
  # ,"percent" = fit$P[, 3] / fit$P[, 2]
)
a_tau <- c(9.025081e-02, 7.947758e-02, 6.818945e-02, 5.685040e-02, 4.891721e-02)
da_tau <- c(7.890561e-04, 1.081341e-04, 1.429821e-04, 9.017466e-05, 1.065579e-04)
dt$a_from_fpi <- c(
  mapply(mean_print, a_tau, da_tau),
  rep(0, length(dt[, 1]) - length(a_tau))
)
cat("$\\chi^2/dof=$ ", fit$chi2, "\n\n")
kable(dt)
gg <- myggplot(repeat_color = 1)
gg <- plot_fit(
  basename = namefit,
  var = "xi",
  id_x = 5,
  data_type = c(
    "A", "B", "C", "D", "E"
    # ,"B", "C", "D"
  ),
  width = 1e-4,
  gg = gg,
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE
)
# filed <- paste0(namefit, "_fit_data.txt")
# df <- read.table(filed, header = FALSE, fill = TRUE)
# idy <- ncol(df) - 2
# kable(df[, c(1,2,3,4,idy,idy+1)],col.names =c("mu","aMpi","afpi","L","afpi-inf","err") )
gg <- gg + geom_vline(xintercept = 0.00678723)
fig <- myplotly(gg, "", "$\\xi$", "$w_0/a$", to_print = F,
                save_pdf = "w0_scale_setting", legend_position = c(0.7, 0.8))
####################################

####################################################################################
# rew analysis fpi
######################################################################

dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all_rew/"
namefit <- paste0(dir, "afpi_max_twist_A12_noC20_rew")
file <- paste0(namefit, "_fit_P.dat")
fit <- read_fit_P_file(file)
cat("\n\n")
dt <- make_table_fit_result(fit)
a_fpi <- fit$P[c(1:5),2]
da_fpi <- fit$P[c(1:5),3]
dt<- data.frame("par"=fit$P[,1], "value"= fit$P[,2],
                "error"=fit$P[,3])

gg <- myggplot()

gg <- gg + geom_pointrange(aes(
  x = a_fpi^2, y = a_fpi-a_w0,
  ymin = a_fpi - da_fpi-a_w0,
  ymax = a_fpi + da_fpi-a_w0,
  color = "w0",
  shape = "w0",
  fill = "w0",
  linewidth = "w0",
  size = "w0"
))

gg<- gg +  scale_size_manual(values = c(0.5, 0.5, 0.5))   
gg<- gg +guides(size = "none") 
gg<- gg +  scale_linewidth_manual(values = c(1.5, 1.5, 1.5))   
gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$a^2(f_\\pi)$", "$a(w_0)-a(f_\\pi)$", to_print = F, 
                save_pdf = "scaling_a_w0_fpi", xrange = c(0, 0.0085))