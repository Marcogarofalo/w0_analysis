library(stringr)
library(Rose)
library(ggplot2)
library(ggthemes)


plot_fit <- function(basename, var, data_type = NULL, gg = NULL, noribbon = FALSE,
         id_x = 1,
         noline = FALSE,
         labelfit = "fit", width = 0.02, size = 1,
         id_color = NULL, id_shape = NULL,
         single_name_for_fit = NULL,
         nolabel_for_fit = FALSE,
         nudge = 0, alpha_line = 1, alpha_ribbon = 0.5,
         stroke = 1,
         filter_data = NULL) {
  filed <- paste0(basename, "_fit_data.txt")
  # Add this check:
  if (!file.exists(filed)) {
    stop(paste0("CRITICAL ERROR: The file '", filed, "' does not exist in the current directory: ", getwd()))
  }
  df <- read.table(filed, header = FALSE, fill = TRUE)
  if (!is.null(filter_data)) {
    last <- ncol(df)
    df <- df[df[, last] %in% filter_data, ]
  }
  
  if (is.null(gg)) gg <- myggplot()
  idy <- ncol(df) - 2
  
  if (is.null(id_color)) {
    color_type <- as.factor(df[, idy + 2])
  } else {
    color_type <- as.factor(df[, id_color])
  }
  
  if (is.null(id_shape)) {
    shape_type <- as.factor(df[, idy + 2])
  } else {
    shape_type <- as.factor(df[, id_shape])
  }
  
  lastr <- nrow(df)
  Nfits <- c(df[1, idy + 2]:df[lastr, idy + 2])
  if (!is.null(data_type)) {
    if (length(data_type) == 1) {
      color_type <- data_type
      shape_type <- data_type
    } else {
      N <- length(which(df[, idy + 2] == 0))
      
      # color_type <- rep(data_type, each = N)
      # shape_type <- rep(data_type, each = N)
      
      color_type <- df[, idy + 2]
      shape_type <- df[, idy + 2]
      
      count <- 1
      for (n in c(1:length(df[, idy + 2]))) {
        if (n != 1) {
          if (df[n, idy + 2] != df[n - 1, idy + 2]) {
            count <- count + 1
          }
        }
        color_type[n] <- data_type[count]
        shape_type[n] <- data_type[count]
      }
    }
  }
  
  
  datalist <- list()
  mycol <- unique(paste0(labelfit, color_type))
  if (!is.null(data_type)) {
    mycol <- Nfits
    count <- 1
    mycol[1] <- paste0(labelfit, data_type[1])
    for (n in c(1:length(df[, idy + 2]))) {
      if (n != 1) {
        if (df[n, idy + 2] != df[n - 1, idy + 2]) {
          count <- count + 1
          mycol[count] <- paste0(labelfit, data_type[count])
        }
      }
    }
  }
  if (length(mycol) != length(Nfits)) {
    mycol <- paste0(labelfit, Nfits)
  }
  
  if (!is.null(single_name_for_fit)) {
    mycol <- rep(single_name_for_fit, length(Nfits))
  }
  
  # if (nolabel_for_fit) {
  #   mycol <- unique(paste0(color_type))
  #   if (length(mycol) != length(Nfits)) {
  #     mycol <- rep(color_type, length(Nfits))
  #   }
  # }
  
  if ((!noribbon) | (!noline)) {
    for (n in Nfits) {
      file <- sprintf("%s_fit_out_n%d_%s.txt", basename, n, var)
      # browser()
      n1 <- n + 1
      if (!file.exists(file)) {
        stop(paste0(
          "CRITICAL ERROR: The file '", filed,
          "' does not exist in the current directory: ", getwd()
        ))
      }
      datalist[[n1]] <- read.table(file,
                                   header = FALSE, fill = TRUE,
                                   col.names = c(paste0("x", n), paste0("fit", n), paste0("fiterr", n))
      )
      if (!noribbon) {
        gg <- gg + geom_ribbon(
          mapping = aes_string(
            x = datalist[[n1]][, 1] + nudge,
            ymin = datalist[[n1]][, 2] - datalist[[n1]][, 3],
            ymax = datalist[[n1]][, 2] + datalist[[n1]][, 3],
            fill = as.factor(mycol[n1]),
            color = as.factor(mycol[n1]),
            shape = as.factor(mycol[n1])
          ),
          alpha = alpha_ribbon
        )
      }
      if (!noline) {
        gg <- gg + geom_line(
          mapping = aes_string(
            x = datalist[[n1]][, 1] + nudge,
            y = datalist[[n1]][, 2],
            fill = as.factor(mycol[n1]),
            color = as.factor(mycol[n1]),
            shape = as.factor(mycol[n1])
          ), alpha = alpha_line
        )
      }
    }
  }
  
  gg <- gg + geom_point(
    data = df,
    mapping = aes(
      x = df[, id_x] + nudge, y = df[, idy],
      color = color_type,
      shape = shape_type,
      fill = color_type
    ),
    size = size, stroke = stroke
  )
  
  gg <- gg + geom_errorbar(
    data = df,
    mapping = aes(
      x = df[, id_x] + nudge, y = df[, idy],
      ymin = df[, idy] - df[, idy + 1],
      ymax = df[, idy] + df[, idy + 1],
      color = color_type,
      shape = shape_type,
      fill = color_type
    ),
    width = width, size = size
  )
  
  return(gg)
}

dft <- NULL
C <- -5
gg <- NULL
path <- "/home/garofalo/analysis/flow/data/fit_all_beta/"
basenames <- c(
  paste0("fit_fpi_C", C, ".000000_N3_a2"),
  paste0("fit_fpi_C", C, ".000000_N3_a2_noBOS_noBtm_noBWTI"),
  # paste0("fit_fpi_C", C, ".000000_N3_a2_a4"),
  # paste0("fit_fpi_C", C, ".000000_N3_a2_Husung0.42"),
  # paste0("fit_fpi_C", C, ".000000_N3_a2_Husung0.21"),
  paste0("fit_fpi_C", C, ".000000_N3_a2_noCOS_noCtm_noCWTI"),
  paste0("fit_fpi_C", C, ".000000_N3_a2_a4WTI"),
  paste0("fit_fpi_C", C, ".000000_N3_a2_a4tm"),
  paste0("fit_fpi_C", C, ".000000_N3_a2_a4OS"),
  paste0("fit_fpi_C", C, ".000000_N3_a2_a4_noa4WTI"),
  paste0("fit_fpi_C", C, ".000000_N3_a2_a4_noa4tm"),
  paste0("fit_fpi_C", C, ".000000_N3_a2_a4_noa4OS"),
  paste0("fit_fpi_C", C, ".000000_N3_a2_HusungWTI"),
  paste0("fit_fpi_C", C, ".000000_N3_a2_Husungtm"),
  paste0("fit_fpi_C", C, ".000000_N3_a2_HusungOS")
)
legend_name <- gsub("fit_fpi_|\\.000000", "", basenames)
legend_name <- paste0("\\verb|",legend_name,"|") 
count <- length(basenames)
df <- data.frame(
  "fit" = rep("", count),
  "res" = rep(0, count),
  "err" = rep(0, count),
  "chi2dof" = rep(0, count),
  "dof" = rep(0, count),
  "Npar" = rep(0, count),
  "Ndat" = rep(0, count),
  "mult" = rep(0, count),
  "k" = rep(0, count)
)

# extract the C value from the string
labels <- list(c("WTI", "tm", "OS"))

for (j in seq_along(basenames)) {
  # fit <- Rose::read_fit_P_file(paste0(path, basenames[j], "_fit_P.dat"))
  # df[j, 1] <- basenames[j]
  # M <- matrix(as.numeric(as.matrix(fit$C)), nrow = nrow(fit$C))
  # e <- eigen(M)
  #
  # lmin <- min(abs(e$values))
  # lmax <- max(abs(e$values))
  # k <- lmax / lmin
  # df[j, -1] <- c(fit$P[1, 2], fit$P[1, 3], fit$chi2dof, fit$dof, fit$npar, fit$ndata, 1, k)

  gg <- plot_fit(paste0(path, basenames[j]), "a2",
    # noline = TRUE,
    data_type = labels[[1]],
    gg = gg,
    id_x = 1,
    single_name_for_fit = legend_name[j],
    # labelfit = "fit",
    width = 0.0001,
    nudge = 0,
    noline = TRUE
  )
}

fpi_FLAG <- 130.5
gg <- gg + geom_hline(yintercept = fpi_FLAG, linetype = "dashed") #+annotate("text", x = 0.0002, y = fpi_FLAG-0.2, label="FLAG")
gg <- gg + geom_point(aes(x = 0, y = fpi_FLAG, color = "FLAG", shape = "FLAG", fill = "FLAG"))

###

title <-""
xlabel <-"$a^2 [\\mbox{fm}^2]$"
ylabel <-"$f_{\\pi} [\\mbox{MeV}]$"
legend_position = c(0.7, 0.98) 

if (!title == "") gg <- gg + ggplot2::ggtitle(title)
if (!xlabel == "") gg <- gg + ggplot2::xlab(xlabel)
if (!ylabel == "") gg <- gg + ggplot2::ylab(ylabel)

if (!is.null(legend_position)) gg <- gg + theme(legend.position = legend_position + c(0.18, 0))

gg<- gg + theme_matplotlib()
ggplot2::ggsave(
  filename = paste0(fpi3reg, ".tex"),
  plot = gg,
  device = tikzDevice::tikz,
  standAlone = TRUE,
  width = width / 100,
  height = height / 100,
  units = "in"
)
