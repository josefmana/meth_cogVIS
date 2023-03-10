# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# remember current directory
dir <- getwd()

# set-up all cores form MCMC computations
options(mc.cores = parallel::detectCores())

# list packages to be used
pkgs <- c("dplyr", "tidyverse", # data wrangling
          "brms", "bayestestR", # statistical models
          "ggplot2", # plotting
          "openxlsx" # saving outputs in .xlsx
          )

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# create folders for models, figures and tables
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present
sapply( c("res","resid"), function(i) if( !dir.exists(i) ) dir.create(i) )


# ---- read data sets ----

dat <- read.table("data/mono_bayes_data.txt", header = T)

# getting to a long format
#vars_orig <- names(dat[26:ncol(dat)])
vars_orig <- names(dat[26:145])
vars_orig <- vars_orig[-c(29:32, # clth
                          37:40, # veg
                          53:56, # ds_tot
                          61:64, # ravlt_b
                          65:68  # ravlt_6
                          )]
#vars_orig <- vars_orig[-which(grepl("veg", vars_orig) | grepl("ds_tot", vars_orig))]
vars <- unique(substr(vars_orig, 1, nchar(vars_orig)-3))
desc <- names(dat[c(2:5, 146:ncol(dat))]) # demo + vep
dum_names <- split(vars_orig, ceiling(seq_along(vars_orig)/4))
dum_dat <- lapply(1:length(dum_names), function(i)
  dat[ , which(colnames(dat) %in% c(desc, dum_names[[i]]))] %>%
    pivot_longer(all_of(dum_names[[i]]),
                 names_to = "wave", names_prefix = paste0(vars[[i]], "_w"),
                 values_to = vars[[i]])
  )
names(dum_dat) <- vars
anal_dat <- dum_dat[[1]]
for (i in vars) {
  anal_dat[[i]] <- dum_dat[[i]][[i]]
}
#for (i in c(vars, "age")) {
#  anal_dat[[i]] <- scale(anal_dat[[i]])
#} # z-transforming to have sensible ropes

# z-transform before concatenate
anal_dat <- anal_dat %>% mutate(wavecon = 2 * (as.numeric(wave) - 1),
                                wave = as.ordered(wave))
for (i in c(vars, desc[-c(1, 4)])) anal_dat[[i]] <- as.numeric(scale(anal_dat[[i]]))
rm(dum_dat, dum_names, vars_orig, desc)

# ---- calculate continuous models ----

fit <- list()
system.time(
  for (i in vars[1:22]) {
    if (i == "gpl" | i == "gpr" | i == "pst_c" | i == "pst_d" | i == "tmt_a" | i == "tmt_b" |
        i == "ravlt_frec15" | i == "ravlt_rec15") {
      fam <- brmsfamily("skew_normal")
    } else {
      fam <- brmsfamily("gaussian")
    }
    for (j in vars[23:28]) {
      fit[[i]][[j]] <- brm(
        as.formula(paste0(i, " ~ age + mo(wave) + ", j, " + (1 + wave | kod)")),
        family = fam, data = anal_dat,
        inits = 0, control = list(adapt_delta = .95), seed = 87542,
        file = paste0(dir, "/fit_mono/", i, "_", j, ".rds")
      )
    }
  }
)

# ---- print parameter values ----
pars <- list()
for (i in vars[23:28]) {
  for (j in vars[1:22]) {
    pars[[i]][[j]] <- describe_posterior(fit[[j]][[i]],
                                         centrality = c("median", "mean"), dispersion = T,
                                         effects = "fixed", ci = .95, ci_method = "hdi",
                                         test = c("p_direction", "rope"),
                                         rope_ci = 1, rope_range = c(-.1, .1),
                                         )
  }
}

for (i in names(pars)) {
  pars[[i]] <- do.call(rbind.data.frame, pars[[i]]) %>% rownames_to_column(var = "test")
  pars[[i]]$test <- substr(pars[[i]]$test, 1, nchar(pars[[i]]$test)-2)
}
write.xlsx(pars, "res/mono_bayes_pars.xlsx", rowNames = F)

# ---- create table ----

tab <- list()
for (i in names(pars)) {
  dum <- pars[[i]][which(grepl(i, pars[[i]]$Parameter)), ]
  tab[[i]] <- data.frame(test = dum$test,
                         median = sprintf("%.2f", round(dum$Median, 2)),
                         `95% ci` = paste0("[", sprintf("%.2f", round(dum$CI_low, 2)), ", ",
                                           sprintf("%.2f", round(dum$CI_high, 2)), "]"),
                         pd = sprintf("%.2f", round(100*dum$pd, 2)),
                         `% in ROPE` = sprintf("%.2f", round(100*dum$ROPE_Percentage, 2))
                         )
}
write.xlsx(tab, "res/mono_bayes_tab.xlsx", rowNames = F)
