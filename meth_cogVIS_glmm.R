# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# remember current directory
dir <- getwd()

# list packages to be used
pkgs <- c("dplyr", "tidyverse", # data wrangling
          "lme4", "lmerTest", # statistical models
          "performance", # model diagnostics
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


# ---- read data ----

dz <- read.table("data/dat_z.txt", header = T) %>% # cognitive data
  slice(-which(is.na(kod))) %>%
  mutate(id = paste0(kod, wave))
dd <- read.table("data/dat_d.txt", header = T) %>% # morphometric data
  slice(-which(is.na(kod))) %>%
  mutate(id = paste0(kod, wave)) %>%
  select(-kod, -wave)
dp <- read.csv("data/dat_p.csv", sep = ",") %>% # baseline perimeter data
  mutate(
    perimeter = ifelse(
      perimeter == "not_detected" , NA , ifelse(
        perimeter == "norm" , "norm" , "pathology"
      )
    )
  )

# join the data
d <- left_join(dd, dz, by = "id") %>%
  left_join( . , dp , by = "kod") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(wave = factor(wave, labels = c("w1", "w2", "w3", "w4")),
         rnfl_nasal = (od_nasal + os_nasal) / 2,
         rnfl_temporal = (od_temporal + os_temporal) / 2) %>%
  select(-id)
str(d)


# ---- observed variables ----

DVs <- list(
  raw = list(
    visual = list(
      attention = c("tmta", "pst_d"),
      executive = c("tmtb", "pst_w", "pst_c"),
      tempo = c("gpr", "gpl")
    ),
    nvisual = list(
      attention = c("ds_f", "ds_b"),
      executive = c("cowa_skp", "sim"),
      tempo = c("ftr", "ftl"),
      memory = c("ravlt_ir", "ravlt_30")
    )
  ),
  z = list(
    visual = list(
      attention = c("tmta_z", "pstd_z"),
      executive = c("tmtb_z", "pstw_z", "pstc_z"),
      tempo = c("gpr_z", "gpl_z")
    ),
    nvisual = list(
      attention = "dsb_z",
      executive = c("skp_z", "sim_z"),
      tempo = c("ftr_z", "ftl_z"),
      memory = c("ir_z", "ravlt30_z")
    )
  )
)

logDVs <- unlist(DVs)[1:7] # list of DVs of response times to be log-transformed
for (i in logDVs) d[[i]] <- log(d[[i]]) # log response times 
for (i in unlist(DVs$raw)) d[[i]] <- as.numeric(scale(d[[i]])) # scale raw scores

# list of independent variables
IVs <- list(
  rnfl = c("os_global", "od_global", "os_nasal", "od_nasal", "os_temporal", "od_temporal"),
  vep = c("op_n1p1", "ol_n1p1", "op_p1", "ol_p1"),
  per = "perimeter"
)

for (i in IVs$rnfl) d[[i]] <- as.numeric(scale(d[[i]])) # scale RNFL

# scaling function for variables from baseline only
prescale <- function(d, var) {
  d2 <- d[ d$wave == "w1" , ]
  d2 <- data.frame(kod = d2$kod, var = as.numeric(scale(d2[[var]])))
  d <- left_join(d, d2, by = "kod")
  return(d$var)
}

for (i in IVs$vep) d[[i]] <- prescale(d, i) # scale VEPs

# prepare predictors
IVs <- unlist(IVs)
cov <- list( fixed = c("age", "edu", "gen", "wave"), rand = "kod" )
d$age <- as.numeric(scale(d$age)) # scale age
d$edu <- prescale(d, "edu") # scale education


# ---- statistical models ----

# formulas
form <- list()
for (i in names(DVs)) {
  for (j in names(DVs[[i]])) {
    for (k in names(DVs[[i]][[j]])) {
      for (l in DVs[[i]][[j]][[k]]) {
        for (m in IVs) {
          form[[i]][[j]][[k]][[l]][[m]] <- as.formula(
            paste0(
              l, " ~ ", paste(cov$fixed, collapse = " + "), " + ", m , " + (1 | kod)"
            )
          )
        }
      }
    }
  }
}

# fit
mod <- list()
for (i in names(DVs)) {
  for (j in names(DVs[[i]])) {
    for (k in names(DVs[[i]][[j]])) {
      for (l in DVs[[i]][[j]][[k]]) {
        for (m in IVs) {
          mod[[i]][[j]][[k]][[l]][[m]] <- lmer(
            formula = form[[i]][[j]][[k]][[l]][[m]], data = d,
          )
        }
      }
    }
  }
}

# residual checks
tab_res <- do.call(rbind.data.frame, strsplit(names(unlist(mod)), split = ".", fixed = T)) %>%
  rename_all(funs(make.names(c("scale", "modality", "domain", "test", "predictor"))))

for (i in names(mod)) {
  for (j in names(mod[[i]])) {
    for (k in names(mod[[i]][[j]])) {
      for (l in names(mod[[i]][[j]][[k]])) {
        for(m in names(mod[[i]][[j]][[k]][[l]])) {
          tab_res[tab_res$scale == i & tab_res$modality == j & tab_res$domain == k & tab_res$test == l & tab_res$predictor == m, 'normality (p-value)'] <- check_normality(mod[[i]][[j]][[k]][[l]][[m]]) %>% as.data.frame()
          tab_res[tab_res$scale == i & tab_res$modality == j & tab_res$domain == k & tab_res$test == l & tab_res$predictor == m, 'homogeneity (p-value)'] <- check_heteroscedasticity(mod[[i]][[j]][[k]][[l]][[m]]) %>% as.data.frame()
          tab_res[tab_res$scale == i & tab_res$modality == j & tab_res$domain == k & tab_res$test == l & tab_res$predictor == m, 'collinearity (VIF > 5)'] <- check_collinearity(mod[[i]][[j]][[k]][[l]][[m]]) %>% as.data.frame() %>% slice(which(VIF > 5)) %>% nrow()
          tab_res[tab_res$scale == i & tab_res$modality == j & tab_res$domain == k & tab_res$test == l & tab_res$predictor == m, 'outliers (#)'] <- check_outliers(mod[[i]][[j]][[k]][[l]][[m]]) %>% as.data.frame() %>% slice(which(Outlier_Cook == 1)) %>% nrow()
          tab_res[tab_res$scale == i & tab_res$modality == j & tab_res$domain == k & tab_res$test == l & tab_res$predictor == m, 'converged?'] <- check_convergence(mod[[i]][[j]][[k]][[l]][[m]]) %>% as.data.frame()
          tab_res[tab_res$scale == i & tab_res$modality == j & tab_res$domain == k & tab_res$test == l & tab_res$predictor == m, 'singular?'] <- check_singularity(mod[[i]][[j]][[k]][[l]][[m]]) %>% as.data.frame()
        }
      }
    }
  }
}

write.xlsx(tab_res, "resid/01_resid_tab.xlsx", rowNames = F)

# get parameters
par <- list()
for (i in names(mod)) {
  for (j in names(mod[[i]])) {
    for (k in names(mod[[i]][[j]])) {
      for (l in names(mod[[i]][[j]][[k]])) {
        for (m in names(mod[[i]][[j]][[k]][[l]])) {
          par[[i]][[j]][[k]][[l]][[m]] <- cbind(
            t(summary(mod[[i]][[j]][[k]][[l]][[m]])[["coefficients"]][ifelse(m == "perimeter", paste0(m, "pathology"), m), ]),
            confint(mod[[i]][[j]][[k]][[l]][[m]], parm = ifelse(m == "perimeter", paste0(m, "pathology"), m))
          )
        }
      }
    }
  }
}

tabs <- list(
  beta = list(),
  t_value = list(),
  p_value = list(),
  beta_ci = list()
)

for (i in names(par)) {
  for (j in names(par[[i]])) {
    for (k in names(par[[i]][[j]])) {
      for (l in names(par[[i]][[j]][[k]])) {
        for (m in names(par[[i]][[j]][[k]][[l]])) {
          tabs$beta[[i]][[j]][[k]][[l]][[m]] <- par[[i]][[j]][[k]][[l]][[m]][ , "Estimate"]
          tabs$t_value[[i]][[j]][[k]][[l]][[m]] <- par[[i]][[j]][[k]][[l]][[m]][ , "t value"]
          tabs$p_value[[i]][[j]][[k]][[l]][[m]] <- par[[i]][[j]][[k]][[l]][[m]][ , "Pr(>|t|)"]
          tabs$beta_ci[[i]][[j]][[k]][[l]][[m]] <- paste0(
            sprintf("%.2f", round(par[[i]][[j]][[k]][[l]][[m]][ , "Estimate"], 2)), " [",
            sprintf("%.2f", round(par[[i]][[j]][[k]][[l]][[m]][ , "2.5 %"], 2)), ", ",
            sprintf("%.2f", round(par[[i]][[j]][[k]][[l]][[m]][ , "97.5 %"], 2)), "]"
          )
        }
      }
    }
  }
}

for (i in names(tabs)) {
  for (j in names(tabs[[i]])) {
    for (k in names(tabs[[i]][[j]])) {
      for (l in names(tabs[[i]][[j]][[k]])) {
        for (m in names(tabs[[i]][[j]][[k]][[l]])) {
          tabs[[i]][[j]][[k]][[l]][[m]] <- do.call(
            cbind.data.frame, tabs[[i]][[j]][[k]][[l]][[m]]
          )
        }
        tabs[[i]][[j]][[k]][[l]] <- do.call(
          rbind.data.frame, tabs[[i]][[j]][[k]][[l]]
        ) %>% rownames_to_column(var = "test")
      }
      tabs[[i]][[j]][[k]] <- do.call(rbind.data.frame, tabs[[i]][[j]][[k]]) %>%
        rownames_to_column(var = "domain") %>%
        mutate(domain = substr(domain, 1, regexpr(".", domain, fixed = T) - 1))
    }
    tabs[[i]][[j]] <- do.call(rbind.data.frame, tabs[[i]][[j]]) %>%
      rownames_to_column(var = "modality") %>%
      mutate(modality = substr(modality, 1, regexpr(".", modality, fixed = T) - 1))
  }
}

for (i in c("raw", "z")) write.xlsx(
  list(
    beta = tabs$beta[[i]],
    beta_ci = tabs$beta_ci[[i]],
    t_value = tabs$t_value[[i]],
    p_value = tabs$p_value[[i]]
  ), file = paste0(dir, "/res/mod_tab_", i, ".xlsx"), rowNames = F
)


# ---- refit without outliers for comparison ----

out <- list()
for (i in names(mod)) {
  for (j in names(mod[[i]])) {
    for (k in names(mod[[i]][[j]])) {
      for (l in names(mod[[i]][[j]][[k]])) {
        for(m in names(mod[[i]][[j]][[k]][[l]])) {
          out[[i]][[j]][[k]][[l]][[m]] <- cbind(
            mod[[i]][[j]][[k]][[l]][[m]]@frame,
            check_outliers(mod[[i]][[j]][[k]][[l]][[m]])
          )
          flag <- unique(out[[i]][[j]][[k]][[l]][[m]][out[[i]][[j]][[k]][[l]][[m]]$Outlier_Cook == 1, ]$kod)
          if (length(flag != 0)) out[[i]][[j]][[k]][[l]][[m]] <- out[[i]][[j]][[k]][[l]][[m]] %>% slice(-which(kod %in% flag))
          else out[[i]][[j]][[k]][[l]][[m]] <- mod[[i]][[j]][[k]][[l]][[m]]@frame
          rm(flag)
        }
      }
    }
  }
}

refit <- list()
for (i in names(DVs)) {
  for (j in names(DVs[[i]])) {
    for (k in names(DVs[[i]][[j]])) {
      for (l in DVs[[i]][[j]][[k]]) {
        for (m in IVs) {
          refit[[i]][[j]][[k]][[l]][[m]] <- lmer(
            formula = form[[i]][[j]][[k]][[l]][[m]],
            data = out[[i]][[j]][[k]][[l]][[m]],
          )
        }
      }
    }
  }
}

tab_res2 <- do.call(rbind.data.frame, strsplit(names(unlist(refit)), split = ".", fixed = T)) %>%
  rename_all(funs(make.names(c("scale", "modality", "domain", "test", "predictor"))))

for (i in names(refit)) {
  for (j in names(refit[[i]])) {
    for (k in names(refit[[i]][[j]])) {
      for (l in names(refit[[i]][[j]][[k]])) {
        for(m in names(refit[[i]][[j]][[k]][[l]])) {
          tab_res2[tab_res2$scale == i & tab_res2$modality == j & tab_res2$domain == k & tab_res2$test == l & tab_res2$predictor == m, 'normality (p-value)'] <- check_normality(refit[[i]][[j]][[k]][[l]][[m]]) %>% as.data.frame()
          tab_res2[tab_res2$scale == i & tab_res2$modality == j & tab_res2$domain == k & tab_res2$test == l & tab_res2$predictor == m, 'homogeneity (p-value)'] <- check_heteroscedasticity(refit[[i]][[j]][[k]][[l]][[m]]) %>% as.data.frame()
          tab_res2[tab_res2$scale == i & tab_res2$modality == j & tab_res2$domain == k & tab_res2$test == l & tab_res2$predictor == m, 'collinearity (VIF > 5)'] <- check_collinearity(refit[[i]][[j]][[k]][[l]][[m]]) %>% as.data.frame() %>% slice(which(VIF > 5)) %>% nrow()
          tab_res2[tab_res2$scale == i & tab_res2$modality == j & tab_res2$domain == k & tab_res2$test == l & tab_res2$predictor == m, 'outliers (#)'] <- check_outliers(refit[[i]][[j]][[k]][[l]][[m]]) %>% as.data.frame() %>% slice(which(Outlier_Cook == 1)) %>% nrow()
          tab_res2[tab_res2$scale == i & tab_res2$modality == j & tab_res2$domain == k & tab_res2$test == l & tab_res2$predictor == m, 'converged?'] <- check_convergence(refit[[i]][[j]][[k]][[l]][[m]]) %>% as.data.frame()
          tab_res2[tab_res2$scale == i & tab_res2$modality == j & tab_res2$domain == k & tab_res2$test == l & tab_res2$predictor == m, 'singular?'] <- check_singularity(refit[[i]][[j]][[k]][[l]][[m]]) %>% as.data.frame()
        }
      }
    }
  }
}

write.xlsx(list(original = tab_res, refit = tab_res2),
           paste0(dir, "/resid/01_resid_tab.xlsx"), rowNames = F)

par2 <- list()

for (i in names(refit)) {
  for (j in names(refit[[i]])) {
    for (k in names(refit[[i]][[j]])) {
      for (l in names(refit[[i]][[j]][[k]])) {
        for (m in names(refit[[i]][[j]][[k]][[l]])) {
          par2[[i]][[j]][[k]][[l]][[m]] <- cbind(
            t(summary(refit[[i]][[j]][[k]][[l]][[m]])[["coefficients"]][ifelse(m == "perimeter", paste0(m, "pathology"), m), ]),
            confint(refit[[i]][[j]][[k]][[l]][[m]], parm = ifelse(m == "perimeter", paste0(m, "pathology"), m), )
          )
        }
      }
    }
  }
}

tabs2 <- list(
  beta = list(),
  t_value = list(),
  p_value = list(),
  beta_ci = list()
)

for (i in names(par2)) {
  for (j in names(par2[[i]])) {
    for (k in names(par2[[i]][[j]])) {
      for (l in names(par2[[i]][[j]][[k]])) {
        for (m in names(par2[[i]][[j]][[k]][[l]])) {
          tabs2$beta[[i]][[j]][[k]][[l]][[m]] <- par2[[i]][[j]][[k]][[l]][[m]][ , "Estimate"]
          tabs2$t_value[[i]][[j]][[k]][[l]][[m]] <- par2[[i]][[j]][[k]][[l]][[m]][ , "t value"]
          tabs2$p_value[[i]][[j]][[k]][[l]][[m]] <- par2[[i]][[j]][[k]][[l]][[m]][ , "Pr(>|t|)"]
          tabs2$beta_ci[[i]][[j]][[k]][[l]][[m]] <- paste0(
            sprintf("%.2f", round(par2[[i]][[j]][[k]][[l]][[m]][ , "Estimate"], 2)), " [",
            sprintf("%.2f", round(par2[[i]][[j]][[k]][[l]][[m]][ , "2.5 %"], 2)), ", ",
            sprintf("%.2f", round(par2[[i]][[j]][[k]][[l]][[m]][ , "97.5 %"], 2)), "]"
          )
        }
      }
    }
  }
}

for (i in names(tabs2)) {
  for (j in names(tabs2[[i]])) {
    for (k in names(tabs2[[i]][[j]])) {
      for (l in names(tabs2[[i]][[j]][[k]])) {
        for (m in names(tabs2[[i]][[j]][[k]][[l]])) {
          tabs2[[i]][[j]][[k]][[l]][[m]] <- do.call(
            cbind.data.frame, tabs2[[i]][[j]][[k]][[l]][[m]]
          )
        }
        tabs2[[i]][[j]][[k]][[l]] <- do.call(
          rbind.data.frame, tabs2[[i]][[j]][[k]][[l]]
        ) %>% rownames_to_column(var = "test")
      }
      tabs2[[i]][[j]][[k]] <- do.call(rbind.data.frame, tabs2[[i]][[j]][[k]]) %>%
        rownames_to_column(var = "domain") %>%
        mutate(domain = substr(domain, 1, regexpr(".", domain, fixed = T) - 1))
    }
    tabs2[[i]][[j]] <- do.call(rbind.data.frame, tabs2[[i]][[j]]) %>%
      rownames_to_column(var = "modality") %>%
      mutate(modality = substr(modality, 1, regexpr(".", modality, fixed = T) - 1))
  }
}

for (i in c("raw", "z")) write.xlsx(
  list(
    beta = tabs2$beta[[i]],
    beta_ci = tabs2$beta_ci[[i]],
    t_value = tabs2$t_value[[i]],
    p_value = tabs2$p_value[[i]]
  ), file = paste0(dir, "/res/refit_tab_", i, ".xlsx"), rowNames = F
)
