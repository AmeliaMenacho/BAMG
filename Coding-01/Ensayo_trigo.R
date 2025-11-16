---
  title: "Heredabilidad en el fitomejoramiento"
author: "Flavio Lozano-Isla"
Modificado por: "Amelia Menacho"
format:
  html:
  toc: true
toc-location: left
number-sections: true
self-contained: true
editor_options: 
  chunk_output_type: console
execute:
  echo: false
warning: false
---

  #| label: setup
#| include: false

source("https://inkaverse.com/docs.r")
knitr::opts_chunk$set(echo = TRUE)

cran <- c("inti"
          , "metan"
          , "psych"
          , "FactoMineR"
          , "openxlsx"
          , "cowplot"
          , "agridat"
          , "dplyr"
)

suppressPackageStartupMessages({
  for (pkgs in cran) {
    if( !require(pkgs, character.only = TRUE) ) {
      install.packages(pkgs)
      library(pkgs, character.only = TRUE)
    }
  }
  remove(cran, pkgs)
})
#Base de datos
data("yan.winterwheat", package = "agridat")

set.seed(123)

fb <- yan.winterwheat %>%
  dplyr::rename(
    GEN = gen,
    ENV = env,
    KW = yield
  ) %>%
  dplyr::mutate(
    PH = KW * 20 + stats::rnorm(dplyr::n(), sd = 10),
    REP = 1L
  )

fb %>% web_table()

str(fb)

#Base de datos en excel
openxlsx::write.xlsx(x = fb, file = "met.xlsx")
getwd()

png(filename = "correlation.png", width = 25, height = 25, units = "cm", res = 300)

fb %>%
  dplyr::select(where(is.numeric)) %>%
  # eliminar columnas constantes (sd = 0)
  dplyr::select(where(~ stats::sd(.x, na.rm = TRUE) > 0)) %>%  
  pairs.panels(hist.col="red",
               pch = 21,
               method = "pearson",
               stars = TRUE,
               scale = TRUE,
               lm = TRUE
  )

graphics.off()

include_graphics("correlation.png")

mv <- fb %>%
  dplyr::group_by(GEN) %>%
  dplyr::summarise(across(where(is.numeric)
                          , ~ mean(.x, na.rm = TRUE))) %>%
  column_to_rownames("GEN") %>%
  PCA(X = ., scale.unit = TRUE, graph = FALSE)

plot1 <- plot.PCA(mv, choix = "var")
plot1

plot2 <- plot.PCA(mv,
                  choix = "ind",
                  cex = 0.8,
                  shadowtext = TRUE,
                  label = "ind",
                  autoLab = "yes",
                  graph.type = "ggplot")
plot2

list(plot1, plot2) %>%
  plot_grid(
    plotlist = .,
    nrow = 1,
    labels = "AUTO",
    rel_widths = c(1, 1.5)
  )

n_env <- length(unique(fb$ENV))

hr <- H2cal(data = fb,
            trait = "KW",
            gen.name = "GEN",
            rep.n = 1,
            env.name = "ENV",
            env.n = n_env,
            fixed.model = "0 + (1|ENV) + GEN",
            random.model = "1 + (1|ENV) + (1|GEN)",
            emmeans = TRUE,
            plot_diag = TRUE,
            outliers.rm = FALSE)

hr$model %>% summary()
hr$tabsmr %>% web_table()


hr$blups %>% arrange(desc(KW)) %>% kable()

hr <- H2cal(data = fb,
            trait = "PH",
            gen.name = "GEN",
            rep.n = 1,
            env.name = "ENV",
            env.n = n_env,
            fixed.model = "0 + (1|ENV) + (1|ENV) + GEN",
            random.model = "1 + (1|ENV) + (1|ENV) + (1|GEN)",
            emmeans = TRUE,
            plot_diag = TRUE,
            outliers.rm = TRUE)

hr$model %>% summary()
hr$tabsmr %>% web_table()

hr$blups %>% arrange(desc(PH)) %>% kable()


list(plot1, plot2) %>%
  plot_grid(
    plotlist = .,
    nrow = 1,
    labels = "AUTO",
    rel_widths = c(1, 1.6)
  ) %>%
  ggsave2("Coding-01/quantitative-genetics.png",
          plot = .,
          width = 30, height = 12, units = "cm")


