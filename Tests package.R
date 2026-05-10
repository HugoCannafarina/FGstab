# Charger le package en dev (sans l'installer, rechargé à chaque modif)
devtools::load_all()

# Regénérer NAMESPACE et le dossier man/ depuis tes roxygen2
devtools::document()

# Lancer les tests
devtools::test()

# Vérification complète style CRAN
devtools::check()


devtools::build_vignettes()

install.packages(c("knitr", "rmarkdown"))

# Construire et prévisualiser la vignette
devtools::install(build_vignettes = FALSE) 
devtools::build_vignettes(install = FALSE)
browseVignettes("FGstab")

knitr::knit("vignettes/getting-started.Rmd.orig",
            output = "vignettes/getting-started.Rmd")

# Générer le site pkgdown
install.packages("pkgdown")
pkgdown::build_site()