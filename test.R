

# Preambles ---------------------------------------------------------------

# install.packages("devtools")
# install.packages("usethis")
# install.packages("rmarkdown")
# install.packages("fs")

library(devtools)
library(roxygen2)
library(knitr)
library(rmarkdown)
library(usethis)
# library(fs)


# use_git_config(
#   user.name = "lmaowisc",
#   user.email = "lmao@biostat.wisc.edu"
# )

# use_git()

# usethis:::use_devtools()


# Create, edit, and load code ---------------------------------------------
## create R file in "./R"
use_r("mi2level.R")
use_r("misc.R")
## Test your function in the new package
### devtools::load_all()
### Ctrl+Shift+L
load_all()



# Check package -----------------------------------------------------------
## 3 types of messages
## • ERRORs: Severe problems - always fix.
## • WARNINGs: Problems that you should fix, and must fix if you’re planning to
## submit to CRAN.
## • NOTEs: Mild problems or, in a few cases, just an observation.
## • When submitting to CRAN, try to eliminate all NOTEs.
check()



# Licenses ----------------------------------------------------------------

use_ccby_license()


# The DESCRIPTION file ----------------------------------------------------

# Type: Package
# Package: poset
# Title: Analysis of Partially Ordered Data
# Version: 1.0
# Author: Lu Mao
# Maintainer: Lu Mao <lmao@biostat.wisc.edu>
#   Description: Win ratio
# License: CC BY 4.0
# URL: https://sites.google.com/view/lmaowisc/
#   Depends:
#   R (>= 3.10)
# Suggests: knitr, rmarkdown
# VignetteBuilder:
#   knitr
# Config/testthat/edition: 3
# Encoding: UTF-8
# LazyData: true
# RoxygenNote: 7.3.1


# Commit changes to git ---------------------------------------------------

# Prerequisites:
# • GitHub account
# • create_github_token() - follow instructions
# • gitcreds::gitcreds_set() - paste PAT
# • git_sitrep() - verify
# use_github()
# - push content to new repository on GitHub

# usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)

edit_r_environ()
usethis::edit_r_environ("project")
###### GitHub #######
create_github_token()
gitcreds::gitcreds_set()
git_sitrep()

use_github()
####################

# Documentation -----------------------------------------------------------

# RStudio: Code > Insert Roxygen Skeleton
## Ctrl+Alt+Shift+R
# Special comments (#') above function
#   definition in R/*.R

#### ----------- an example ---
#' Multiplicative win ratio (MWR) regression analysis
#'
#' @description Fit a multiplicative win ratio (MWR) regression model to
#' partially ordered outcome against covariates
#' @return An object of class \code{MWRreg} with the following components:
#' \item{beta}{A vector of estimated regression coefficients.}
#' @seealso \code{\link{wprod}}
#' @export
#' @importFrom utils combn
#' @importFrom stats complete.cases
#' @aliases MWRreg
#' @keywords MWRreg
#' @references Mao, L. (2024). Win ratio for partially ordered data.
#' Under revision.
#' @examples
#' set.seed(12345)
#' n <- 200
#' Z <- cbind(rnorm(n),rnorm(n))
#' \dontrun{
#'   use_git()
#' }
####

### Steps ###
# Go to function definition: Ctrl+.(type function name)
# • Cursor in function definition
# • Insert roxygen skeleton (Ctrl+Alt+Shift+R)
# • Complete the roxygen fields
# • document() (Ctrl+Shift+D) - create .rd files in ./man
# • ?myfunction

document()

# document() updates the NAMESPACE file with directives from Roxygen
# comments in your R code about export() and import()


# Package-level documentation ---------------------------------------------

use_package_doc()
#> ✔ Writing 'R/mypackage-package.R'
#> • Modify ‘R/mypackage-package.R’
document()
# ?poset

check() # again



# Install package to your library -----------------------------------------

install()




use_pkgdown_github_pages()
build_readme()

# Run once to configure your package to use pkgdown
usethis::use_pkgdown()
pkgdown::build_site()

usethis::use_pkgdown_github_pages()


# Work space  ------------------------------------------------------------



?mtroc

data(uts)

?uts

### VIGNETTES #######

use_vignette("win-ratio-methodology")
use_article("Two-level-mi")
use_vignette("Two-level-mi")

use_vignette("test")
uts <- readRDS("USdat.rds")

use_data(uts, overwrite = TRUE)
