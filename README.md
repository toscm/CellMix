# README

This is a copy of the R package [CellMix_1.6.2.tar.gz](http://web.cbio.uct.ac.za/~renaud/CRAN/src/contrib/CellMix_1.6.2.tar.gz), downloaded from <http://web.cbio.uct.ac.za/~renaud/CRAN/#>.

The only changes made to the package were:

1. Package number increased from `1.6.2` to `1.6.3`
2. Switched dependency on package `BiocInstaller` to `BiocManager`
3. Updated [DESCRIPTION](DESCRIPTION), so all dependencies can be resolved in R (>= 4.2)
4. The original README from version `1.6.2` was renamed to [README_ORIGINAL.md](README_ORIGINAL.md) and this README was added instead.

## Installation

Run the command below from and interactive R session:

```R
if (!"devtools" %in% installed.packages()[, "Package"]) {
  install.packages("devtools", repos = c(CRAN = "https://cloud.r-project.org"))
}
if (!"csSAM" %in% installed.packages()[, "Package"]) {
  devtools::install_url("https://cran.r-project.org/src/contrib/Archive/csSAM/csSAM_1.2.4.tar.gz")
}
devtools::install_github("toscm/CellMix")
```