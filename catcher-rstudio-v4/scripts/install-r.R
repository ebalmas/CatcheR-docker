options(timeout=1800, Ncpus=as.integer(Sys.getenv("NCPUS","4")))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
repo <- Sys.getenv("CRAN_REPO")
options(repos=c(CRAN=repo))
read_packages <- function(file) scan(file,what=character(),quiet=TRUE,comment.char="#")
stage <- commandArgs(trailingOnly=TRUE)[1]
cran <- read_packages("/opt/catcher-build/packages-cran.txt")
bioc <- read_packages("/opt/catcher-build/packages-bioconductor.txt")
# Install CRAN first from the recorded Linux binary snapshot.
if(stage=="cran") {
missing <- setdiff(cran,rownames(installed.packages()))
if(length(missing)) install.packages(missing,dependencies=NA)
missing <- setdiff(cran,rownames(installed.packages()))
if(length(missing)) stop("Missing CRAN packages: ",paste(missing,collapse=", "))
}
if(stage=="bioc") {
BiocManager::install(version=Sys.getenv("BIOCONDUCTOR_VERSION"),ask=FALSE,update=FALSE)
options(repos=BiocManager::repositories(version=Sys.getenv("BIOCONDUCTOR_VERSION")))
missing <- setdiff(bioc,rownames(installed.packages()))
if(length(missing)) BiocManager::install(missing,ask=FALSE,update=FALSE,
                                       dependencies=NA,Ncpus=getOption("Ncpus"))
missing <- setdiff(bioc,rownames(installed.packages()))
if(length(missing)) stop("Missing Bioconductor packages: ",paste(missing,collapse=", "))
}
if(stage=="github") {
# Fixed source archives: installing these must not track a moving GitHub branch.
for(pkg in c("scMAGeCK","SeuratWrappers","catcheR")) {
  tarball <- paste0("/opt/catcher-vendor/",pkg,".tar.gz")
  status <- system2(file.path(R.home("bin"),"R"),
    c("CMD","INSTALL","--no-multiarch",shQuote(tarball)))
  if(status!=0L || !requireNamespace(pkg,quietly=TRUE))
    stop("Source installation failed: ",pkg)
}
}
