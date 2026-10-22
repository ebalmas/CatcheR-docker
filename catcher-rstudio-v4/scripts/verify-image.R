options(rgl.useNULL=TRUE)
root <- "/opt/catcher-build"
read_packages <- function(file) scan(file,what=character(),quiet=TRUE,comment.char="#")
pkgs <- unique(c(read_packages(file.path(root,"packages-cran.txt")),
                 read_packages(file.path(root,"packages-bioconductor.txt")),
                 "scMAGeCK","SeuratWrappers","catcheR"))
loaded <- vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)
if(any(!loaded)) stop("Cannot load required packages: ",paste(pkgs[!loaded],collapse=", "))
manifest <- as.data.frame(installed.packages()[,c("Package","Version","Built","LibPath")])
write.table(manifest,file.path(root,"installed-r-packages.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,"R-sessionInfo.txt"))
# Package and HDF5 smoke checks; these do not alter catcheR source or validate biology.
stopifnot(packageVersion("Seurat") >= "5.0.0")
required <- c("catcheR_tenxcatch","catcheR_parsecatch","catcheR_pipcatch",
              "catcheR_threshold_sweep","catcheR_pip_nocatch","catcheR_parse_seurat")
stopifnot(all(vapply(required,function(n)is.function(getExportedValue("catcheR",n)),logical(1))))
h5 <- tempfile(fileext=".h5")
f <- hdf5r::H5File$new(h5,"w");f[["values"]] <- 1:3;f$close_all()
f <- hdf5r::H5File$new(h5,"r");stopifnot(identical(f[["values"]][],1:3));f$close_all();unlink(h5)
m <- Matrix::Matrix(matrix(c(1,0,2,3,1,0),nrow=3),sparse=TRUE)
rownames(m) <- paste0("g",1:3);colnames(m) <- c("cell1","cell2")
obj <- Seurat::CreateSeuratObject(m)
stopifnot(ncol(obj)==2L,nrow(obj)==3L)
r <- tempfile(fileext=".rds");saveRDS(obj,r);stopifnot(ncol(readRDS(r))==2L);unlink(r)
message("IMAGE CHECK OK: ",length(pkgs)," required namespaces; catcheR, Seurat and HDF5.")
