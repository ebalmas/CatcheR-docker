suppressPackageStartupMessages(library(catcheR))
suppressPackageStartupMessages(library(Matrix))
write_v4_fastq <- function(seqs, path) {
  con <- gzfile(path, "wt"); on.exit(close(con))
  writeLines(as.vector(rbind(paste0("@read", seq_along(seqs)), seqs, "+", vapply(nchar(seqs), function(n) strrep("I",n), ""))), con)
}

run_v4_toy <- function(chem, root, plots = FALSE, called.cells = NULL) {
  scratch <- file.path(root, "scratch"); dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  rc <- catcheR:::.cat_rc
  barcode <- "GATTACAA"; ref <- "GGCGCGTTCATCTGGGGGAGCCG"; const <- "ACGCGTAGCTCGCTGATCAG"
  cassette <- paste0(rc(ref), rc(rep(c("AAAAAA", "CCCCCC"),c(8,2))), rc(barcode), const)
  cassette <- c(cassette, rep(rc("TACGCGTTCATCTGGGGGAGCCG"), 3))
  umi <- vapply(seq_along(cassette), function(i) paste0(strrep("A",8), paste(c("A","C","G","T")[(i %/% c(16,4,1)) %% 4 + 1],collapse="")), "")
  if (chem == "10X") {
    cb <- c(rep("AAAACCCCGGGGTTTT",10),rep("TTTTGGGGCCCCAAAA",3))
    cr <- paste0(cb, "A", umi)
    r1 <- cr; r2 <- cassette
  } else if (chem == "Parse") {
    b1 <- c(rep("AAAAAAAA",10),rep("TTTTTTTT",3))
    cr <- paste0(substr(umi,1,10), "GGGGGGGG", "ATGAGGGGTCAG", "CCCCCCCC", "TCCAACCACCTC", b1, "TTTTTTTTTT")
    r1 <- paste0("N", cassette); r2 <- cr
    writeLines(c("# toy whitelist", "bci,sequence,uid,well,plate,stype,wind,bc_round",
      "1,AAAAAAAA,p1,A1,1,R,1,1", "2,TTTTTTTT,p2,A2,1,R,2,1",
      "1,CCCCCCCC,p3,A1,1,L,1,2", "1,GGGGGGGG,p4,A1,1,L,1,3"), file.path(scratch,"barcode_data.csv"))
  } else {
    cr <- paste0(c(rep("AAAACCCC",10),rep("CCCCAAAA",3)), "ATG", "GGGTTT", "GAG", "ACACAC", "TCGAG", "TTTTGGGG", "ACGACG", strrep("T",12))
    r1 <- cr; r2 <- cassette
  }
  write.table(data.frame(barcode, "GENE_1"), file.path(scratch,"rc_barcodes_genes.csv"), sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)
  f1 <- paste0(chem,"_r1.fastq.gz"); f2 <- paste0(chem,"_r2.fastq.gz")
  write_v4_fastq(r1, file.path(scratch,f1)); write_v4_fastq(r2, file.path(scratch,f2))
  args <- list(folder=scratch, output_dir=file.path(root,"output",chem), subdir=FALSE, fastq.read1=f1, fastq.read2=f2, sample_name="toy", UMI.count=1, plots=plots, verbose=FALSE)
  if (!is.null(called.cells)) args$called.cells <- called.cells
  if (chem=="10X") args$cell.whitelist <- unique(cb)
  if (chem=="Parse") args$max.edit.distance <- 0
  fn <- switch(chem, `10X`=catcheR_tenxcatch, Parse=catcheR_parsecatch, PIP=catcheR_pipcatch)
  do.call(fn,args)
}




gzip_fixture <- function(path) {
  out <- paste0(path,".gz")
  input <- file(path,"rb"); data <- readBin(input,"raw",n=file.info(path)$size);close(input)
  con <- gzfile(out,"wb");writeBin(data,con);close(con)
  out
}
root <- tempfile("catcher-container-smoke-");dir.create(root)
on.exit_cleanup <- function() unlink(root,recursive=TRUE)
for(chem in c("10X","Parse","PIP")) {
  folder <- file.path(root,chem);dir.create(folder)
  x <- run_v4_toy(chem,folder)
  fn <- switch(chem,"10X"=catcheR_tenx_filtercatch,Parse=catcheR_parse_filtercatch,PIP=catcheR_pip_filtercatch)
  args <- list(x=x,UMI.count=2,plots=FALSE,verbose=FALSE)
  if(chem!="PIP") args$min.cell.UMI <- 0
  f <- do.call(fn,args)
  stopifnot(length(f$valid_cells)==1L)
  empty <- catcheR_pip_nocatch(x,filtered=f,verbose=FALSE)
  stopifnot(nrow(empty$empty_cells)==1L)
  ids <- c(unique(x$tab$cellID),unique(x$empty_reads$cellID))
  mat <- Matrix::Matrix(matrix(1:6,nrow=3),sparse=TRUE)
  mtx <- file.path(folder,"matrix.mtx");features <- file.path(folder,"features.tsv");cells <- file.path(folder,"barcodes.tsv")
  if(chem=="Parse") {
    mtx <- file.path(folder,"count_matrix.mtx");features <- file.path(folder,"all_genes.csv");cells <- file.path(folder,"cell_metadata.csv")
    Matrix::writeMM(Matrix::t(mat),mtx)
    write.csv(data.frame(gene_id=paste0("g",1:3),gene_name=paste0("G",1:3)),features,row.names=FALSE)
    write.csv(data.frame(bc_wells=ids,sample="s1",gene_count=3,tscp_count=c(6,15)),cells,row.names=FALSE)
    mtx <- gzip_fixture(mtx);features <- gzip_fixture(features);cells <- gzip_fixture(cells)
    obj <- catcheR_parse_seurat(x,filtered=f,empty=empty,mtx=mtx,genes=features,cells=cells,verbose=FALSE)
  } else {
    Matrix::writeMM(mat,mtx)
    write.table(data.frame(paste0("g",1:3),paste0("G",1:3)),features,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    writeLines(if(chem=="10X") paste0(ids,"-1") else ids,cells)
    mtx <- gzip_fixture(mtx);features <- gzip_fixture(features);cells <- gzip_fixture(cells)
    seurat_fn <- if(chem=="10X") catcheR_tenx_seurat else catcheR_pip_seurat
    obj <- seurat_fn(x,filtered=f,empty=empty,mtx=mtx,features=features,cells=cells,verbose=FALSE)
  }
  stopifnot(ncol(obj)==2L,sum(obj$ValidCell)==1L,sum(obj$is_empty)==1L)
  if(chem=="10X") {
    # Cell Ranger v3 HDF5 layout: exercise the actual .h5 reader and barcode join.
    h5path <- file.path(folder,"filtered_feature_bc_matrix.h5")
    counts <- as(mat,"dgCMatrix")
    h5 <- hdf5r::H5File$new(h5path,"w")
    g <- h5$create_group("matrix")
    g[["data"]] <- counts@x;g[["indices"]] <- counts@i;g[["indptr"]] <- counts@p
    g[["shape"]] <- dim(counts);g[["barcodes"]] <- paste0(ids,"-1")
    feat <- g$create_group("features")
    feat[["id"]] <- paste0("g",1:3);feat[["name"]] <- paste0("G",1:3)
    feat[["feature_type"]] <- rep("Gene Expression",3)
    feat[["genome"]] <- rep("GRCh38",3)
    h5$close_all()
    objh5 <- catcheR_tenx_seurat(x,filtered=f,empty=empty,mtx=h5path,verbose=FALSE)
    stopifnot(ncol(objh5)==2L,sum(objh5$ValidCell)==1L,sum(objh5$is_empty)==1L)
    message("10X HDF5 -> Seurat: OK")
  }
  message("PIPELINE OK: ",chem," FASTQ -> extraction -> filter -> nocatch -> Seurat")
}
on.exit_cleanup()
