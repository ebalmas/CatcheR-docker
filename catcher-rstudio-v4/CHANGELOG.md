# Docker catcheR V4 — 21 ottobre 2026

- Ambiente RStudio Server con R 4.5.2, Bioconductor 3.22, Python 3 e Bash.
- Inclusi tutti i 100 pacchetti R richiesti dalle note del vecchio Docker,
  con le dipendenze obbligatorie e i pacchetti necessari a catcheR.
- catcheR V4 fissato al commit 724ef8900a6edb31136a7b0ccf4d3d5e3ba18dbe;
  nessuna modifica ai suoi sorgenti.
- Copia completa del repository in catcheR_V4/, incluse docs/ e vignettes/;
  nel container è consultabile da /home/rstudio/catcheR_V4.
- Unica cartella condivisa EB_analysis, accessibile come /home/rstudio/project,
  con scripts/, scratch/ e Output/.
- Predisposti i percorsi richiesti dagli script legacy tramite collegamenti.
- Rimossi i requisiti di bcl2fastq, Cell Ranger e micromamba.
- Aggiunta la libreria di sistema ImageMagick necessaria a SpatialExperiment/Banksy.
- Script per costruzione, avvio, arresto e verifiche; versioni e risultati delle
  prove raccolti in validation/.
