# Verifiche Docker catcheR V4 — 21 ottobre 2026

Immagine locale: `hedgelab/catcher-rstudio:v4-724ef89`  
ID: `sha256:77e7f3f01227d3d738ef4e1cbbc140a221a5cf268ebebe24282a36d2dde68db6`  
Piattaforma: Linux amd64; dimensione locale circa 5.12 GB.

La build è terminata con codice 0. L'ambiente e i flussi sintetici indicati
sotto hanno superato i controlli. Lo script legacy di test ha invece un errore
di sintassi nello snapshot originale: non è stata dichiarata una validazione
completa della release catcheR, né sono state modificate le sue funzioni.

## Ambiente installato

- R 4.5.2, Bioconductor 3.22, Python 3.12.3, GNU Bash 5.2.21.
- Tutti i 100 pacchetti R richiesti dalle vecchie note sono nella configurazione.
- Le 111 librerie richieste esplicitamente dal nuovo ambiente si caricano;
  il manifest contiene 472 pacchetti R complessivi, incluse le dipendenze.
- Seurat 5.4.0; catcheR 0.4.0; scMAGeCK 0.99.1;
  SeuratWrappers 0.4.0.
- catcheR è fissato al commit `724ef8900a6edb31136a7b0ccf4d3d5e3ba18dbe`.
  Verificata l'identità SHA-256 di tutti i 66 file sorgente copiati nell'immagine.
- I programmi 10X Cell Ranger e bcl2fastq sono assenti; micromamba non è installato.

## Risultati

| Controllo | Esito |
|---|---|
| Build e caricamento delle 111 librerie richieste | Passato |
| Scrittura/lettura HDF5 e creazione/salvataggio oggetto Seurat | Passato |
| 10X: FASTQ gz → estrazione → filtering → nocatch → Seurat da terzetto MTX/TSV gz | Passato |
| 10X: ingresso Seurat da file HDF5 in formato Cell Ranger v3 | Passato |
| Parse: FASTQ gz → estrazione → filtering → nocatch → Seurat da MTX/CSV gz | Passato |
| PIPseq: FASTQ gz → estrazione → filtering → nocatch → Seurat da MTX/TSV gz | Passato |
| Interfaccia Python, compilazione sintattica Python e sintassi degli script Bash | Passato |
| Collegamenti dei percorsi legacy richiesti dagli script | Passato |
| Suite originale testthat | Passata: 46 aspettative |
| Pagina di accesso RStudio | HTTP 200 |
| catcheR/Seurat accessibili all'utente rstudio; docs e vignettes consultabili | Passato |
| Unico bind mount EB_analysis → /home/rstudio/project | Passato |
| Scrittura in Output, visibilità dal computer e persistenza dopo arresto | Passato |
| Percorso condiviso contenente spazi | Passato |
| Script originale tests/test_catcheR.R | Non eseguibile: errore di sintassi alla riga 217 |

## Problema conservato nello snapshot originale

`catcheR_V4/tests/test_catcheR.R:217` contiene `else` su una riga separata
dall'istruzione `if` precedente. R interrompe il parsing con `unexpected 'else'`.
L'errore è stato osservato sia nella prova preliminare sia nell'immagine finale.
Il file rimane identico a quello pubblicato sul branch V4.

Per questo `check-image.sh` termina con codice **1** quando arriva a questo
ultimo script, dopo avere completato con successo pipeline, Python/Bash e
testthat. L'avvio di RStudio e l'uso delle funzioni catcheR non dipendono
dall'esecuzione di questo file di test. Non è stata saltata o corretta
silenziosamente la prova che fallisce.

## Ambito delle prove

I dati sono piccoli e sintetici; nessun input reale di EB_analysis è stato
letto o modificato. Per consentire la prova su poche letture, il fixture passa
esplicitamente `min.cell.UMI = 0` per 10X e Parse. I parametri predefiniti del
pacchetto restano quelli dello snapshot V4, senza modifiche.

Le prove verificano l'ambiente e i collegamenti tra gli step. Non costituiscono
una verifica biologica dei campioni, di ogni metodo delle librerie installate,
dell'intero workflow shell legacy o di ogni configurazione di filtering/QC.
Il controllo HTTP verifica la pagina di accesso; la console R è stata provata
come utente rstudio nel container, senza simulare una sessione interattiva
completa del browser. Il container temporaneo di verifica è stato arrestato;
nessuna password di prova viene consegnata o salvata nei file.

## Tracce disponibili

- `docker-build.log`: build completata e verifica delle librerie.
- `check-image.log`: prove finali e errore dello script legacy.
- `check-exit-codes.json`: codici di uscita distinti.
- `rstudio-integration.json` e `rstudio-check.log`: avvio, mount e persistenza.
- `installed-r-packages.tsv` e `R-sessionInfo.txt`: versioni realmente installate.
- `image.json`, `source-integrity.json`, `source-files-sha256.json`: identità dell'immagine e dei sorgenti.
- `old-library-coverage.json`: copertura delle librerie delle vecchie note.

L'immagine è disponibile in Docker Desktop su questo computer. Lo ZIP contiene
i file per ricostruirla e i sorgenti, non un export binario dell'immagine.
Non è stata pubblicata su Docker Hub.
