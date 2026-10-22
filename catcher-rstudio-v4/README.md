# catcheR V4 con RStudio, R, Python e Bash

Guida operativa: [SOP — utilizzo e organizzazione delle cartelle](SOP_catcheR_V4_Docker.md).

Ambiente RStudio Server impostato come il Docker CNV, con R 4.5.2.
Include tutte le librerie R richieste in modify_docker.md, SeuratWrappers e
catcheR V4 al commit 724ef8900a6edb31136a7b0ccf4d3d5e3ba18dbe.
Il codice di catcheR è quello pubblicato su GitHub, senza correzioni aggiunte.

Non include bcl2fastq, Cell Ranger o micromamba. Si parte da FASTQ ed export di
espressione già prodotti dalle rispettive piattaforme; non esegue demultiplexing
BCL o allineamento dei trascrittomi.

## 1. Costruire l'immagine

Avvia Docker Desktop. Dalla cartella di questi file:

```bash
bash build-image.sh
```

Se l'immagine è già stata costruita su questo computer, passa al punto 2.
La build completa installa molti pacchetti e può richiedere tempo. Il tag locale è
`hedgelab/catcher-rstudio:v4-724ef89`. Il programma usa linux/amd64 come
il Docker CNV; su Mac Apple Silicon Docker esegue l'emulazione.

## 2. Avviare RStudio sul proprio esperimento

Condividi una sola cartella, **EB_analysis**, che contiene tutti i file della tua analisi:

```bash
bash start-catcher.sh "/percorso/EB_analysis"
```

Lo script chiede una password, senza mostrarla. Apri:

**http://localhost:8081**

Utente: **rstudio**. Password: quella appena scelta.

Per usare un'altra porta, per esempio 8082:

```bash
bash start-catcher.sh "/percorso/EB_analysis" 8082
```

L'intera cartella EB_analysis è visibile in RStudio come
`/home/rstudio/project`. I file salvati lì restano sul computer anche quando
il container si ferma. scratch contiene i FASTQ, i riferimenti e gli export
delle chemistries; Output raccoglie gli output secondo il contratto dello
snapshot V4 installato. I file che già esistono non vengono sostituiti dallo script.

```text
EB_analysis/                 # unica cartella condivisa
  scripts/                   # i tuoi script e template Rmd
  scratch/                   # input dell'analisi
  Output/                    # risultati
```

Lo script crea queste tre sottocartelle se mancano. Puoi aggiungerne altre;
tutto ciò che scrivi in EB_analysis è visibile anche dal computer.
Non occorrono mount separati per script, input e output.

La cartella di questo Docker include `catcheR_V4/`: è la copia completa del
repository al commit indicato, con `docs/`, `vignettes/`, `inst/rmd/`, test e README.
Dentro RStudio la trovi nella home come `catcheR_V4` (percorso reale
`/opt/catcheR/source`). Rimane accessibile anche quando condividi EB_analysis;
per conservare le tue modifiche agli script, copiali in `project/scripts/`.
I dati di EB_analysis sono esclusi dal contesto di costruzione dell'immagine.

## 3. Usare catcheR

Nella console R:

```r
setwd("/home/rstudio/project")
library(catcheR)
library(Seurat)
packageVersion("catcheR")
packageVersion("Seurat")
```

I template originali sono disponibili tramite il pacchetto, e anche in
`/opt/catcher/templates`. Per copiarne uno nell'esperimento senza
sovrascrivere un file esistente:

```r
dir.create("scripts", showWarnings = FALSE)
file.copy(
  system.file("rmd", "catcheR_parsecatch_analysis_EB.Rmd", package = "catcheR"),
  "scripts/catcheR_parsecatch_analysis_EB.Rmd",
  overwrite = FALSE
)
```

Scegli il template tenx, parse o pip e modifica i parametri per i tuoi file.
Nel blocco YAML iniziale imposta la cartella di analisi con un percorso assoluto,
anche se il documento è in `scripts/`:

```yaml
params:
  analysis_dir: "/home/rstudio/project"
  scratch: "scratch"
  output: "Output"
```

Mantieni gli altri parametri del template e adattali all'esperimento.
Così scratch e Output si riferiscono sempre alla radice di EB_analysis,
anche durante il Knit da RStudio. I README e `docs/data-flow.md` della copia
V4 descrivono i file richiesti dalle diverse chemistries.
Gli export RNA di Cell Ranger sono leggibili anche senza il programma Cell
Ranger installato. I formati restano quelli richiesti dalla versione V4 scelta.

Dal menu Terminal di RStudio sono disponibili Bash, GNU sed/awk/grep, gzip e
Python 3. Il Python di estrazione plasmidica usa la libreria standard. Il venv
`/opt/catcher-python` è separato da R: non servono le vecchie correzioni PATH
di micromamba. I percorsi legacy `/data/scratch` e `/home/*` richiesti
dagli script shell sono predisposti come collegamenti, senza cambiare i sorgenti.

## 4. Fermare il container

```bash
bash stop-catcher.sh
```

La cartella dell'esperimento resta intatta. Il container temporaneo viene
rimosso; i pacchetti incorporati nell'immagine restano pronti per il prossimo
avvio. Per mantenere altri pacchetti in modo riproducibile, aggiungili alla
configurazione e ricostruisci l'immagine.

## 5. Verificare l'ambiente

```bash
bash check-image.sh
```

La verifica usa piccoli dati sintetici e le suite di test incluse nello snapshot
V4. Non usa né modifica gli input di EB_analysis. I risultati della costruzione
e delle prove eseguite su questo computer sono raccolti in `validation/`.

Lo script legacy originale `catcheR_V4/tests/test_catcheR.R` contiene un errore
di sintassi alla riga 217 (`else` su una nuova riga). È conservato identico al
repository: il controllo completo lo segnala con codice di uscita 1.
Il rapporto `validation/VALIDATION.md` distingue questo problema dello script
di test dai risultati delle verifiche dell'ambiente e dei flussi catcheR.

## Librerie e versioni

- packages-cran.txt: librerie CRAN del vecchio ambiente e dipendenze catcheR.
- packages-bioconductor.txt: librerie Bioconductor, fissate alla release 3.22
  compatibile con R 4.5; terra, sf e ggrastr sono nella lista CRAN.
- packages-github.txt: scMAGeCK, SeuratWrappers e catcheR con commit fissati.
- sources.json: origine e checksum degli archivi sorgente incorporati.
- /opt/catcher-build/installed-r-packages.tsv: tutte le versioni realmente installate.
- /opt/catcher-build/R-sessionInfo.txt: ambiente R della verifica.

scMAGeCK non è disponibile nel repository Bioconductor 3.22: si usa la sua fonte
GitHub ufficiale. SeuratWrappers include il pacchetto e le dipendenze obbligatorie;
i suoi metodi che richiedono altri software opzionali mantengono i rispettivi
requisiti. I dati di esempio e le librerie non sostituiscono la prova della tua
analisi biologica.

L'immagine non viene pubblicata automaticamente su Docker Hub. Le opzioni
CATCHER_IMAGE, CATCHER_CONTAINER, DOCKER_PLATFORM, DOCKER_BUILDER e NCPUS
permettono di personalizzare nome e build senza modificare gli script.

## Riferimenti

- [Rocker RStudio](https://rocker-project.org/images/versioned/rstudio)
- [Bioconductor: versioni R supportate](https://bioconductor.org/about/release-announcements/)
- [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
- [scMAGeCK](https://github.com/weili-lab/scMAGeCK)
- [catcheR V4](https://github.com/ebalmas/catcheR/tree/V4)
