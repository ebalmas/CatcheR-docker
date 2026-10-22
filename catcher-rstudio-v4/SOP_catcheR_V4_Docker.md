# SOP — Utilizzo del Docker catcheR V4 con RStudio

**Versione:** 1.0 · **Data:** 21 ottobre 2026  
**Immagine concordata:** `hedgelab/catcherv4-hedgelab_1.0:base`  
**Sorgenti catcheR incorporati:** commit `724ef8900a6edb31136a7b0ccf4d3d5e3ba18dbe`.

## 1. Scopo

Avviare RStudio nell'ambiente catcheR V4, mantenendo sul computer tutti gli
script, gli input e gli output dell'esperimento. Si condivide una sola cartella:
quella che contiene i propri `.Rmd`, `scratch/` e `Output/`.

Il nome della cartella è libero: `analysis_EB`, `EB_analysis` o un altro nome.
In questa SOP si usa **analysis_EB**. Dentro il container la stessa cartella
apparirà sempre come **`/home/rstudio/project`**.

La procedura riguarda l'ambiente e la gestione dei file. La scelta dei parametri
biologici e di filtering va effettuata per il singolo esperimento.

## 2. Organizzazione delle cartelle sul computer

### Cartella del Docker

Conservare i file del Docker insieme, per esempio nella cartella già preparata
`make_docker/catcher-rstudio-v4/`:

```text
catcher-rstudio-v4/
  Dockerfile
  build-image.sh
  start-catcher.sh
  stop-catcher.sh
  check-image.sh
  README.md
  SOP_catcheR_V4_Docker.md
  packages-*.txt
  sources.json
  scripts/                      # strumenti per costruire/verificare il Docker
  vendor/                       # archivi sorgente usati dalla build
  catcheR_V4/                    # copia dei sorgenti GitHub
    docs/
    vignettes/
    inst/rmd/
  validation/                   # versioni e risultati delle verifiche
```

**`start-catcher.sh` resta qui.** Non occorre copiarlo dentro ogni esperimento.
La cartella `scripts/` del Docker contiene strumenti dell'ambiente: è distinta
dall'eventuale cartella `scripts/` dell'esperimento.

### Cartella dell'esperimento

Mantenere la propria organizzazione abituale:

```text
AB04_NANOG/
  AB04/
  AB04_home/
  analysis_EB/                  # condividere questa cartella intera
    n1_AB04_QC.Rmd
    n2_AB04_DEG.Rmd
    n3_AB04_MiloR.Rmd
    altri_script.R
    old_codes/
    scratch/
    Output/
    summary.pptx
```

Gli script possono stare direttamente in `analysis_EB/`. Non è obbligatorio
spostarli in una sottocartella. Si possono conservare anche documenti, vecchi
script e altre sottocartelle: saranno visibili attraverso la stessa condivisione.

| Elemento | Funzione | Persistenza |
|---|---|---|
| `.R`, `.Rmd`, `old_codes/` | I propri script e report | File sul computer |
| `scratch/` | Input dell'analisi e, se previsti, file di passaggio catcheR | File sul computer |
| `Output/` | Risultati salvati dall'analisi | File sul computer |
| `scripts/` nell'esperimento | Eventuale raccolta personale di script | Facoltativa |

Lo script di avvio crea `scratch/`, `Output/` e `scripts/` se mancano. Non sposta
gli script già presenti. `scripts/` può essere lasciata vuota.

**La condivisione è in lettura e scrittura:** modificare o cancellare un file in
`/home/rstudio/project` modifica o cancella lo stesso file sul computer.

## 3. Preparazione iniziale dell'immagine — una volta

Avviare Docker Desktop e attendere che sia pronto. Aprire il Terminale nella
cartella `catcher-rstudio-v4`.

I comandi seguenti vanno eseguiti nel **Terminale del computer**. Usare `bash`
per gli script forniti.

### A. Costruire con il nome concordato

Se la build non è ancora stata eseguita:

```bash
export CATCHER_IMAGE="hedgelab/catcherv4-hedgelab_1.0:base"
bash build-image.sh
```

Attendere il completamento senza errori. La compilazione può richiedere tempo.
Non serve ricostruire l'immagine per ogni esperimento o per modificare uno script.

### B. Se è già stata costruita con il nome predefinito

Se era stato eseguito `bash build-image.sh` senza impostare `CATCHER_IMAGE`,
il nome predefinito è `hedgelab/catcher-rstudio:v4-724ef89`. Per assegnare il
nome concordato, senza ricompilare:

```bash
docker tag hedgelab/catcher-rstudio:v4-724ef89 \
  hedgelab/catcherv4-hedgelab_1.0:base
```

Se si è già costruita l'immagine con il nome concordato, questo passaggio non serve.

### C. Pubblicare su Docker Hub — quando desiderato

Dopo una build completata:

```bash
docker login --username hedgelab
docker push hedgelab/catcherv4-hedgelab_1.0:base
```

Il push pubblica l'immagine; non include i dati personali della cartella
`analysis_EB`. Il Dockerfile incorpora i componenti esplicitamente previsti
nella cartella di build. Non aggiungere dati dell'esperimento a tali componenti.

### D. Su un altro computer, oppure dopo aver rimosso l'immagine locale

Quando l'immagine è stata pubblicata su Docker Hub:

```bash
docker pull --platform linux/amd64 hedgelab/catcherv4-hedgelab_1.0:base
```

Serve anche una copia della cartella con gli script di avvio. In alternativa,
se l'immagine non è stata pubblicata, ricostruirla dai file della cartella Docker.

La piattaforma prevista è `linux/amd64`, come per l'ambiente CNV. Su Mac Apple
Silicon viene eseguita attraverso Docker Desktop.

## 4. Avvio quotidiano di un esperimento

1. Avviare Docker Desktop.
2. Individuare la cartella che contiene davvero gli `.Rmd`, `scratch/` e `Output/`.
3. Aprire il Terminale nella cartella del Docker.
4. Sostituire il percorso di esempio con il percorso completo dell'esperimento:

```bash
export CATCHER_IMAGE="hedgelab/catcherv4-hedgelab_1.0:base"
bash start-catcher.sh "/percorso/completo/AB04_NANOG/analysis_EB"
```

Le virgolette vanno mantenute, soprattutto se il percorso contiene spazi.
Controllare il percorso prima di avviare: se si indica una cartella inesistente,
lo script può crearla invece di aprire l'esperimento desiderato.

5. Scegliere la password quando richiesta. Durante la digitazione non compaiono
   caratteri: è normale.
6. Aprire nel browser **http://localhost:8081**.
7. Accedere con utente **rstudio** e la password appena scelta.

`CATCHER_IMAGE` va impostata nuovamente quando si apre un nuovo Terminale,
perché lo script mantiene il nome predefinito originale.

### Avvio da un'altra cartella del Terminale

È possibile richiamare lo script usando il suo percorso completo, senza copiarlo
nell'esperimento. Sostituire entrambi i percorsi:

```bash
export CATCHER_IMAGE="hedgelab/catcherv4-hedgelab_1.0:base"
bash "/percorso/catcher-rstudio-v4/start-catcher.sh" \
  "/percorso/completo/AB04_NANOG/analysis_EB"
```

## 5. Controllo iniziale dentro RStudio

Nella **console R di RStudio**:

```r
setwd("/home/rstudio/project")
getwd()
list.files()

library(catcheR)
library(Seurat)
packageVersion("catcheR")
packageVersion("Seurat")
```

Verificare che `list.files()` mostri i propri script, `scratch/` e `Output/`.
Se compaiono altri file o una cartella quasi vuota, controllare quale cartella
è stata condivisa prima di procedere con l'analisi.

| Sul computer | Dentro RStudio/container |
|---|---|
| `…/analysis_EB/` | `/home/rstudio/project/` |
| `…/analysis_EB/scratch/` | `/home/rstudio/project/scratch/` |
| `…/analysis_EB/Output/` | `/home/rstudio/project/Output/` |
| `…/analysis_EB/n1_AB04_QC.Rmd` | `/home/rstudio/project/n1_AB04_QC.Rmd` |

I percorsi macOS che iniziano con `/Users/…` non sono i percorsi da usare dentro
questo container. Per i file condivisi usare `/home/rstudio/project/…`, oppure
percorsi relativi a questa cartella dopo aver impostato la directory di lavoro.

## 6. Utilizzare i propri script o copiare un template

Aprire i propri `.Rmd` direttamente da `/home/rstudio/project` tramite il pannello
Files. Non devono essere nella sottocartella `scripts/`.

Per iniziare da un template catcheR originale:

| Chimica | Template |
|---|---|
| 10X, flusso nativo R | `catcheR_tenxcatch_analysis_EB.Rmd` |
| Parse | `catcheR_parsecatch_analysis_EB.Rmd` |
| PIPseq | `catcheR_pipcatch_analysis_EB.Rmd` |

Esempio nella console R, scegliendo il template 10X e copiandolo nella radice
dell'esperimento:

```r
template <- "catcheR_tenxcatch_analysis_EB.Rmd"
file.copy(
  system.file("rmd", template, package = "catcheR"),
  file.path("/home/rstudio/project", template),
  overwrite = FALSE
)
```

Se il file esiste già, non viene sovrascritto. Usare un altro nome per creare una
nuova copia, oppure aprire la propria copia esistente.

I sorgenti di riferimento e la documentazione sono consultabili in:

```text
/home/rstudio/catcheR_V4/docs/
/home/rstudio/catcheR_V4/vignettes/
/home/rstudio/catcheR_V4/inst/rmd/
```

Questa copia risiede nell'immagine. Per conservare le proprie modifiche,
lavorare sulle copie salvate in `/home/rstudio/project`.

## 7. Impostare cartelle e input nei template catcheR

Nel blocco YAML già presente all'inizio del proprio `.Rmd`, modificare i
parametri corrispondenti, mantenendo gli altri:

```yaml
params:
  analysis_dir: "/home/rstudio/project"
  scratch: "scratch"
  output: "Output"
```

Non aggiungere un secondo blocco `params:`. La voce `output` dentro `params`
indica la cartella dei risultati catcheR; è distinta dalla voce YAML principale
`output: html_document`, che sceglie il formato del report.

Impostare anche `sample_name`, i nomi dei FASTQ, i riferimenti barcode-gene e
gli altri parametri della propria chimica. I nomi inclusi nei template sono
esempi e devono essere adattati ai propri file.

Salvare gli input in `scratch/` e indicare nei parametri i nomi o i percorsi
relativi corretti. Se si conservano sottocartelle, per esempio `scratch/10X/`,
aggiornare anche `scratch: "scratch/10X"`, oppure i percorsi dei file secondo
il template; non basta spostare i file in Finder.

### Distinguere i due tipi di input

| Chimica | Identificazione dei barcode cellulari nell'arricchimento | Espressione per l'ingresso in Seurat |
|---|---|---|
| 10X | FASTQ, riferimento barcode-gene; lista di cellule/whitelist configurata nel template, anche da `.h5` | `.h5` Cell Ranger, oppure terzetto MTX/features/barcodes |
| Parse | FASTQ, riferimento barcode-gene e `barcode_data.csv` | `count_matrix.mtx.gz`, `all_genes.csv.gz`, `cell_metadata.csv.gz` |
| PIPseq | FASTQ e riferimento barcode-gene; riconoscimento dei segmenti cellulari in R1 | Terzetto DRAGEN/PIPseeker MTX/features/barcodes, anche con prefisso del campione |

Il file di espressione non sostituisce i FASTQ dell'arricchimento CATCHER.
Usare gli export del campione e della chimica corrispondenti. Per i dettagli
consultare `catcheR_V4/docs/data-flow.md` e la guida della propria chimica.

### Caso 10X: si dispone di un solo file `.h5`

Il lettore nativo accetta un `.h5` Cell Ranger come lista per i barcode e come
matrice per Seurat. Nel parametro `cell_whitelist` indicare il nome del proprio
`.h5` presente in `scratch/`.

**Il template 10X di questo snapshot cerca automaticamente il terzetto MTX
nel blocco Seurat.** Mettere soltanto un `.h5` in scratch non attiva quel blocco.
Nella propria copia del `.Rmd`, sostituire il contenuto del chunk `seurat` con
una chiamata esplicita, adattando il nome del file:

```r
gex_h5 <- file.path(SCRATCH, "filtered_feature_bc_matrix.h5")
stopifnot(file.exists(gex_h5))

obj <- catcheR_tenx_seurat(
  res, filtered = flt, empty = noc,
  mtx = gex_h5,
  sample_id = params$sample_name,
  project = params$sample_name,
  rename = TRUE,
  subset.valid = TRUE,
  output_dir = OUT,
  verbose = TRUE
)
```

Questo esempio usa gli oggetti `res`, `flt`, `noc`, `SCRATCH` e `OUT` creati
nei blocchi precedenti del template. Non eseguirlo prima di quei blocchi.
La modifica riguarda la copia personale dello script, non i sorgenti del pacchetto.

## 8. Eseguire, salvare e controllare i risultati

1. Salvare il proprio `.Rmd` dentro `analysis_EB`.
2. Verificare i percorsi e i nomi del campione.
3. Eseguire i blocchi in ordine, oppure usare Knit quando il documento è pronto.
4. Controllare gli output e il manifesto degli input scritto da catcheR.
5. Salvare gli oggetti e i file necessari al downstream dentro la cartella condivisa.

Nei flussi nativi catcheR l'ordine è estrazione/QC → filtering → nocatch → Seurat.
Il passaggio fra questi step avviene tramite gli oggetti R. Alcuni file di
passaggio possono essere copiati anche in scratch quando `to_scratch: true`.

Gli output catcheR seguono la configurazione del template. Il report HTML
prodotto da Knit può essere scritto accanto al `.Rmd`: non tutti i file generati
vanno automaticamente in `Output/`.

Per salvare esplicitamente un oggetto personale:

```r
saveRDS(obj, "/home/rstudio/project/Output/mio_oggetto_seurat.rds")
```

Per più campioni della stessa chimica usare cartelle di output separate, per
esempio `output: "Output/AB019a"` e `output: "Output/AB020a"` nelle rispettive
copie del template. I nomi di alcune figure sono condivisi: cambiare solo il
nome del `.Rmd` non separa i risultati. Per conservare prove diverse dello stesso
campione, usare cartelle di output distinte anche tra quelle prove.

Al termine salvare il codice e, per tracciare l'ambiente:

```r
writeLines(capture.output(sessionInfo()),
           "/home/rstudio/project/Output/sessionInfo.txt")
```

## 9. Fermare e riaprire l'ambiente

Prima di fermarlo, salvare script e oggetti necessari. Il contenuto della memoria
R non viene conservato automaticamente dal container.

Dal Terminale del computer, nella cartella Docker:

```bash
bash stop-catcher.sh
```

In alternativa, da qualsiasi cartella:

```bash
docker stop catcher-rstudio-v4
```

Il container temporaneo viene rimosso. I file in `analysis_EB` e l'immagine
Docker restano disponibili. Per riprendere, ripetere l'avvio del punto 4.

I file salvati esclusivamente fuori da `/home/rstudio/project` non sono coperti
da questa condivisione. Anche eventuali pacchetti aggiunti soltanto nel container
temporaneo vanno predisposti nell'immagine per essere disponibili ai prossimi avvii.

### Passare a un altro esperimento

Fermare il container, quindi rilanciare `start-catcher.sh` indicando il percorso
della nuova cartella. Il nome del container predefinito è unico, quindi lo stesso
avvio non può essere duplicato mentre è già in esecuzione.

### Scegliere un'altra porta

Per usare, per esempio, la porta 8082:

```bash
export CATCHER_IMAGE="hedgelab/catcherv4-hedgelab_1.0:base"
bash start-catcher.sh "/percorso/completo/analysis_EB" 8082
```

Aprire quindi **http://localhost:8082**.

## 10. Problemi comuni

| Situazione | Controllo o azione |
|---|---|
| Docker non risponde | Avviare Docker Desktop e attendere che sia pronto. |
| Immagine non trovata | Controllare `CATCHER_IMAGE` e i nomi con `docker image ls`; costruire, assegnare il tag o scaricare l'immagine. |
| Il container esiste già | Riaprire RStudio se è quello desiderato; altrimenti fermarlo prima di avviare un altro esperimento. |
| La porta 8081 è occupata | Scegliere una porta diversa con il secondo argomento dello script. |
| In RStudio non compaiono i propri script | Controllare `/home/rstudio/project` e la cartella selezionata nell'avvio. |
| R segnala che `/Users/…` non esiste | Usare il percorso interno `/home/rstudio/project/…` per i file condivisi. |
| Si usa soltanto un `.h5` e Seurat viene saltato | Adattare il blocco Seurat 10X come indicato al punto 7. |
| I report di campioni diversi si sostituiscono | Impostare cartelle di output distinte nei rispettivi script. |

## 11. Verifiche disponibili e limite noto

Dalla cartella Docker:

```bash
export CATCHER_IMAGE="hedgelab/catcherv4-hedgelab_1.0:base"
bash check-image.sh
```

Sul Docker preparato sono passati il caricamento delle librerie richieste,
i flussi sintetici 10X/Parse/PIPseq verso Seurat, la lettura HDF5, Python/Bash,
l'avvio della pagina RStudio e la persistenza dei file con un unico mount.
La suite `testthat` inclusa è passata.

Lo script originale `catcheR_V4/tests/test_catcheR.R` contiene un errore di
sintassi alla riga 217. Il controllo completo lo esegue per ultimo e termina
con codice 1 per questo problema. Il file è rimasto identico al repository.
Il dettaglio è in [validation/VALIDATION.md](validation/VALIDATION.md).
Queste prove non costituiscono una validazione biologica sui propri campioni.

## 12. Promemoria operativo

**Una volta:** costruire o scaricare l'immagine e conservare la cartella Docker.

**Ogni sessione:** avviare Docker Desktop → indicare l'immagine concordata →
avviare lo script con la propria cartella di analisi → aprire RStudio → lavorare
in `/home/rstudio/project` → salvare → fermare il container.

**Le cartelle dell'esperimento restano proprie. Gli script possono stare nella
radice. La condivisione comprende tutta analysis_EB.**
