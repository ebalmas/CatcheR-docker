#!/usr/bin/env bash
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
image="${CATCHER_IMAGE:-hedgelab/catcher-rstudio:v4-724ef89}"
platform="${DOCKER_PLATFORM:-linux/amd64}"
docker run --rm --platform "$platform" --entrypoint Rscript \
  --mount "type=bind,src=$here/scripts,dst=/validation,readonly" \
  "$image" /validation/smoke-pipelines.R
docker run --rm --platform "$platform" --entrypoint python3 \
  "$image" /usr/local/lib/R/site-library/catcheR/python/plasmid_final_corrected.py --help >/dev/null
docker run --rm --platform "$platform" --entrypoint bash "$image" -c \
  'set -e; for f in /usr/local/lib/R/site-library/catcheR/scripts/*.sh; do bash -n "$f"; done; test -e /home/grep.sh; test -e /home/plasmid_inter2.R; ! command -v bcl2fastq; ! command -v cellranger; echo "Python/Bash e percorsi legacy: OK"'

docker run --rm --platform "$platform" --entrypoint Rscript "$image" -e \
  'library(catcheR); testthat::test_dir("/opt/catcheR/source/tests/testthat", reporter="summary", stop_on_failure=TRUE)'
echo "Verifiche Docker, pipeline sintetiche e testthat: OK. Ora eseguo lo script legacy originale (vedi validation/VALIDATION.md)."
docker run --rm --platform "$platform" --entrypoint Rscript "$image" \
  -e 'source("/opt/catcheR/source/tests/test_catcheR.R"); if(.n_fail > 0L) stop("Legacy suite failures: ",.n_fail)' 
