#!/usr/bin/env bash
set -euo pipefail
docker stop "${CATCHER_CONTAINER:-catcher-rstudio-v4}"
