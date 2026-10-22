#!/usr/bin/env bash
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
image="${CATCHER_IMAGE:-hedgelab/catcher-rstudio:v4-724ef89}"
platform="${DOCKER_PLATFORM:-linux/amd64}"
builder="${DOCKER_BUILDER:-$(docker context show)}"
docker build --builder "$builder" --platform "$platform" \
  --build-arg NCPUS="${NCPUS:-4}" --progress=plain --tag "$image" "$here"
