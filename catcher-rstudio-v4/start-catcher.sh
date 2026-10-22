#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Uso: bash start-catcher.sh CARTELLA_EB_ANALYSIS [PORTA, default 8081]" >&2
  exit 2
fi
project="$1"
port="${2:-8081}"
case "$port" in ''|*[!0-9]*) echo "La porta deve essere numerica." >&2; exit 2;; esac
if ((port < 1024 || port > 65535)); then
  echo "Scegli una porta fra 1024 e 65535." >&2; exit 2
fi
image="${CATCHER_IMAGE:-hedgelab/catcherv4-hedgelab_1.0:base}"
name="${CATCHER_CONTAINER:-catcher-rstudio-v4}"
platform="${DOCKER_PLATFORM:-linux/amd64}"
docker info >/dev/null
docker pull --platform "$platform" "$image"
if docker container inspect "$name" >/dev/null 2>&1; then
  echo "Il container $name esiste già. Per fermarlo: docker stop $name" >&2
  exit 1
fi
password="${RSTUDIO_PASSWORD:-}"
if [[ -z "$password" ]]; then
  if [[ ! -t 0 ]]; then
    echo "Imposta RSTUDIO_PASSWORD oppure avvia da un terminale interattivo." >&2
    exit 2
  fi
  read -r -s -p "Scegli la password di RStudio: " password
  printf '\n'
fi
if [[ -z "$password" ]]; then echo "La password non può essere vuota." >&2; exit 2; fi
mkdir -p "$project/scripts" "$project/scratch" "$project/Output"
project="$(cd "$project" && pwd -P)"
user_args=()
if [[ "$(uname -s)" == Linux && "$(id -u)" != 0 ]]; then
  user_args=(-e "USERID=$(id -u)" -e "GROUPID=$(id -g)")
fi
export PASSWORD="$password"
docker run -d --rm --platform "$platform" --name "$name" \
  --shm-size=1g -p "127.0.0.1:$port:8787" -e PASSWORD \
  ${user_args[@]+"${user_args[@]}"} \
  --mount "type=bind,src=$project,dst=/home/rstudio/project" \
  "$image" >/dev/null
unset PASSWORD password
echo "RStudio: http://localhost:$port"
echo "Utente: rstudio"
echo "Esperimento: $project"
echo "Per fermare: docker stop $name"