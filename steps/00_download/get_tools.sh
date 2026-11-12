#!/bin/bash
set -euo pipefail

REPO_URL=${SEQ2SCFV_REPO_URL:-"git@github.com:dimostzim/seq2scfv-unofficial-updated.git"}
REPO_REF=${SEQ2SCFV_REPO_REF:-"main"}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_DIR="${ROOT_DIR}/seq2scfv"
IMAGE_NAME=${SEQ2SCFV_IMAGE_NAME:-"seq2scfv"}
HTTP_PROXY=${HTTP_PROXY:-""}
HTTPS_PROXY=${HTTPS_PROXY:-""}
NO_PROXY=${NO_PROXY:-""}
DOG_REF_REPO_URL=${DOG_REF_REPO_URL:-"https://github.com/MikolajKocikowski/covid-scfv-publication.git"}
DOG_REF_REPO_REF=${DOG_REF_REPO_REF:-"main"}
DOG_REF_DIR=${DOG_REF_DIR:-"${ROOT_DIR}/covid-scfv-publication"}
OUT_DIR=${OUT_DIR:-"${ROOT_DIR}/germlines/dog"}

usage() {
  cat <<'EOF'
Usage:
  bash scfv_analysis/setup/setup_annotation.sh

Environment overrides:
  SEQ2SCFV_REPO_URL
  SEQ2SCFV_REPO_REF
  SEQ2SCFV_IMAGE_NAME
  DOG_REF_REPO_URL
  DOG_REF_REPO_REF
  DOG_REF_DIR
  OUT_DIR
  HTTP_PROXY, HTTPS_PROXY, NO_PROXY
EOF
}

if [[ $# -gt 0 ]]; then
  case "$1" in
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
fi

if [ -d "${REPO_DIR}/.git" ]; then
  git -C "${REPO_DIR}" fetch --all --prune
  git -C "${REPO_DIR}" checkout "${REPO_REF}"
  git -C "${REPO_DIR}" pull --ff-only || true
else
  git clone "${REPO_URL}" "${REPO_DIR}"
  git -C "${REPO_DIR}" checkout "${REPO_REF}"
fi

# Build the Docker image using the repo Dockerfile, but point it to the fork.
docker build --network=host -t "${IMAGE_NAME}" \
  -f "${REPO_DIR}/Docker/Dockerfile" \
  --build-arg SEQ2SCFV_REPO="${REPO_URL}" \
  --build-arg SEQ2SCFV_REF="${REPO_REF}" \
  ${HTTP_PROXY:+--build-arg HTTP_PROXY="${HTTP_PROXY}"} \
  ${HTTPS_PROXY:+--build-arg HTTPS_PROXY="${HTTPS_PROXY}"} \
  ${NO_PROXY:+--build-arg NO_PROXY="${NO_PROXY}"} \
  "${REPO_DIR}/Docker"

# Clone/update dog reference repo and stage files for pipeline use.
if [ -d "${DOG_REF_DIR}/.git" ]; then
  git -C "${DOG_REF_DIR}" fetch --all --prune
  git -C "${DOG_REF_DIR}" checkout "${DOG_REF_REPO_REF}"
  git -C "${DOG_REF_DIR}" pull --ff-only || true
else
  git clone "${DOG_REF_REPO_URL}" "${DOG_REF_DIR}"
  git -C "${DOG_REF_DIR}" checkout "${DOG_REF_REPO_REF}"
fi

SRC="${DOG_REF_DIR}/for_igblast"
if [ ! -d "${SRC}/databases/dog" ] || [ ! -f "${SRC}/optional_file/dog_gl.aux" ] || [ ! -d "${SRC}/internal_data/dog" ]; then
  echo "ERROR: Missing expected dog reference files in ${SRC}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"
rm -rf "${OUT_DIR}/databases" "${OUT_DIR}/optional_file" "${OUT_DIR}/internal_data"
mkdir -p "${OUT_DIR}/databases" "${OUT_DIR}/optional_file" "${OUT_DIR}/internal_data"
cp -a "${SRC}/databases/dog" "${OUT_DIR}/databases/"
cp -a "${SRC}/optional_file/dog_gl.aux" "${OUT_DIR}/optional_file/"
cp -a "${SRC}/internal_data/dog" "${OUT_DIR}/internal_data/"

echo "scFv analysis setup complete: image '${IMAGE_NAME}', dog germlines '${OUT_DIR}'"
