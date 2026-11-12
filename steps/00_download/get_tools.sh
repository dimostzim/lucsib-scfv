#!/bin/bash
# Download Mikolaj Kocikowski's tools at pinned commits, build the Docker images,
# stage the dog germline references and create the R and Python conda environments.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

usage() {
  cat <<'EOF'
Usage:
  bash steps/00_download/get_tools.sh --tools-dir DIR --seq2scfv-url URL --seq2scfv-commit SHA
    --dog-ref-url URL --dog-ref-commit SHA --qc-image NAME --seq2scfv-image NAME
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tools-dir) TOOLS_DIR="$2"; shift 2 ;;
    --seq2scfv-url) SEQ2SCFV_URL="$2"; shift 2 ;;
    --seq2scfv-commit) SEQ2SCFV_COMMIT="$2"; shift 2 ;;
    --dog-ref-url) DOG_REF_URL="$2"; shift 2 ;;
    --dog-ref-commit) DOG_REF_COMMIT="$2"; shift 2 ;;
    --qc-image) QC_IMAGE="$2"; shift 2 ;;
    --seq2scfv-image) SEQ2SCFV_IMAGE="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p "${TOOLS_DIR}"

checkout() {
  local url="$1" commit="$2" dir="$3"
  if [[ ! -d "${dir}/.git" ]]; then
    git clone "${url}" "${dir}"
  fi
  git -C "${dir}" fetch --quiet origin
  git -C "${dir}" checkout --quiet "${commit}"
  echo "${dir}: $(git -C "${dir}" rev-parse HEAD)"
}

proxy_args=()
for var in HTTP_PROXY HTTPS_PROXY NO_PROXY http_proxy https_proxy no_proxy; do
  if [[ -n "${!var:-}" ]]; then
    proxy_args+=(--build-arg "${var}=${!var}")
  fi
done

echo "[1/5] seq2scfv source"
checkout "${SEQ2SCFV_URL}" "${SEQ2SCFV_COMMIT}" "${TOOLS_DIR}/seq2scfv"

echo "[2/5] seq2scfv Docker image"
docker build --network=host -t "${SEQ2SCFV_IMAGE}" \
  -f "${REPO_DIR}/envs/seq2scfv.Dockerfile" \
  --build-arg SEQ2SCFV_REPO="${SEQ2SCFV_URL}" \
  --build-arg SEQ2SCFV_COMMIT="${SEQ2SCFV_COMMIT}" \
  "${proxy_args[@]}" \
  "${TOOLS_DIR}/seq2scfv/Docker"

echo "[3/5] QC Docker image"
docker build --network=host -t "${QC_IMAGE}" \
  -f "${REPO_DIR}/envs/qc.Dockerfile" \
  "${proxy_args[@]}" \
  "${REPO_DIR}/envs"

echo "[4/5] Dog germline references"
checkout "${DOG_REF_URL}" "${DOG_REF_COMMIT}" "${TOOLS_DIR}/covid-scfv-publication"
SRC="${TOOLS_DIR}/covid-scfv-publication/for_igblast"
OUT_DIR="${TOOLS_DIR}/germlines/dog"
rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}/databases" "${OUT_DIR}/optional_file" "${OUT_DIR}/internal_data"
cp -a "${SRC}/databases/dog" "${OUT_DIR}/databases/"
cp -a "${SRC}/optional_file/dog_gl.aux" "${OUT_DIR}/optional_file/"
cp -a "${SRC}/internal_data/dog" "${OUT_DIR}/internal_data/"

echo "[5/5] Conda environments"
for env in r python; do
  env_dir="${TOOLS_DIR}/envs/${env}"
  yml_sha="$(sha256sum "${REPO_DIR}/envs/${env}.yml" | cut -d' ' -f1)"
  # recreate the environment when its yml changed since it was built
  if [[ -d "${env_dir}" && "$(cat "${env_dir}/.envs_yml_sha256" 2>/dev/null)" != "${yml_sha}" ]]; then
    rm -rf "${env_dir}"
  fi
  if [[ ! -d "${env_dir}" ]]; then
    conda env create --quiet -p "${env_dir}" -f "${REPO_DIR}/envs/${env}.yml"
    echo "${yml_sha}" > "${env_dir}/.envs_yml_sha256"
  fi
  conda list -p "${env_dir}" > "${TOOLS_DIR}/envs/${env}.packages.txt"
done

# record the package versions inside both images
mkdir -p "${TOOLS_DIR}/images"
docker run --rm "${QC_IMAGE}" bash -c 'conda list -n nanopore_qc' > "${TOOLS_DIR}/images/${QC_IMAGE}.packages.txt"
docker run --rm "${SEQ2SCFV_IMAGE}" bash -c 'conda list; pip freeze; igblastn -version' \
  > "${TOOLS_DIR}/images/${SEQ2SCFV_IMAGE}.packages.txt"

echo "Tools ready in ${TOOLS_DIR}"
