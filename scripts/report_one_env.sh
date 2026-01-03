#!/usr/bin/env bash
set -euo pipefail

# Report reproducibility artifacts for ONE already-installed micromamba env described by a YAML file.
#
# Usage:
#   bash scripts/report_one_env.sh path/to/environment.yml
#
# Outputs:
#   repro/platform.txt
#   repro/envs/<stem>-conda-env-export.yml     (fully pinned env export with versions/builds)
#   repro/envs/<stem>-conda-list.txt           (package list w/ versions)
#   repro/envs/<stem>-pip-freeze.txt           (if pip present)
#   repro/README-deps.md                       (table; creates if missing, else appends row)

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${REPO_ROOT}/repro"
ENV_OUT_DIR="${OUT_DIR}/envs"
README_SNIPPET="${OUT_DIR}/README-deps.md"
PLATFORM_TXT="${OUT_DIR}/platform.txt"

MM="micromamba"

if [[ $# -ne 1 ]]; then
  echo "[ERROR] Provide exactly one argument: path/to/environment.yml"
  exit 1
fi

REL_PATH="$1"
if [[ "${REL_PATH}" = /* ]]; then
  ABS_PATH="${REL_PATH}"
  DISPLAY_PATH="${REL_PATH}"
else
  ABS_PATH="${REPO_ROOT}/${REL_PATH#./}"
  DISPLAY_PATH="./${REL_PATH#./}"
fi

if [[ ! -f "${ABS_PATH}" ]]; then
  echo "[ERROR] File not found: ${ABS_PATH}"
  exit 1
fi

if ! command -v "$MM" >/dev/null 2>&1; then
  echo "[ERROR] micromamba not found on PATH."
  exit 1
fi

mkdir -p "${ENV_OUT_DIR}"

echo "[INFO] Repo root: ${REPO_ROOT}"
echo "[INFO] Using: $MM ($(command -v "$MM"))"
"$MM" --version || true
echo "[INFO] Env YAML: ${DISPLAY_PATH}"

extract_env_name() {
  local f="$1"
  local line
  line="$(grep -m1 -E '^[[:space:]]*name:[[:space:]]*' "$f" || true)"
  [[ -z "${line}" ]] && return 1
  echo "${line}" | sed -E 's/^[[:space:]]*name:[[:space:]]*//; s/[[:space:]]+#.*$//; s/[[:space:]]+$//'
}

path_to_stem() {
  local p="$1"
  p="${p#./}"
  echo "${p}" | sed 's/\//__/g'
}

env_name="$(extract_env_name "${ABS_PATH}" || true)"
if [[ -z "${env_name}" ]]; then
  echo "[ERROR] Could not parse 'name:' from ${DISPLAY_PATH}. Add a 'name: <env>' line."
  exit 1
fi
echo "[INFO] Parsed env name: ${env_name}"

if ! "$MM" env list | awk '{print $1}' | grep -qx "${env_name}"; then
  echo "[ERROR] Env '${env_name}' is not installed (micromamba env list did not find it)."
  echo "       Install it manually first, e.g.:"
  echo "       micromamba create -n ${env_name} -f ${DISPLAY_PATH} -y"
  exit 1
fi
echo "[INFO] Env exists; proceeding with exports (no install/update)."

# Platform capture (write if missing)
if [[ ! -f "${PLATFORM_TXT}" ]]; then
  {
    echo "## Platform"
    echo ""
    echo "### Date"
    date -Is
    echo ""
    echo "### OS release (/etc/os-release)"
    [[ -f /etc/os-release ]] && cat /etc/os-release || echo "N/A"
    echo ""
    echo "### Kernel (uname -a)"
    uname -a
    echo ""
    echo "### libc (ldd --version)"
    command -v ldd >/dev/null 2>&1 && (ldd --version 2>&1 | head -n 2) || echo "N/A"
    echo ""
    echo "### CPU (lscpu)"
    command -v lscpu >/dev/null 2>&1 && lscpu || echo "N/A"
    echo ""
    echo "### Memory (free -h)"
    command -v free >/dev/null 2>&1 && free -h || echo "N/A"
    echo ""
    echo "### GPU (nvidia-smi)"
    command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi || echo "N/A"
  } > "${PLATFORM_TXT}"
  echo "[INFO] Wrote platform info: ${PLATFORM_TXT}"
else
  echo "[INFO] Platform info already exists: ${PLATFORM_TXT} (leaving unchanged)"
fi

stem="$(path_to_stem "${DISPLAY_PATH}")"
out_export="${ENV_OUT_DIR}/${stem}-conda-env-export.yml"
out_list="${ENV_OUT_DIR}/${stem}-conda-list.txt"
out_lock="${ENV_OUT_DIR}/${stem}-conda-explicit.lock"
out_pip="${ENV_OUT_DIR}/${stem}-pip-freeze.txt"

echo "[INFO] Writing: ${out_export}"
"$MM" env export -n "${env_name}" > "${out_export}"

echo "[INFO] Writing: ${out_list}"
"$MM" list -n "${env_name}" > "${out_list}"

echo "[INFO] Writing: ${out_lock}"
"$MM" env export -n "${env_name}" --explicit > "${out_lock}"

if "$MM" run -n "${env_name}" python -c "import importlib.util,sys; sys.exit(0 if importlib.util.find_spec('pip') else 1)" >/dev/null 2>&1; then
  echo "[INFO] Writing: ${out_pip}"
  "$MM" run -n "${env_name}" python -m pip freeze \
  | grep -v '@ file://' \
  > "${out_pip}"
else
  rm -f "${out_pip}" || true
  echo "[INFO] pip not present; no pip freeze written."
fi

# README snippet
if [[ ! -f "${README_SNIPPET}" ]]; then
  cat > "${README_SNIPPET}" <<'MD'
## Reproducibility: software + OS

This repo uses **micromamba/conda** environments. For each environment YAML in the repository, we record:

- Fully pinned export (includes versions and build strings): `*-conda-env-export.yml`
- Package list (includes versions): `*-conda-list.txt`
- Explicit lockfile (fully pinned, rebuildable): `*-conda-explicit.lock`
- pip freeze (if pip is present): `*-pip-freeze.txt`

### Operating system / platform

- Platform details: `repro/platform.txt`

### Environments discovered in subdirectories

| Env YAML path | Conda env name | Exports |
|---|---:|---|
MD
fi

exports="repro/envs/${stem}-conda-env-export.yml<br>repro/envs/${stem}-conda-list.txt<br>repro/envs/${stem}-conda-explicit.lock"
if [[ -f "${out_pip}" ]]; then
  exports="${exports}<br>repro/envs/${stem}-pip-freeze.txt"
fi

row="| \`${DISPLAY_PATH}\` | \`${env_name}\` | ${exports} |"

if grep -Fq "| \`${DISPLAY_PATH}\` |" "${README_SNIPPET}"; then
  echo "[INFO] README already has an entry for ${DISPLAY_PATH}; leaving README unchanged."
else
  echo "${row}" >> "${README_SNIPPET}"
  echo "[INFO] Appended README row to: ${README_SNIPPET}"
fi

echo "[INFO] Done."
