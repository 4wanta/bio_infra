#!/bin/bash
set -euo pipefail

REGION="${AWS_REGION:-ap-northeast-1}"
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
REGISTRY="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

login_ecr() {
  aws ecr get-login-password --region "${REGION}" \
    | docker login --username AWS --password-stdin "${REGISTRY}"
}

ensure_repo() {
  local name="$1"
  aws ecr create-repository --repository-name "${name}" --region "${REGION}" 2>/dev/null \
    || echo "Repository ${name} already exists (or create skipped)."
}

build_push() {
  local repo="$1"
  local subdir="$2"
  local context="${ROOT_DIR}/docker/${subdir}"
  ensure_repo "${repo}"
  docker build -t "${repo}:latest" "${context}"
  docker tag "${repo}:latest" "${REGISTRY}/${repo}:latest"
  docker push "${REGISTRY}/${repo}:latest"
  echo "Pushed ${REGISTRY}/${repo}:latest"
}

echo "Logging in to ECR..."
login_ecr

build_push bio-ngs-qc qc
build_push bio-ngs-mapping mapping
build_push bio-ngs-count count

echo "Done. Update nextflow.config awsbatch profile if account/region differs."
