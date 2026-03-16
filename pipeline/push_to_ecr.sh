#!/bin/bash
# 設定項目
REPO_NAME="bio-analysis-tools"
REGION="ap-northeast-1"
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
IMAGE_URI="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${REPO_NAME}:latest"

echo "Step 1: Creating ECR repository (if not exists)..."
aws ecr create-repository --repository-name ${REPO_NAME} --region ${REGION} || echo "Repository already exists."

echo "Step 2: Logging in to ECR..."
aws ecr get-login-password --region ${REGION} | docker login --username AWS --password-stdin ${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com

echo "Step 3: Building Docker image..."
docker build -t ${REPO_NAME} .

echo "Step 4: Tagging and Pushing image to ECR..."
docker tag ${REPO_NAME}:latest ${IMAGE_URI}
docker push ${IMAGE_URI}

echo "--------------------------------------------------"
echo "Success! Your image URI is:"
echo ${IMAGE_URI}
echo "--------------------------------------------------"
