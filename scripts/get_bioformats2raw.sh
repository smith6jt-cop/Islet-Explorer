#!/usr/bin/env bash
set -euo pipefail
# Download bioformats2raw and raw2ometiff JARs locally for conversions without Docker.

BF2RAW_VERSION="0.8.0"
RAW2OME_VERSION="0.8.0"
DEST_DIR="$(dirname "$0")/tools"
mkdir -p "$DEST_DIR"

echo "Downloading bioformats2raw ${BF2RAW_VERSION} ..."
curl -L -o "$DEST_DIR/bioformats2raw-${BF2RAW_VERSION}.jar" \
  "https://github.com/glencoesoftware/bioformats2raw/releases/download/v${BF2RAW_VERSION}/bioformats2raw-${BF2RAW_VERSION}.jar"

echo "Downloading raw2ometiff ${RAW2OME_VERSION} ..."
curl -L -o "$DEST_DIR/raw2ometiff-${RAW2OME_VERSION}.jar" \
  "https://github.com/glencoesoftware/raw2ometiff/releases/download/v${RAW2OME_VERSION}/raw2ometiff-${RAW2OME_VERSION}.jar"

echo "Done. JARs saved under: $DEST_DIR"
echo "Ensure Java (JRE) is installed, e.g., sudo apt-get install -y openjdk-17-jre"

