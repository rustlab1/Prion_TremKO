#!/bin/bash
set -e

echo "Downloading input data into $(pwd)..."

echo "Downloading data.csv..."
wget https://raw.githubusercontent.com/rustlab1/PER_HYP/main/85bad014-dee0-45b0-87de-6d71e4bf14b2/Input/data.csv

echo "Downloading sxaf055.pdf..."
wget https://raw.githubusercontent.com/rustlab1/PER_HYP/main/85bad014-dee0-45b0-87de-6d71e4bf14b2/Input/sxaf055.pdf

echo "Download complete. Input data folder is ready."

