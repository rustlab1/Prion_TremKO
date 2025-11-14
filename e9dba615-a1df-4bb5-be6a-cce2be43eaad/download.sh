#!/bin/bash
set -e

echo "Downloading input data into $(pwd)..."

echo "Downloading data.csv..."
wget https://raw.githubusercontent.com/rustlab1/Prion_TremKO/main/e9dba615-a1df-4bb5-be6a-cce2be43eaad/input_data/data.csv

echo "Downloading paper.pdf..."
wget https://raw.githubusercontent.com/rustlab1/Prion_TremKO/main/e9dba615-a1df-4bb5-be6a-cce2be43eaad/input_data/paper.pdf

echo "Download complete. Input data folder is ready."
