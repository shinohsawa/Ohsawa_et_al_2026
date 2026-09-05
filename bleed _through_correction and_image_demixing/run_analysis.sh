#!/bin/bash
# Unix/Linux/Mac shell script to run the demixing analysis

echo "========================================"
echo "Bleedthrough Correction Pipeline"
echo "========================================"
echo

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate demixing

if [ $? -ne 0 ]; then
    echo "ERROR: Could not activate conda environment 'demixing'"
    echo "Please create the environment first:"
    echo "  conda env create -f environment.yml"
    exit 1
fi

echo "Environment activated: demixing"
echo

# Run the analysis
python demix_images.py

if [ $? -ne 0 ]; then
    echo
    echo "ERROR: Analysis failed!"
    echo "Check the error messages above."
    exit 1
fi

echo
echo "========================================"
echo "Analysis complete!"
echo "========================================"
echo "Check the output folders for results."
echo
