@echo off
REM Windows batch script to run the demixing analysis
REM Double-click this file to run the pipeline

echo ========================================
echo Bleedthrough Correction Pipeline
echo ========================================
echo.

REM Activate conda environment
call conda activate demixing
if errorlevel 1 (
    echo ERROR: Could not activate conda environment 'demixing'
    echo Please create the environment first:
    echo   conda env create -f environment.yml
    pause
    exit /b 1
)

echo Environment activated: demixing
echo.

REM Run the analysis
python demix_images.py
if errorlevel 1 (
    echo.
    echo ERROR: Analysis failed!
    echo Check the error messages above.
    pause
    exit /b 1
)

echo.
echo ========================================
echo Analysis complete!
echo ========================================
echo Check the output folders for results.
echo.
pause
