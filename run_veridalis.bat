@echo off
rem Veridalis launcher for Windows.
rem Double-click to start the app. Run "run_veridalis.bat update" to reinstall the packages first.
setlocal EnableExtensions
cd /d "%~dp0"
title Veridalis

if not exist "app.py" goto no_app
if not exist "requirements.txt" goto no_app

set "FORCE=0"
if /i "%~1"=="update" set "FORCE=1"

rem ---- find a usable Python, version 3.11 or newer ----
set "PY="
set "PYCHECK=import sys; sys.exit(sys.version_info < (3, 11))"

where py >nul 2>nul
if errorlevel 1 goto try_python
py -3 -c "%PYCHECK%" >nul 2>nul
if not errorlevel 1 set "PY=py -3"
if defined PY goto have_python

:try_python
python -c "%PYCHECK%" >nul 2>nul
if not errorlevel 1 set "PY=python"
if defined PY goto have_python

for /d %%D in ("%LOCALAPPDATA%\Programs\Python\Python3*" "%ProgramFiles%\Python3*") do (
  if exist "%%~D\python.exe" if not defined PY (
    "%%~D\python.exe" -c "%PYCHECK%" >nul 2>nul
    if not errorlevel 1 set PY="%%~D\python.exe"
  )
)
if defined PY goto have_python
goto no_python

:have_python
set "VPY=.venv\Scripts\python.exe"

rem ---- remove a broken environment, then create one if needed ----
if exist "%VPY%" (
  "%VPY%" -c "import sys" >nul 2>nul
  if errorlevel 1 rmdir /s /q ".venv"
)
if not exist "%VPY%" (
  echo Setting up Veridalis for the first time. This can take a few minutes.
  %PY% -m venv .venv
  if errorlevel 1 goto venv_fail
)

rem ---- install packages when missing, or when "update" was given ----
set "NEEDS=%FORCE%"
"%VPY%" -c "import streamlit, folium, pandas, numpy" >nul 2>nul
if errorlevel 1 set "NEEDS=1"
if "%NEEDS%"=="1" (
  echo Installing the required packages. This needs an internet connection.
  "%VPY%" -m pip install --disable-pip-version-check -r requirements.txt
  if errorlevel 1 goto pip_fail
)

rem ---- skip Streamlit's first-run email prompt ----
if not exist "%USERPROFILE%\.streamlit\credentials.toml" call :write_credentials

echo.
echo Starting Veridalis. Your browser will open in a moment.
echo Keep this window open while you use the app. Press Ctrl+C or close the window to stop it.
echo.
"%VPY%" -m streamlit run app.py --browser.gatherUsageStats false
echo.
echo Veridalis has stopped.
pause
goto end

:write_credentials
mkdir "%USERPROFILE%\.streamlit" >nul 2>nul
> "%USERPROFILE%\.streamlit\credentials.toml" echo [general]
>> "%USERPROFILE%\.streamlit\credentials.toml" echo email = ""
exit /b 0

:no_app
echo Could not find app.py next to this file.
echo If you downloaded a ZIP, right-click it, choose Extract All, and run this file from the extracted folder, not from inside the ZIP.
pause
exit /b 1

:no_python
echo Python 3.11 or newer was not found on this computer.
echo.
echo 1. Install Python from the page that is about to open.
echo 2. During installation, tick the box "Add python.exe to PATH".
echo 3. Run this file again.
start "" "https://www.python.org/downloads/"
pause
exit /b 1

:venv_fail
echo Could not create the Python environment.
echo Try installing Python 3.12 from python.org, then run this file again.
pause
exit /b 1

:pip_fail
echo Package installation failed.
echo Check your internet connection and run this file again.
echo If the messages above mention building pandas, numpy or pyarrow, your Python version is probably too new for them.
echo In that case install Python 3.12, delete the folder named .venv next to this file, and run this file again.
pause
exit /b 1

:end
endlocal
