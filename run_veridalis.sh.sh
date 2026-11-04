#!/usr/bin/env bash
# Veridalis launcher for macOS and Linux.
#   bash run_veridalis.sh           start the app
#   bash run_veridalis.sh update    reinstall the packages, then start the app
set -u
cd "$(dirname "$0")" || exit 1

if [ ! -f app.py ] || [ ! -f requirements.txt ]; then
  echo "Could not find app.py next to this file."
  echo "If you downloaded a ZIP, extract it and run this file from the extracted folder."
  exit 1
fi

FORCE=0
[ "${1:-}" = "update" ] && FORCE=1

# ---- find a usable Python, version 3.11 or newer ----
PY=""
for candidate in python3.13 python3.12 python3.11 python3.14 python3 python; do
  if command -v "$candidate" >/dev/null 2>&1 &&
     "$candidate" -c 'import sys; sys.exit(sys.version_info < (3, 11))' >/dev/null 2>&1; then
    PY="$candidate"
    break
  fi
done
if [ -z "$PY" ]; then
  echo "Python 3.11 or newer was not found."
  echo "Install it from https://www.python.org/downloads/ and run this file again."
  exit 1
fi

VPY=".venv/bin/python"

# ---- remove a broken environment, then create one if needed ----
if [ -d .venv ] && ! "$VPY" -c 'import sys' >/dev/null 2>&1; then
  rm -rf .venv
fi
if [ ! -x "$VPY" ]; then
  echo "Setting up Veridalis for the first time. This can take a few minutes."
  if ! "$PY" -m venv .venv; then
    rm -rf .venv
    echo "Could not create the Python environment."
    echo "On Ubuntu or Debian, run: sudo apt install python3-venv"
    echo "Then run this file again."
    exit 1
  fi
fi

# ---- install packages when missing, or when "update" was given ----
if [ "$FORCE" = "1" ] || ! "$VPY" -c 'import streamlit, folium, pandas, numpy' >/dev/null 2>&1; then
  echo "Installing the required packages. This needs an internet connection."
  if ! "$VPY" -m pip install --disable-pip-version-check -r requirements.txt; then
    echo
    echo "Package installation failed."
    echo "Check your internet connection and run this file again."
    echo "If the messages above mention building pandas, numpy or pyarrow, your Python version is probably too new."
    echo "In that case install Python 3.12, delete the .venv folder and run this file again."
    exit 1
  fi
fi

# ---- skip Streamlit's first-run email prompt ----
if [ ! -f "$HOME/.streamlit/credentials.toml" ]; then
  mkdir -p "$HOME/.streamlit"
  printf '[general]\nemail = ""\n' > "$HOME/.streamlit/credentials.toml"
fi

echo
echo "Starting Veridalis. Your browser will open in a moment."
echo "Keep this terminal open while you use the app. Press Ctrl+C to stop it."
echo
exec "$VPY" -m streamlit run app.py --browser.gatherUsageStats false
