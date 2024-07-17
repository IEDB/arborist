#!/bin/sh

python3 -m venv _venv
. _venv/bin/activate
sudo chown -R $USER _venv
pip install --upgrade pip
pip install -r requirements.txt