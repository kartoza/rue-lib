#!/bin/bash

echo "Starting GeoPackage Viewer..."
echo ""
echo "Make sure you have installed dependencies:"
echo "  pip install -r requirements.txt"
echo ""
echo "The viewer will be available at:"
echo "  http://localhost:5555"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

cd "$(dirname "$0")"
python server.py
