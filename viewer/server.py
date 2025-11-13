#!/usr/bin/env python3
"""
Simple GeoPackage viewer server
"""

import json
import os
from pathlib import Path

from flask import Flask, jsonify, send_from_directory
from flask_cors import CORS
from osgeo import ogr, osr

app = Flask(__name__)
CORS(app)

# Path to outputs directory
OUTPUTS_DIR = Path(__file__).parent.parent


def get_geopackages():
    """Get all GeoPackage files from outputs directory."""
    gpkg_files = []
    for root, _dirs, files in os.walk(OUTPUTS_DIR):
        for file in files:
            if file.endswith(".gpkg"):
                rel_path = os.path.relpath(os.path.join(root, file), OUTPUTS_DIR)
                gpkg_files.append(rel_path)
    return gpkg_files


def get_layers(gpkg_path):
    """Get all layers from a GeoPackage."""
    full_path = OUTPUTS_DIR / gpkg_path
    if not full_path.exists():
        return []

    ds = ogr.Open(str(full_path))
    if ds is None:
        return []

    layers = []
    for i in range(ds.GetLayerCount()):
        layer = ds.GetLayerByIndex(i)
        layers.append({"name": layer.GetName(), "feature_count": layer.GetFeatureCount()})

    ds = None
    return layers


def layer_to_geojson(gpkg_path, layer_name):
    """Convert a GeoPackage layer to GeoJSON."""
    full_path = OUTPUTS_DIR / gpkg_path
    if not full_path.exists():
        return None

    ds = ogr.Open(str(full_path))
    if ds is None:
        return None

    layer = ds.GetLayerByName(layer_name)
    if layer is None:
        ds = None
        return None

    # Get source CRS
    source_srs = layer.GetSpatialRef()

    # Create target CRS (WGS84 - EPSG:4326)
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(4326)

    # Create coordinate transformation
    transform = None
    if source_srs and not source_srs.IsSame(target_srs):
        transform = osr.CoordinateTransformation(source_srs, target_srs)

    # Create GeoJSON structure
    features = []
    for feature in layer:
        geom = feature.GetGeometryRef()
        if geom:
            # Clone geometry to avoid modifying original
            geom = geom.Clone()

            # Transform to WGS84 if needed
            if transform:
                geom.Transform(transform)

            # Get properties
            properties = {}
            for i in range(feature.GetFieldCount()):
                field_name = feature.GetFieldDefnRef(i).GetName()
                properties[field_name] = feature.GetField(i)

            # Convert geometry to GeoJSON
            geom_json = json.loads(geom.ExportToJson())

            features.append({"type": "Feature", "geometry": geom_json, "properties": properties})

    ds = None

    geojson = {"type": "FeatureCollection", "features": features}

    return geojson


@app.route("/")
def index():
    """Serve the main HTML page."""
    return send_from_directory(".", "index.html")


@app.route("/api/geopackages")
def api_geopackages():
    """Get list of all GeoPackages."""
    gpkgs = get_geopackages()
    return jsonify(gpkgs)


@app.route("/api/layers/<path:gpkg_path>")
def api_layers(gpkg_path):
    """Get list of layers in a GeoPackage."""
    layers = get_layers(gpkg_path)
    return jsonify(layers)


@app.route("/api/layer/<path:gpkg_path>/<layer_name>")
def api_layer_data(gpkg_path, layer_name):
    """Get GeoJSON data for a specific layer."""
    geojson = layer_to_geojson(gpkg_path, layer_name)
    if geojson is None:
        return jsonify({"error": "Layer not found"}), 404
    return jsonify(geojson)


if __name__ == "__main__":
    print("Starting GeoPackage Viewer...")
    print(f"Serving GeoPackages from: {OUTPUTS_DIR}")
    print("Open http://localhost:5555 in your browser")
    app.run(debug=True, host="0.0.0.0", port=5555)
