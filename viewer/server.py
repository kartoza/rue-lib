#!/usr/bin/env python3
"""Simple GeoPackage viewer server with processing helpers."""

import json
import os
import sys
from datetime import datetime
from pathlib import Path

from flask import Flask, jsonify, request, send_from_directory
from flask_cors import CORS
from osgeo import ogr, osr

ROOT_DIR = Path(__file__).resolve().parent.parent

# Ensure the package source directory is importable when running from viewer/
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from rue_lib.site.runner import SiteConfig, generate_parcels  # noqa: E402
from rue_lib.streets.config import StreetConfig  # noqa: E402
from rue_lib.streets.runner import generate_streets  # noqa: E402

app = Flask(__name__)
CORS(app)

# Path to outputs directory
OUTPUTS_DIR = Path(__file__).parent.parent
RUNS_BASE_DIR = OUTPUTS_DIR / "outputs" / "viewer_runs"


def ensure_feature_collection(data, label):
    """Validate GeoJSON FeatureCollection payloads."""

    if not isinstance(data, dict) or data.get("type") != "FeatureCollection":
        raise ValueError(f"{label} must be a GeoJSON FeatureCollection")

    features = data.get("features")
    if not isinstance(features, list) or not features:
        raise ValueError(f"{label} must contain at least one feature")

    # Filter out empty geometries while validating
    valid_features = []
    for feature in features:
        geom = feature.get("geometry") if isinstance(feature, dict) else None
        if not geom:
            continue
        valid_features.append(
            {
                "type": "Feature",
                "geometry": geom,
                "properties": feature.get("properties", {}) or {},
            }
        )

    if not valid_features:
        raise ValueError(f"{label} features must include valid geometries")

    return {"type": "FeatureCollection", "features": valid_features}


def combine_roads(arterial_fc, secondary_fc):
    """Merge arterial and secondary roads into a single collection with attributes."""

    combined = []

    for feature in arterial_fc["features"]:
        props = feature.get("properties", {}) or {}
        props.setdefault("road_type", "road_art")
        props.setdefault("type", "road_art")
        combined.append({"type": "Feature", "geometry": feature["geometry"], "properties": props})

    for feature in secondary_fc["features"]:
        props = feature.get("properties", {}) or {}
        props.setdefault("road_type", "road_sec")
        props.setdefault("type", "road_sec")
        combined.append({"type": "Feature", "geometry": feature["geometry"], "properties": props})

    if not combined:
        raise ValueError("At least one arterial or secondary road feature is required")

    return {"type": "FeatureCollection", "features": combined}


def save_geojson(data, path):
    """Persist GeoJSON dict to disk."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh)


def relative_to_root(path):
    """Return POSIX-style path relative to the repository root."""

    return str(Path(path).resolve().relative_to(OUTPUTS_DIR))


def run_processing_steps(boundary_fc, arterial_fc, secondary_fc, params=None):
    """Execute Step 1 and Step 2 sequentially using uploaded inputs and parameters."""

    boundary = ensure_feature_collection(boundary_fc, "Boundary")
    arterial = ensure_feature_collection(arterial_fc, "Arterial roads")
    secondary = ensure_feature_collection(secondary_fc, "Secondary roads")

    roads = combine_roads(arterial, secondary)

    timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    run_dir = RUNS_BASE_DIR / timestamp
    run_dir.mkdir(parents=True, exist_ok=True)

    boundary_path = run_dir / "boundary.geojson"
    roads_path = run_dir / "roads.geojson"

    save_geojson(boundary, boundary_path)
    save_geojson(roads, roads_path)

    # Extract parameters with defaults
    params = params or {}

    # Step 2 parameters
    arterial_width = params.get("road_arterial_width_m", 20.0)
    secondary_width = params.get("road_secondary_width_m", 15.0)
    local_width = params.get("road_local_width_m", 12.0)
    arterial_setback = params.get("arterial_setback_depth", 60.0)
    secondary_setback = params.get("secondary_setback_depth", 60.0)
    pref_depth = params.get("off_grid_partitions_preferred_depth", 140.0)
    pref_width = params.get("off_grid_partitions_preferred_width", 140.0)
    grid_depth_arterial = params.get("on_grid_partition_depth_arterial_roads", 40.0)
    grid_depth_secondary = params.get("on_grid_partition_depth_secondary_roads", 30.0)

    step1_dir = run_dir / "step1"
    site_cfg = SiteConfig(
        site_path=str(boundary_path),
        roads_path=str(roads_path),
        output_dir=str(step1_dir),
        geopackage_path=str(step1_dir / "output.gpkg"),
        road_arterial_width_m=arterial_width,
        road_secondary_width_m=secondary_width,
    )
    step1_output = generate_parcels(site_cfg)

    step2_dir = run_dir / "step2"
    street_cfg = StreetConfig(
        parcel_path=str(boundary_path),
        roads_path=str(roads_path),
        output_dir=str(step2_dir),
        road_arterial_width_m=arterial_width,
        road_secondary_width_m=secondary_width,
        geopackage_path=str(step2_dir / "output.gpkg"),
        road_locals_width_m=local_width,
        arterial_setback_depth=arterial_setback,
        secondary_setback_depth=secondary_setback,
        off_grid_partitions_preferred_depth=pref_depth,
        off_grid_partitions_preferred_width=pref_width,
        on_grid_partition_depth_arterial_roads=grid_depth_arterial,
        on_grid_partition_depth_secondary_roads=grid_depth_secondary,
    )
    step2_output = generate_streets(street_cfg)

    return {
        "run_directory": relative_to_root(run_dir),
        "boundary_path": relative_to_root(boundary_path),
        "roads_path": relative_to_root(roads_path),
        "step1_output": relative_to_root(step1_output),
        "step2_output": relative_to_root(step2_output),
    }


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


@app.route("/api/save-inputs", methods=["POST"])
def api_save_inputs():
    """Save boundary and roads to the outputs folder."""

    payload = request.get_json(silent=True)
    if not payload:
        return jsonify({"error": "Invalid or missing JSON payload"}), 400

    try:
        boundary = payload.get("boundary")
        arterial = payload.get("arterial")
        secondary = payload.get("secondary")

        # Create saved_inputs directory
        saved_dir = OUTPUTS_DIR / "saved_inputs"
        saved_dir.mkdir(parents=True, exist_ok=True)

        saved_files = {}

        if boundary:
            boundary_fc = ensure_feature_collection(boundary, "Boundary")
            boundary_path = saved_dir / "boundary.geojson"
            save_geojson(boundary_fc, boundary_path)
            saved_files["boundary"] = relative_to_root(boundary_path)

        if arterial:
            arterial_fc = ensure_feature_collection(arterial, "Arterial roads")
            arterial_path = saved_dir / "arterial_roads.geojson"
            save_geojson(arterial_fc, arterial_path)
            saved_files["arterial"] = relative_to_root(arterial_path)

        if secondary:
            secondary_fc = ensure_feature_collection(secondary, "Secondary roads")
            secondary_path = saved_dir / "secondary_roads.geojson"
            save_geojson(secondary_fc, secondary_path)
            saved_files["secondary"] = relative_to_root(secondary_path)

        return jsonify({"message": "Inputs saved successfully", "files": saved_files})
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    except Exception as exc:
        print(f"Error saving inputs: {exc}")
        return jsonify({"error": "Failed to save inputs"}), 500


@app.route("/api/load-saved-inputs", methods=["GET"])
def api_load_saved_inputs():
    """Load the most recent saved inputs."""

    saved_dir = OUTPUTS_DIR / "saved_inputs"
    if not saved_dir.exists():
        return jsonify({"exists": False})

    result = {"exists": True, "files": {}}

    boundary_path = saved_dir / "boundary.geojson"
    arterial_path = saved_dir / "arterial_roads.geojson"
    secondary_path = saved_dir / "secondary_roads.geojson"

    try:
        if boundary_path.exists():
            with boundary_path.open("r", encoding="utf-8") as f:
                result["files"]["boundary"] = json.load(f)

        if arterial_path.exists():
            with arterial_path.open("r", encoding="utf-8") as f:
                result["files"]["arterial"] = json.load(f)

        if secondary_path.exists():
            with secondary_path.open("r", encoding="utf-8") as f:
                result["files"]["secondary"] = json.load(f)

        return jsonify(result)
    except Exception as exc:
        print(f"Error loading saved inputs: {exc}")
        return jsonify({"error": "Failed to load saved inputs"}), 500


@app.route("/api/run-steps", methods=["POST"])
def api_run_steps():
    """Run Step 1 and Step 2 with boundary and road inputs."""

    payload = request.get_json(silent=True)
    if not payload:
        return jsonify({"error": "Invalid or missing JSON payload"}), 400

    try:
        boundary = payload.get("boundary")
        arterial = payload.get("arterial")
        secondary = payload.get("secondary")
        params = payload.get("params", {})

        if not boundary or not arterial or not secondary:
            return jsonify(
                {"error": "Boundary, arterial roads, and secondary roads are required"}
            ), 400

        result = run_processing_steps(boundary, arterial, secondary, params)
        return jsonify({"message": "Processing complete", **result})
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    except Exception as exc:  # pragma: no cover - catch-all for unexpected issues
        print(f"Error running processing steps: {exc}")
        return jsonify({"error": "Failed to process inputs"}), 500


if __name__ == "__main__":
    print("Starting GeoPackage Viewer...")
    print(f"Serving GeoPackages from: {OUTPUTS_DIR}")
    print("Open http://localhost:5555 in your browser")
    app.run(debug=True, host="0.0.0.0", port=5555)
