"""Geospatial functionality examples for rue-lib."""

# Note: This requires GDAL and GeoPandas to be installed
# In Nix environment: nix develop
# Otherwise: pip install gdal geopandas

try:
    from rue_lib.geo import (
        create_sample_geodataframe,
        get_driver_count,
        get_gdal_version,
    )

    # Display GDAL version
    print(f"GDAL Version: {get_gdal_version()}")

    # Display available OGR drivers
    print(f"Available OGR Drivers: {get_driver_count()}")

    # Create and display sample GeoDataFrame
    print("\nSample GeoDataFrame:")
    gdf = create_sample_geodataframe()
    print(gdf)

    print("\nGeometry details:")
    print(gdf.geometry)

except ImportError as e:
    print(f"Error: {e}")
    print("\nTo use geospatial features, install with:")
    print("  pip install rue-lib[geo]")
    print("\nOr use the Nix development environment:")
    print("  nix develop")
