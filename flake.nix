{
  description =
    "rue-lib: Python library for the Rapid Urbanisation Explorer project";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils, }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; };
        };

        # Python environment with all dependencies
        pythonEnv = pkgs.python313.withPackages (ps:
          with ps; [
            # Core dependencies
            gdal
            geopandas
            scipy
            rich
            shapely
            # TUI dependencies
            textual
            matplotlib
            fiona
            pillow
            # Development dependencies
            pytest
            pytest-cov
            ruff
            bandit
            pip
            setuptools
            wheel
            build
            twine
          ]);

      in {
        # Development shell
        devShells.default = pkgs.mkShell {
          buildInputs = [
            pythonEnv
            # Image viewers for TUI
            pkgs.fim # Frame buffer image viewer
            pkgs.feh # X11 image viewer (fallback)

            pkgs.dbeaver-bin # Database viewer
            pkgs.chafa
            # Additional tools
            pkgs.gdal
            pkgs.geos
            pkgs.proj
            pkgs.git
            pkgs.pre-commit
            pkgs.cspell
            pkgs.gum

            # Nix tools
            pkgs.nixpkgs-fmt
          ];

          shellHook = ''
            # Set Python path for rue-lib and rue-tui
            export PYTHONPATH="${self}/src:${self}:$PYTHONPATH"

            # Set GDAL environment
            export GDAL_DATA="${pkgs.gdal}/share/gdal"
            export PROJ_LIB="${pkgs.proj}/share/proj"

            # Create convenient aliases
            alias rue-tui="python -m rue_tui"
            alias run-tui="python scripts/run_tui.py"

            # Ensure pre-commit is installed
            if [ -f .git/hooks/pre-commit ]; then
              echo "Pre-commit hooks already installed"
            else
              echo "Installing pre-commit hooks..."
              pre-commit install
            fi
            echo "Python: $(python --version)"
            echo "GDAL: $(gdal-config --version)"
            echo "Neovim: $(nvim --version | head -n 1)"
            # Colors and styling
            CYAN='\033[38;2;83;161;203m'
            GREEN='\033[92m'
            RED='\033[91m'
            RESET='\033[0m'
            ORANGE='\033[38;2;237;177;72m'
            GRAY='\033[90m'
            # Clear screen and show welcome banner
            clear
            echo -e "$RESET$ORANGE"
            chafa .github/assets/banner.png --size=30x80 --colors=256 | sed 's/^/                  /'
            # Quick tips with icons
            echo -e "$RESET$ORANGE \n__________________________________________________________________\n"
            echo -e "        🌈 Your Dev Environment is prepared."
            echo -e ""
            echo -e "Quick Commands:$RESET"
            echo -e "   $GRAY▶$RESET  $CYAN./scripts/vscode.sh$RESET     - VSCode preconfigured for python dev"
            echo -e "   $GRAY▶$RESET  $CYAN./scripts/checks.sh$RESET     - Run pre-commit checks"
            echo -e "   $GRAY▶$RESET  $CYAN./scripts/clean.sh$RESET      - Cleanup dev folder"
            echo -e "   $GRAY▶$RESET  $CYAN nix flake show$RESET        - Show available configurations"
            echo -e "   $GRAY▶$RESET  $CYAN nix flake check$RESET       - Run all flake checks"
            echo -e "   $GRAY▶$RESET  $CYAN nix run .#test$RESET         - Run all tests"
            echo -e "   $GRAY▶$RESET  $CYAN nix run .#test-cov$RESET     - Run all tests with coverage"
            echo -e "   $GRAY▶$RESET  $CYAN nix run .#run-examples$RESET - Run all examples and open in QGIS"
            echo -e "$RESET$ORANGE \n__________________________________________________________________\n"
            echo "To run QGIS with your profile, use one of these commands:"
            echo -e "$RESET$ORANGE \n__________________________________________________________________\n"
            echo ""
            echo "  scripts/start_qgis.sh"
            echo "  scripts/start_qgis_ltr.sh"
            echo "  scripts/start_qgis_master.sh"
            echo ""
            echo -e "$RESET$ORANGE \n__________________________________________________________________\n"
            echo "To run the RUE text user interface wizard:"
            echo -e "$RESET$ORANGE \n__________________________________________________________________\n"
            echo ""
            echo "  scripts/rue_tui.py"

          '';
        };

        # Package definition
        packages.default = pkgs.python313Packages.buildPythonPackage {
          pname = "rue-lib";
          version = "0.1.0";
          src = ./.;
          format = "pyproject";

          propagatedBuildInputs = with pkgs.python313Packages; [
            gdal
            geopandas
            scipy
            rich
            shapely
          ];

          nativeBuildInputs = with pkgs.python313Packages; [ setuptools wheel ];

          nativeCheckInputs = with pkgs.python313Packages; [
            pytest
            pytest-cov
          ];

          checkPhase = ''
            pytest tests/
          '';
        };

        # App to run tests
        apps = {
          test = {
            type = "app";
            program = toString (pkgs.writeShellScript "run-tests" ''
              export PYTHONPATH="${self}/src:$PYTHONPATH"
              export GDAL_DATA="${pkgs.gdal}/share/gdal"
              export PROJ_LIB="${pkgs.proj}/share/proj"
              cd ${self}
              ${pythonEnv}/bin/pytest tests/ -v
            '');
          };

          test-cov = {
            type = "app";
            program = toString (pkgs.writeShellScript "run-tests-coverage" ''
              export PYTHONPATH="${self}/src:$PYTHONPATH"
              export GDAL_DATA="${pkgs.gdal}/share/gdal"
              export PROJ_LIB="${pkgs.proj}/share/proj"
              # Create a writable temp directory for coverage data
              TEMP_DIR=$(mktemp -d)
              cp -r ${self}/* "$TEMP_DIR/" 2>/dev/null || true
              cd "$TEMP_DIR"
              export COVERAGE_FILE="$TEMP_DIR/.coverage"
              ${pythonEnv}/bin/pytest tests/ --cov=rue_lib --cov-report=html --cov-report=term
              echo "Coverage report written to $TEMP_DIR/htmlcov/"
            '');
          };

          rue-tui = {
            type = "app";
            program = toString (pkgs.writeShellScript "rue-tui" ''
              export PYTHONPATH="${self}/src:${self}:$PYTHONPATH"
              export GDAL_DATA="${pkgs.gdal}/share/gdal"
              export PROJ_LIB="${pkgs.proj}/share/proj"
              cd ${self}

              echo "🏙️ Launching RUE TUI..."
              echo "======================="

              # Create demo data if needed
              if [ ! -f "data/site.geojson" ] || [ ! -f "data/roads.geojson" ]; then
                echo "📁 Creating demo data..."
                mkdir -p data
                ${pythonEnv}/bin/python rue_tui/demo_data.py
              fi

              # Launch the TUI
              ${pythonEnv}/bin/python -m rue_tui
            '');
          };

          run-examples = {
            type = "app";
            program = toString (pkgs.writeShellScript "run-examples" ''
              export PYTHONPATH="${self}/src:$PYTHONPATH"
              export GDAL_DATA="${pkgs.gdal}/share/gdal"
              export PROJ_LIB="${pkgs.proj}/share/proj"
              cd ${self}

              # Create unique writable output directory
              OUTPUT_DIR="/tmp/rue-examples-$(date +%Y%m%d-%H%M%S)"
              mkdir -p "$OUTPUT_DIR"

              echo "🌈 Running RUE-lib examples in sequence..."
              echo "========================================="
              echo "Output directory: $OUTPUT_DIR"
              echo ""

              # Step 1: Generate Parcels
              echo "🏗️  Running Step 1: Generate Parcels"
              echo "-----------------------------------"
              ${pythonEnv}/bin/python examples/step1_generate_parcels.py --output-dir "$OUTPUT_DIR" --geopackage "$OUTPUT_DIR/output.gpkg"

              # Step 2: Generate Streets
              echo ""
              echo "🛣️  Running Step 2: Generate Streets"
              echo "----------------------------------"
              ${pythonEnv}/bin/python examples/step2_generate_streets.py --parcels "$OUTPUT_DIR/parcels.geojson" --output-dir "$OUTPUT_DIR" --geopackage "$OUTPUT_DIR/output.gpkg"

              # Step 3: Generate Clusters
              echo ""
              echo "🏘️  Running Step 3: Generate Clusters"
              echo "------------------------------------"
              ${pythonEnv}/bin/python examples/step3_generate_cluster.py --input "$OUTPUT_DIR/all_grids_merged.geojson" --output-dir "$OUTPUT_DIR" --geopackage "$OUTPUT_DIR/output.gpkg"

              # Open output geopackage in QGIS
              echo ""
              echo "🗺️  Opening output geopackage in QGIS LTR..."
              echo "-------------------------------------------"
              if [ -f "$OUTPUT_DIR/output.gpkg" ]; then
                nix run github:qgis/QGIS#qgis-ltr -- --profile RUE "$OUTPUT_DIR/output.gpkg" &
                echo "✅ QGIS LTR opened with output.gpkg"
              else
                echo "❌ Output geopackage not found at $OUTPUT_DIR/output.gpkg"
                echo "Looking for alternative outputs..."
                find "$OUTPUT_DIR" -name "*.gpkg" -type f | head -5
              fi

              echo ""
              echo "🎉 All examples completed!"
              echo "========================="
              echo "📁 Outputs saved to: $OUTPUT_DIR"
            '');
          };
        };
      });
}
