{
  description = "rue-lib: Python library for the Rapid Urbanisation Explorer project";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    timvim = {
      url = "github:timlinux/timvim";
      flake = false;
    };
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      timvim,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = {
            allowUnfree = true;
          };
        };

        # Python environment with all dependencies
        pythonEnv = pkgs.python313.withPackages (
          ps: with ps; [
            # Core dependencies
            gdal
            geopandas
            rich
            shapely
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
          ]
        );

        # Neovim with timvim configuration
        neovimWithTimvim = pkgs.neovim.override {
          configure = {
            customRC = ''
              " Load timvim configuration
              set runtimepath^=${timvim}
              set runtimepath+=${timvim}/after
              lua << EOF
              -- Set up timvim if init.lua exists
              local timvim_init = "${timvim}/init.lua"
              if vim.fn.filereadable(timvim_init) == 1 then
                dofile(timvim_init)
              end
              EOF
            '';
            packages.myPlugins = with pkgs.vimPlugins; {
              start = [
                # Essential plugins
                vim-nix
                nvim-lspconfig
                nvim-treesitter
              ];
            };
          };
        };

      in
      {
        # Development shell
        devShells.default = pkgs.mkShell {
          buildInputs = [
            pythonEnv
            neovimWithTimvim

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
            export PYTHONPATH="$pythonWithPackages/lib/python*/site-packages:${self}/src:$PYTHONPATH"

            # Set GDAL environment
            export GDAL_DATA="${pkgs.gdal}/share/gdal"
            export PROJ_LIB="${pkgs.proj}/share/proj"

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
            echo -e "   $GRAY▶$RESET  $CYAN./scripts/vscode.sh$RESET  - VSCode preconfigured for python dev"
            echo -e "   $GRAY▶$RESET  $CYAN./scripts/checks.sh$RESET  - Run pre-commit checks"
            echo -e "   $GRAY▶$RESET  $CYAN./scripts/clean.sh$RESET  - Cleanup dev dolder o "
            echo -e "   $GRAY▶$RESET  $CYAN nix flake show$RESET    - Show available configurations"
            echo -e "   $GRAY▶$RESET  $CYAN nix flake check$RESET   - Run all flake checks"
            echo -e "   $GRAY▶$RESET  $CYAN nix run .#test $RESET   - Run all tests"
            echo -e "   $GRAY▶$RESET  $CYAN nix run .#test-cov $RESET   - Run all tests with coverage"
            echo -e "$RESET$ORANGE \n__________________________________________________________________\n"
            echo "To run QGIS with your profile, use one of these commands:"
            echo -e "$RESET$ORANGE \n__________________________________________________________________\n"
            echo ""
            echo "  scripts/start_qgis.sh"
            echo "  scripts/start_qgis_ltr.sh"
            echo "  scripts/start_qgis_master.sh"
            echo ""

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
            rich
            shapely
          ];

          nativeBuildInputs = with pkgs.python313Packages; [
            setuptools
            wheel
          ];

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
            program = toString (
              pkgs.writeShellScript "run-tests" ''
                export PYTHONPATH="${self}/src:$PYTHONPATH"
                export GDAL_DATA="${pkgs.gdal}/share/gdal"
                export PROJ_LIB="${pkgs.proj}/share/proj"
                cd ${self}
                ${pythonEnv}/bin/pytest tests/ -v
              ''
            );
          };

          test-cov = {
            type = "app";
            program = toString (
              pkgs.writeShellScript "run-tests-coverage" ''
                export PYTHONPATH="${self}/src:$PYTHONPATH"
                export GDAL_DATA="${pkgs.gdal}/share/gdal"
                export PROJ_LIB="${pkgs.proj}/share/proj"
                cd ${self}
                ${pythonEnv}/bin/pytest tests/ --cov=rue_lib --cov-report=html --cov-report=term
              ''
            );
          };
        };
      }
    );
}
