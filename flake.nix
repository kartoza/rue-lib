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

  outputs = { self, nixpkgs, flake-utils, timvim }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = {
            allowUnfree = true;
          };
        };

        # Python environment with all dependencies
        pythonEnv = pkgs.python311.withPackages (ps: with ps; [
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
        ]);

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
            
            # Additional tools
            pkgs.gdal
            pkgs.geos
            pkgs.proj
            pkgs.git
            pkgs.pre-commit
            pkgs.cspell
            
            # Nix tools
            pkgs.nixpkgs-fmt
          ];

          shellHook = ''
            echo "=================================="
            echo "rue-lib development environment"
            echo "=================================="
            echo ""
            echo "Python: $(python --version)"
            echo "GDAL: $(gdal-config --version)"
            echo "Neovim: $(nvim --version | head -n 1)"
            echo ""
            echo "Available commands:"
            echo "  nix run .#test       - Run Python tests"
            echo "  nix run .#test-cov   - Run tests with coverage"
            echo "  pytest               - Run tests directly"
            echo "  pre-commit install   - Install pre-commit hooks"
            echo "  pre-commit run --all-files - Run all pre-commit checks"
            echo "  nvim                 - Open Neovim with timvim config"
            echo ""
            
            # Set up Python path
            export PYTHONPATH="${self}/src:$PYTHONPATH"
            
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
          '';
        };

        # Package definition
        packages.default = pkgs.python311Packages.buildPythonPackage {
          pname = "rue-lib";
          version = "0.1.0";
          src = ./.;
          format = "pyproject";

          propagatedBuildInputs = with pkgs.python311Packages; [
            gdal
            geopandas
            rich
            shapely
          ];

          nativeBuildInputs = with pkgs.python311Packages; [
            setuptools
            wheel
          ];

          checkInputs = with pkgs.python311Packages; [
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
              cd ${self}
              ${pythonEnv}/bin/pytest tests/ --cov=rue_lib --cov-report=html --cov-report=term
            '');
          };
        };
      }
    );
}
