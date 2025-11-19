#!/usr/bin/env bash
echo "🪛 Running QGIS with the RUE profile:"
echo "--------------------------------"
echo "Do you want to enable debug mode?"
choice=$(gum choose "🪲 Yes" "🐞 No")
case $choice in
"🪲 Yes") developer_mode=1 ;;
"🐞 No") developer_mode=0 ;;
esac
echo "Do you want to enable experimental features?"
choice=$(gum choose "🪲 Yes" "🐞 No")
case $choice in
"🪲 Yes") RUE_EXPERIMENTAL=1 ;;
"🐞 No") RUE_EXPERIMENTAL=0 ;;
esac

# Running on local used to skip tests that will not work in a local dev env
RUE_LOG=$HOME/RUE.log
RUE_TEST_DIR="$(pwd)/test" # Set test directory relative to project root
rm -f "$RUE_LOG"
#nix-shell -p \
#  This is the old way using default nix packages with overrides
#  'qgis.override { extraPythonPackages = (ps: [ ps.pyqtwebengine ps.jsonschema ps.debugpy ps.future ps.psutil ]);}' \
#  --command "RUE_LOG=${GEEST_LOG} GEEST_DEBUG=${developer_mode} RUNNING_ON_LOCAL=1 qgis --profile RUE"

# This is the new way, using Ivan Mincis nix spatial project and a flake
# see flake.nix for implementation details
RUE_LOG=${GEEST_LOG} \
  RUE_DEBUG=${developer_mode} \
  RUE_EXPERIMENTAL=${GEEST_EXPERIMENTAL} \
  RUE_TEST_DIR=${GEEST_TEST_DIR} \
  RUNNING_ON_LOCAL=1 \
  nix run github:qgis/QGIS#qgis -- --profile RUE2
#nix run .#default -- --profile RUE2
