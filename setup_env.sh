#!/usr/bin/env bash
# Create an isolated venv for the BioEmu import tooling on the VM.
# Usage:  ./setup_env.sh [/path/to/venv]
# Then:   source <venv>/bin/activate   (or call <venv>/bin/python bioemu_import.py ...)
#
# mdr-process is NOT installed here — it is expected on PATH in exouser's env.
set -euo pipefail

VENV="${1:-$HOME/bioemu-venv}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 -m venv "$VENV"
# shellcheck disable=SC1091
source "$VENV/bin/activate"
pip install --upgrade pip
pip install -r "$HERE/requirements.txt"

echo
echo "venv ready: $VENV"
echo "sanity check:"
python - <<'PY'
import parmed, numpy, sys
print("  parmed", parmed.__version__, "| numpy", numpy.__version__,
      "| python", "%d.%d" % sys.version_info[:2])
PY
echo
echo "next:  $VENV/bin/python $HERE/bioemu_import.py --root /opt/bioemu init"
