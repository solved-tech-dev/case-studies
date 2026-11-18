"""Run a relocated Code_Saturne installation with a temporary build configuration."""

import configparser
import os
from pathlib import Path
import sys
import tempfile


prefix = Path(os.environ["CS_PREFIX"])
build = Path(os.environ.get("CS_LINGLONGAMG_BUILD", Path(__file__).resolve().parents[1] / "build"))
build.mkdir(parents=True, exist_ok=True)

text = (prefix / "lib/code_saturne_build.cfg").read_text()
old_root = os.environ.get("CS_ORIGINAL_BUILD_ROOT")
if old_root:
    text = text.replace(old_root, os.environ["CS_BUILD_ROOT"])

runtime = os.environ.get("CS_LINGLONGAMG_GCC_RUNTIME")
if runtime:
    # Some compiler runtime packages supply the versioned libraries alone.
    text = text.replace("-lgfortran", "-l:libgfortran.so.5").replace("-lquadmath", "-l:libquadmath.so.0")

parsed = configparser.ConfigParser(interpolation=None)
parsed.read_string(text)
if runtime:
    parsed["compilers"]["ldflags"] = f"-L{runtime} -Wl,-rpath,{runtime} " + parsed["compilers"]["ldflags"]

site = next((prefix / "lib").glob("python*/site-packages"))
sys.path.insert(0, str(site))
from code_saturne.base.cs_package import package
from code_saturne.base.cs_script import master_script

with tempfile.TemporaryDirectory(dir=build) as temp:
    config = Path(temp) / "code_saturne_build.cfg"
    with config.open("w") as output:
        parsed.write(output)
    pkg = package(config_file=str(config), install_prefix=str(prefix))
    raise SystemExit(master_script(sys.argv[1:], pkg).execute())
