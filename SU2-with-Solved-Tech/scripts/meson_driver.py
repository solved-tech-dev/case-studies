#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-or-later
"""Run Meson from the prepared SU2 tree and set generated-state timestamps."""
import os
from pathlib import Path
import sys
import time

source = Path(os.environ["SU2_SOURCE_DIR"])
sys.path.insert(0, str(source / "externals/meson"))
from mesonbuild import coredata, mesonmain

upstream_save = coredata.save


def save_with_client_time(obj, build_dir):
    filename = upstream_save(obj, build_dir)
    now = time.time()
    os.utime(filename, (now, now))
    return filename


coredata.save = save_with_client_time
if __name__ == "__main__":
    sys.exit(mesonmain.main())
