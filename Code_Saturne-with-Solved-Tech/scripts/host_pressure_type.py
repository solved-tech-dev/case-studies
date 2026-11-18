"""Identify the host adapter type that selects host ParCSR."""
import re
from pathlib import Path
import sys

source = Path(sys.argv[1]).read_text()
library_strings = set(Path(sys.argv[2]).read_text().splitlines())
pattern = (
    r'else if \(strcmp\(cs_sles_get_type\(sc\), "(cs_sles_\w+_t)"\) == 0\) '
    r'\{\s*strncpy\(external_type, "HYPRE_ParCSR", 31\);\s*need_external = true;'
)
types = set(re.findall(pattern, source))
if len(types) != 1:
    raise SystemExit(f"Expected one host pressure adapter type in {sys.argv[1]}, found {len(types)}")
(host_type,) = types
if host_type not in library_strings:
    raise SystemExit("The host source registration type does not match installed libsaturne")
print(host_type)
