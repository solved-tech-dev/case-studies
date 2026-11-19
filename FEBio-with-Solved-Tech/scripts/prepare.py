"""Build the article's three FEBio inputs from a pinned TestSuite checkout."""

import argparse
import hashlib
import runpy
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "bp10": (
        "bp10-12",
        (12,),
        "6d4c09fa6178fc7abc0c81ff59adcba5fb478670c1532fd48d074e4e56705aef",
        "1d264456635fbe7bf0a7c0bd3b112718df08f16b3d0447e8672556dd779d4fc4",
    ),
    "tu01": (
        "tu01-8",
        (8,),
        "85d8cd2d09ad9d5568ea9131e59a6e6562e7090ab859e56f2e22165aaf65fef5",
        "9b807de5023489b883d94a5b6c69763970dd5ab3fc9cdb03627c7e3e7bc55914",
    ),
    "mi06": (
        "mi06-refined",
        (),
        "3e63e170ca8cccd0f4c9f47e10a8ef7719f5ea59b586dd34c74a201efb593f8d",
        "78ea699c0b79754d6200e463b1e9f40bbab29740bd4eb29115a3594abdc4945f",
    ),
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(case, suite):
    name, factors, source_hash, deck_hash = CASES[case]
    source = suite / "Verify" / (case + ".feb")
    folder = ROOT / "cases" / case
    if digest(source) != source_hash:
        raise ValueError(f"Unexpected upstream input: {source}")
    output = folder / (name + ".feb")
    if output.exists():
        if digest(output) == deck_hash:
            print(f"Verified {output}")
            return
        raise FileExistsError(f"Preserving modified input: {output}")
    with TemporaryDirectory(dir=folder) as temporary:
        generated = Path(temporary) / output.name
        runpy.run_path(str(folder / "refine.py"))["build"](source, generated, *factors)
        if digest(generated) != deck_hash:
            raise ValueError(f"Generated {case} differs from the validated input")
        generated.rename(output)
    print(f"Verified {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", type=Path, help="FEBio TestSuite checkout")
    parser.add_argument(
        "cases", nargs="*", metavar="CASE", help="bp10, tu01, mi06 (default: all)"
    )
    args = parser.parse_args()
    if not __debug__ or set(args.cases) - CASES.keys():
        parser.error("Use standard Python and select bp10, tu01 or mi06.")
    for case in args.cases or CASES:
        prepare(case, args.suite)
