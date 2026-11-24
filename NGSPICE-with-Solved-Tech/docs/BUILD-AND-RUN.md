# Build and run on Linux

Use Linux x86-64 and Python 3.11 or later. Building ngspice requires Git, a C compiler, GNU Make, Autoconf, Automake, Libtool, Flex and Bison. On Debian and Ubuntu the build packages are `build-essential git autoconf automake libtool flex bison`. CUDA and the PanguLU SDK are needed only for the solver replay.

Run the examples from the package root. Each invocation writes to a new case output directory, so existing run outputs remain available.

## Native ngspice run and capture

```bash
JOBS=8 ./tools/build-ngspice.sh build
case_id=pdn-interleaved-segments-dc
python3 tools/generate_cases.py "$case_id" --output-dir netlists
python3 tools/run_case.py capture "$case_id" \
    --ngspice build/ngspice-build/src/ngspice \
    --netlist "netlists/$case_id.cir" \
    --output-dir "runs/$case_id/klu"
```

The build script checks out the commit in [upstream.lock.json](../upstream.lock.json) and applies the capture hook in its own build directory. If the build host is offline, provide a local Git repository containing that commit as the second argument. The native output includes `ngspice.log`, the physical operating-point plot `field.raw`, the original equations and reference vectors under `systems/`, and the capture verification record `capture.json`.

With a standard KLU-enabled ngspice binary, `ngspice -n -b -r field.raw netlists/<case>.cir` runs the same deck directly. Use the build from this package when you need the equation captures.

## PanguLU replay

Obtain the PanguLU 26.46 ABI 2 SDK for CUDA 12.9, its runtime dependencies and a valid license for the GPU. Follow the SDK's licensing instructions. Set the SDK and CUDA paths for your installation. If the CUDA toolkit uses a split layout, set `CUDA_RUNTIME_INCLUDE`, `CUDA_RUNTIME_LIB_DIR` and, when needed, `CUDA_COMPILER_INCLUDE`. Add SDK and dependent library directories to `LD_LIBRARY_PATH` as required by your installation.

```bash
export PANGULU_ROOT=/path/to/pangulu-26.46
export CUDA_HOME=/path/to/cuda-12.9
./tools/build-pangulu.sh build
case_id=pdn-interleaved-segments-dc
CUDA_VISIBLE_DEVICES=0 python3 tools/run_case.py replay "$case_id" \
    --pangulu-replay build/pangulu-replay \
    --capture-dir "runs/$case_id/klu" \
    --output-dir "runs/$case_id/pangulu" --host-threads 1
```

The reference PanguLU configuration sets one host OpenMP thread and allocates 40 host CPUs to the GPU job. The program uses the current visible CUDA device. `pangulu.log` shows the CUDA device, synchronized solver phases and checks for every original solve. `result.json` is written only after all three solves pass. It reports native solver time and replay wall time separately.

Generate all five exact decks with `python3 tools/generate_cases.py all --output-dir netlists` in a fresh deck directory, then repeat capture and replay with the case IDs listed in [the case documentation](../cases/README.md). The capture and replay steps can run on separate hosts sharing the capture output. Preserve the complete capture directory when transferring it. See [the native capture validation](../validation/README.md) for the reference checks.
