#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/environment.sh"
: "${CS_SOURCE_DIR:?Set CS_SOURCE_DIR to the source tree used to build CS_PREFIX}"
: "${LINGLONGAMG_LICENSE_FILE:?Set LINGLONGAMG_LICENSE_FILE to the device authorization file}"

sdk_lib=$LINGLONGAMG_ROOT/lib
runtime=$CS_LINGLONGAMG_BUILD/linglongamg
compiler=${CXX:-g++}
adapter_dir=$CS_LINGLONGAMG_DIR/overlay/adapter
host_source="${CS_LINGLONGAMG_HOST_SLES_SOURCE:-$CS_SOURCE_DIR/src/alge/cs_sles_default.cpp}"

# The host and the adapter must agree on the registration type that selects ParCSR.
host_type_file=$(mktemp)
trap 'rm -f -- "$host_type_file"' EXIT
strings "$CS_PREFIX/lib/libsaturne.so" > "$host_type_file"
host_type=$(python3 "$CS_LINGLONGAMG_DIR/scripts/host_pressure_type.py" "$host_source" "$host_type_file")

mkdir -p "$runtime/lib"
install -m 755 "$sdk_lib/liblinglongamg.so" "$runtime/lib/"
install -m 644 "$sdk_lib/delivery.receipt" "$runtime/lib/"
install -m 600 "$LINGLONGAMG_LICENSE_FILE" "$runtime/license.key"

common=(-std=c++17 -O2 -g -Wall -Wextra -DHAVE_CONFIG_H
  -I"$CS_PREFIX/include/code_saturne" -I"$CS_PREFIX/include" -I"$CS_SOURCE_DIR/src"
  -I"$sdk_lib" -I"$CUDA_HOME/include" -I"$HYPRE_ROOT/include"
  -L"$CS_PREFIX/lib" -L"$runtime/lib" -L"$CUDA_HOME/lib64" -L"$HYPRE_ROOT/lib"
  -Wl,-rpath-link,"$CS_LINGLONGAMG_LIBPATH" -Wl,-rpath,"$CS_LINGLONGAMG_LIBPATH")
if [[ -n ${MPI_HOME:-} ]]; then
  common+=(-I"$MPI_HOME/include")
fi

"$compiler" "${common[@]}" -D"CS_SLES_PRESSURE_BACKEND_TYPE=\"$host_type\"" \
  -fPIC -shared "$adapter_dir/cs_linglongamg_adapter.cpp" \
  -lsaturne -lHYPRE -llinglongamg -lcudart -Wl,--no-undefined \
  -o "$CS_LINGLONGAMG_ADAPTER_LIBRARY"
"$compiler" "${common[@]}" -fPIC -shared "$adapter_dir/cs_pressure_timing.cpp" \
  -lsaturne -ldl -Wl,--no-undefined -o "$CS_LINGLONGAMG_BUILD/libcs_pressure_timing.so"

hostname > "$CS_LINGLONGAMG_BUILD/build-host.txt"
"$compiler" --version > "$CS_LINGLONGAMG_BUILD/compiler.txt"
sha256sum "$adapter_dir/cs_linglongamg_adapter.cpp" "$adapter_dir/cs_pressure_timing.cpp" \
  "$runtime/lib/liblinglongamg.so" "$CS_LINGLONGAMG_ADAPTER_LIBRARY" \
  "$CS_LINGLONGAMG_BUILD/libcs_pressure_timing.so" > "$CS_LINGLONGAMG_BUILD/build.sha256"
printf 'Built LinglongAMG adapter and native pressure timer in %s\n' "$CS_LINGLONGAMG_BUILD"
