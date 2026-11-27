#!/usr/bin/env bash

package_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
work_root="${OPM_SOLVED_TECH_WORK_ROOT:-$package_dir/work}"

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

info() {
    printf '[opm-solved-tech] %s\n' "$*"
}

validate_route() {
    case "${1:-}" in
        linglongamg|xianzhikrylov) ;;
        *)
            die "Invalid --route='${1:-}'. Expected one of: linglongamg, xianzhikrylov."
            ;;
    esac
}

find_linglong_header() {
    local root=$1
    for candidate in \
        "${LINGLONGAMG_INCLUDE_DIR:-}/linglongamg.h" \
        "$root/include/linglongamg.h" \
        "$root/lib/linglongamg.h" \
        "$root/linglongamg.h"; do
        if [[ -n $candidate && -r $candidate ]]; then
            printf '%s\n' "$(dirname -- "$candidate")"
            return 0
        fi
    done
    return 1
}

find_linglong_library() {
    local root=$1
    for candidate in \
        "${LINGLONGAMG_LIBRARY:-}" \
        "$root/lib/liblinglongamg.so" \
        "$root/lib64/liblinglongamg.so" \
        "$root/liblinglongamg.so"; do
        if [[ -n $candidate && -r $candidate ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

find_xianzhi_include() {
    local root=$1
    for candidate in \
        "${XIANZHIKRYLOV_INCLUDE_DIR:-}" \
        "$root/include" \
        "$root"; do
        if [[ -n $candidate && -r "$candidate/XianzhiKrylovAPI.h" && -r "$candidate/XianzhiKrylovData.h" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

find_xianzhi_library() {
    local root=$1
    local candidate
    if [[ -n ${XIANZHIKRYLOV_LIBRARY:-} && -r ${XIANZHIKRYLOV_LIBRARY:-} ]]; then
        printf '%s\n' "$XIANZHIKRYLOV_LIBRARY"
        return 0
    fi
    for dir in "$root/lib" "$root/lib64" "$root"; do
        if [[ -r "$dir/libXianzhiKrylov.so" ]]; then
            printf '%s\n' "$dir/libXianzhiKrylov.so"
            return 0
        fi
        candidate=$(compgen -G "$dir/libXianzhiKrylov.so.*" | head -n 1 || true)
        if [[ -n $candidate && -r $candidate ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

find_license_file() {
    local override=$1
    local root=$2
    if [[ -n $override ]]; then
        [[ -r $override ]] && { printf '%s\n' "$override"; return 0; }
        return 1
    fi
    for candidate in \
        "$root/license.key" \
        "$root/license/license.key"; do
        if [[ -r $candidate ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    local match
    for glob_pat in \
        "$root/licenses/*.key" \
        "$root/*.lic" \
        "$root/license/*.lic"; do
        match=$(compgen -G "$glob_pat" | head -n 1 || true)
        if [[ -n $match && -r $match ]]; then
            printf '%s\n' "$match"
            return 0
        fi
    done
    return 1
}
