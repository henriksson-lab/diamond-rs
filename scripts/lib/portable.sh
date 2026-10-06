#!/usr/bin/env bash

# Small portability helpers shared by target-native validation scripts.
# GitHub's Linux and Windows bash environments provide sha256sum, while the
# stock macOS image provides shasum. Keep the emitted value identical.
sha256_file() {
    if (($# != 1)); then
        echo "error: sha256_file expects exactly one path" >&2
        return 2
    fi
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum -- "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 -- "$1" | awk '{print $1}'
    else
        echo "error: neither sha256sum nor shasum is available" >&2
        return 2
    fi
}
