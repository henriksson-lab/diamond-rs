#!/usr/bin/env bash
set -euo pipefail

zig_bin=${ZIG:-zig}
exec "$zig_bin" c++ -target x86_64-linux-musl "$@"
