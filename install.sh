#!/usr/bin/env bash
# Compatibility shim — use alkyl.sh install instead
exec bash "$(dirname "${BASH_SOURCE[0]}")/alkyl.sh" install "$@"
