#!/usr/bin/env bash
# Compatibility shim — use alkyl.sh setup-key perplexity <KEY> instead
exec bash "$(dirname "${BASH_SOURCE[0]}")/alkyl.sh" setup-key perplexity "$@"
