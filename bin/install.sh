#!/bin/sh
set -eu
bundle_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
command_dir="${HOME}/.local/bin"
mkdir -p "$command_dir"
target="$command_dir/myloc"
if [ -e "$target" ] || [ -L "$target" ]; then
    if [ ! -L "$target" ]; then
        printf 'Refusing to replace existing file: %s\n' "$target" >&2
        exit 1
    fi
fi
ln -sfn "$bundle_dir/bin/myloc" "$target"
printf 'Installed %s\n' "$target"
printf 'Ensure ~/.local/bin is in PATH, then run: myloc --version\n'
printf 'Keep the complete bundle; rerun install.sh if you move it.\n'
