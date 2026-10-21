#!/bin/bash
# Bundle the exact checkout and public Git dependencies before Slurm dispatch.
set -euo pipefail
out=$(realpath -m "${1:?usage: ci/bundle.sh OUT.tar.zst}")
src=$(git rev-parse --show-toplevel)
cache=${SCI_CACHE:-${HOME}/.cache/sci-bundle}
stage=$(mktemp -d "${TMPDIR:-/tmp}/sci-bundle.XXXXXX")
trap 'rm -rf "$stage"' EXIT
mkdir -p "$cache/mirrors" "$stage/mirrors"
while read -r slug; do
    [[ -z $slug || $slug == \#* ]] && continue
    mkdir -p "$cache/locks"
    exec {lock}> "$cache/locks/${slug//\//_}"
    flock "$lock"
    mirror=$cache/mirrors/$slug.git
    if [ ! -d "$mirror" ]; then
        mkdir -p "$(dirname "$mirror")"
        git clone -q --mirror "https://github.com/$slug.git" "$mirror"
    else
        git -C "$mirror" fetch -q --prune origin '+refs/heads/*:refs/heads/*' '+refs/tags/*:refs/tags/*'
    fi
    mkdir -p "$stage/mirrors/${slug%/*}"
    cp -a "$mirror" "$stage/mirrors/$slug.git"
    exec {lock}>&-
done < "$src/ci/dependencies.txt"
# Use Git's transport to copy reachable objects, without borrowed alternates
# or unrelated local branches from a developer's reference clone.
git clone -q --no-local --single-branch --no-tags "$src" "$stage/src"
if [ "$SCI_TIER" = golden-record ]; then
    git -C "$stage/src" fetch -q --no-tags "$src" "${SCI_BASE_SHA:?base revision missing}"
    git check-ref-format "refs/tags/${SCI_STABLE_REF:?stable release missing}"
    git -C "$stage/src" fetch -q --no-tags "$src" \
        "refs/tags/$SCI_STABLE_REF:refs/tags/$SCI_STABLE_REF"
fi
git -C "$stage/src" remote set-url origin "https://github.com/${GITHUB_REPOSITORY:?GITHUB_REPOSITORY missing}.git"
git -C "$src" rev-parse HEAD > "$stage/revision.txt"
printf '%s\n' "${SCI_TIER:?SCI_TIER missing}" > "$stage/tier.txt"
: > "$stage/task.env"
for key in SCI_BASE_SHA SCI_STABLE_REF SCI_LIBNEO_REF; do
    printf 'export %s=%q\n' "$key" "${!key:-}" >> "$stage/task.env"
done
while IFS= read -r key; do
    printf 'export %s=%q\n' "$key" "${!key}" >> "$stage/task.env"
done < <(compgen -A variable SCI_MATRIX_)
if [ "$SCI_TIER" = golden-record ]; then
    # Only the private golden subset crosses the campus dispatch boundary.
    # Credentials and the private repository's Git database never enter it.
    "$src/ci/prepare-golden.sh" "$stage/golden-data"
    git -C "$src" cat-file -e "${SCI_BASE_SHA:?base revision missing}^{commit}"
    git -C "$src" cat-file -e "${SCI_STABLE_REF:?stable release missing}^{commit}"
fi
tar -C "$stage" -cf - . | zstd -q -T2 -o "$out"
sha256sum "$out"
