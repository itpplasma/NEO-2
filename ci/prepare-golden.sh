#!/bin/bash
# Prepare private runtime inputs on the trusted, campus-connected dispatcher.
set -euo pipefail
out=${1:?usage: ci/prepare-golden.sh OUT}
work=$(mktemp -d "${RUNNER_TEMP:-/tmp}/neo2-data.XXXXXX")
trap 'rm -rf "$work"' EXIT
if [ -n "${GITLAB_ACCESS_TOKEN:-}" ]; then
    cat > "$work/askpass" <<'ASKPASS'
#!/bin/sh
case "$1" in
  *Username*) printf '%s\n' oauth2 ;;
  *Password*) printf '%s\n' "$GITLAB_ACCESS_TOKEN" ;;
esac
ASKPASS
    chmod 700 "$work/askpass"
    export GIT_ASKPASS=$work/askpass GIT_TERMINAL_PROMPT=0
fi
export GIT_LFS_SKIP_SMUDGE=1
git clone --no-checkout https://gitlab.tugraz.at/plasma/data.git "$work/data"
git -C "$work/data" lfs install --local --skip-repo
branch=${SCI_DATA_BRANCH:-main}
if ! git -C "$work/data" show-ref --verify --quiet "refs/remotes/origin/$branch"; then
    if [ "${SCI_GOLDEN_CHANGE:-0}" = 1 ]; then
        echo "golden-change needs plasma/data branch $branch with regenerated records" >&2
        exit 1
    fi
    branch=main
fi
git -C "$work/data" sparse-checkout set TESTS/NEO-2/golden_record
git -C "$work/data" checkout --detach "origin/$branch"
git -C "$work/data" lfs pull --include='TESTS/NEO-2/golden_record/**' --exclude=''
git -C "$work/data" lfs checkout 'TESTS/NEO-2/golden_record/**'
mkdir -p "$out/TESTS/NEO-2"
cp "$work/data/activate.sh" "$out/"
cp -a "$work/data/TESTS/NEO-2/golden_record" "$out/TESTS/NEO-2/"
git -C "$work/data" rev-parse HEAD > "$out/revision.txt"
(
    cd "$out"
    find activate.sh TESTS -type f -print0 | sort -z | xargs -0 sha256sum
) > "$out/inputs.sha256"
# A missing LFS object must fail preparation, never become a solver input.
if grep -rl '^version https://git-lfs.github.com/spec/v1$' "$out/TESTS"; then
    echo 'Unresolved Git LFS pointer in golden inputs' >&2
    exit 1
fi
