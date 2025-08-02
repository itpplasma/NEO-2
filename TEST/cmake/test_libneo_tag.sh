#!/bin/bash

set -e

GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

[[ -n "$CODE" ]] && { SAVED_CODE="$CODE"; unset CODE; }
trap '[[ -n "$SAVED_CODE" ]] && export CODE="$SAVED_CODE"' EXIT

LIBNEO_TAG="v1.0.0"
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
TEST_DIR="test_builds"

log() { echo -e "${BLUE}$1${NC}"; }
pass() { echo -e "${GREEN}✓ $1${NC}"; }
fail() { echo -e "${RED}✗ $1${NC}"; }

run_test() {
    local name="$1"
    local make_args="$2"
    local expected_msg="$3"
    
    log "\nTest: $name"
    (
        mkdir -p "$TEST_DIR/$name"
        cd "$TEST_DIR/$name"
        git clone -q --depth 1 --branch "$CURRENT_BRANCH" https://github.com/itpplasma/NEO-2.git .
        
        if make clean >/dev/null 2>&1 && make $make_args 2>&1 | tee build.log; then
            pass "Build completed"
	    (grep -q "$expected_msg" build.log) && pass "Found: $expected_msg" || fail "Missing: $expected_msg"
            
            for exe in "NEO-2-QL/neo_2_ql.x" "NEO-2-PAR/neo_2_par.x"; do
                [[ -f "build/$exe" ]] && pass "$exe created" || fail "$exe not found"
            done
            
            for dir in build/libneo build/_deps/libneo-src; do
                [[ -d "$dir/.git" ]] && { cd "$dir" && git rev-parse HEAD && break; } || true
            done 2>/dev/null
        else
            fail "Build failed"
            return 1
        fi
    )
}

log "Testing LIBNEO_GIT_TAG functionality\n"

DEFAULT_COMMIT=$(run_test "default" "" "Using libneo branch")
TAGGED_COMMIT=$(run_test "with_tag" "LIBNEO_GIT_TAG=$LIBNEO_TAG" "Using libneo tag $LIBNEO_TAG")

if [[ -n "$DEFAULT_COMMIT" && -n "$TAGGED_COMMIT" && "$DEFAULT_COMMIT" != "$TAGGED_COMMIT" ]]; then
    pass "Different libneo versions fetched"
else
    fail "Could not verify different libneo versions"
fi

rm -rf "$TEST_DIR"
log "\nTest completed!"
