#!/bin/bash
CLONE_URL="https://github.com/itpplasma/NEO-2.git"

TMP_DIR=${1:-"$(mktemp -d)"}
if [ ! -d "$TMP_DIR" ]; then
    echo "Create temporary dir"
    echo $TMP_DIR
    mkdir $TMP_DIR
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CUR_VER=$(git -C "$SCRIPT_DIR" describe --tags --always --dirty)
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)

if [ "$REF_VER" == "$CUR_VER" ]; then
    echo "Reference version and current version are the same. Exiting."
    exit 1
fi

PROJECT_ROOT_CUR=$(cd "$SCRIPT_DIR/../.." && pwd)

echo "Temporary directory: $TMP_DIR"
echo "Current version: $CUR_VER"
echo "Project root (current): $PROJECT_ROOT_CUR"

main() {
    build "$PROJECT_ROOT_CUR"
    export CURRENT_NEO2_PAR=$PROJECT_ROOT_CUR/build/NEO-2-PAR/neo_2_par.x
    export CURRENT_NEO2_QL=$PROJECT_ROOT_CUR/build/NEO-2-QL/neo_2_ql.x
    get_test_data_and_run
}

get_test_data_and_run() {
    cd "$TMP_DIR"
    local TEST_DIR=TESTS/NEO-2/golden_record/
    # Fetch if not exists
    if [ ! -d "data" ]; then
        if [ -z "$GITLAB_ACCESS_TOKEN" ]; then
            echo "Error: GITLAB_ACCESS_TOKEN is not set"
            exit 1
        fi
        echo "Fetching test data..."
        git clone https://oauth2:${GITLAB_ACCESS_TOKEN}@gitlab.tugraz.at/plasma/data.git
    fi
    cd data
    git checkout $CURRENT_BRANCH
    git pull
    git config lfs.fetchinclude "$TEST_DIR"
    git config lfs.fetchexclude ""
    git lfs pull
    git lfs checkout $TEST_DIR
    cd $TEST_DIR
    make
    #wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc > /dev/null 2>&1
}

build() {
    local PROJECT_ROOT="$1"
    echo "Building NEO-2 in $PROJECT_ROOT"
    cd $PROJECT_ROOT
    make
}

main "$@"
