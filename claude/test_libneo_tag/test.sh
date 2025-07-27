#!/bin/bash

# Test script for LIBNEO_GIT_TAG functionality

set -e

echo "Testing LIBNEO_GIT_TAG build functionality"
echo "=========================================="

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Save and unset CODE environment variable
SAVED_CODE=""
if [ -n "$CODE" ]; then
    echo -e "${BLUE}INFO: CODE environment variable is set to: $CODE${NC}"
    echo "Temporarily unsetting CODE for this test..."
    SAVED_CODE="$CODE"
    unset CODE
fi

# Ensure CODE is restored on exit
restore_code() {
    if [ -n "$SAVED_CODE" ]; then
        export CODE="$SAVED_CODE"
        echo -e "${BLUE}Restored CODE environment variable${NC}"
    fi
}
trap restore_code EXIT

# Get the latest commit hash from libneo
echo -e "${BLUE}Fetching latest libneo commit hash...${NC}"
LIBNEO_HEAD=$(git ls-remote https://github.com/itpplasma/libneo.git HEAD | cut -f1)
echo "Latest libneo commit: $LIBNEO_HEAD"

# Use v0.0.3 tag for comparison
echo -e "${BLUE}Using libneo v1.0.0 tag for comparison...${NC}"
LIBNEO_OLD="v1.0.0"
echo "Using tag for reference: $LIBNEO_OLD"

# Get current branch
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
echo -e "${BLUE}Current NEO-2 branch: $CURRENT_BRANCH${NC}"

# Create test build directories
TEST_DIR="test_builds"
mkdir -p $TEST_DIR

echo -e "\n${BLUE}Test 1: Default build (no LIBNEO_GIT_TAG)${NC}"
echo "Cloning NEO-2 repository..."
cd $TEST_DIR
git clone https://github.com/itpplasma/NEO-2.git default
cd default
git checkout $CURRENT_BRANCH
echo "Running make (standard build process)..."
# Run make and capture output
make clean 2>/dev/null || true
if make 2>&1 | tee build_output.txt; then
    echo -e "${GREEN}✓ Default build completed successfully${NC}"
    
    # Check if it uses branch matching
    if grep -q "Using libneo branch" build_output.txt; then
        echo -e "${GREEN}✓ Default build uses branch matching as expected${NC}"
        grep "Using libneo" build_output.txt | head -1
    else
        echo -e "${RED}✗ Default build did not show branch matching${NC}"
    fi
    
    # Check if executable was created
    if [ -f "build/NEO-2-QL/neo_2_ql.x" ]; then
        echo -e "${GREEN}✓ neo_2_ql.x executable created${NC}"
    else
        echo -e "${RED}✗ neo_2_ql.x executable not found${NC}"
    fi
else
    echo -e "${RED}✗ Default build failed${NC}"
    echo "Cannot proceed with tag version test"
    exit 1
fi

cd ../..

echo -e "\n${BLUE}Test 2: Build with LIBNEO_GIT_TAG${NC}"
echo "Cloning NEO-2 repository..."
cd $TEST_DIR
git clone https://github.com/itpplasma/NEO-2.git with_tag
cd with_tag
git checkout $CURRENT_BRANCH

echo "Running make with LIBNEO_GIT_TAG=$LIBNEO_OLD..."
# Run make with LIBNEO_GIT_TAG
make clean 2>/dev/null || true
if make LIBNEO_GIT_TAG=$LIBNEO_OLD 2>&1 | tee build_output.txt; then
    echo -e "${GREEN}✓ Tagged build completed successfully${NC}"
    
    # Check if it uses the specified tag
    if grep -q "Using libneo tag $LIBNEO_OLD" build_output.txt; then
        echo -e "${GREEN}✓ Build with tag uses specified commit as expected${NC}"
        grep "Using libneo" build_output.txt | head -1
    else
        echo -e "${RED}✗ Build did not show tag usage${NC}"
    fi
    
    # Check if executable was created
    if [ -f "build/NEO-2-QL/neo_2_ql.x" ]; then
        echo -e "${GREEN}✓ neo_2_ql.x executable created${NC}"
    else
        echo -e "${RED}✗ neo_2_ql.x executable not found${NC}"
    fi
else
    echo -e "${RED}✗ Tagged build failed${NC}"
    echo "Cannot complete test"
    exit 1
fi

cd ../..

echo -e "\n${BLUE}Test 3: Verify different libneo versions${NC}"
# Check if the two builds actually fetched different versions
if [ -d "$TEST_DIR/default/build/libneo" ] || [ -d "$TEST_DIR/default/build/_deps/libneo-src/.git" ]; then
    if [ -d "$TEST_DIR/default/build/libneo/.git" ]; then
        DEFAULT_COMMIT=$(cd $TEST_DIR/default/build/libneo && git rev-parse HEAD 2>/dev/null || echo "N/A")
    else
        DEFAULT_COMMIT=$(cd $TEST_DIR/default/build/_deps/libneo-src && git rev-parse HEAD 2>/dev/null || echo "N/A")
    fi
    
    if [ -d "$TEST_DIR/with_tag/build/_deps/libneo-src/.git" ]; then
        TAGGED_COMMIT=$(cd $TEST_DIR/with_tag/build/_deps/libneo-src && git rev-parse HEAD)
    else
        TAGGED_COMMIT="N/A"
    fi
    
    echo "Default build libneo commit: $DEFAULT_COMMIT"
    echo "Tagged build libneo commit: $TAGGED_COMMIT"
    
    if [ "$DEFAULT_COMMIT" != "N/A" ] && [ "$TAGGED_COMMIT" != "N/A" ] && [ "$DEFAULT_COMMIT" != "$TAGGED_COMMIT" ]; then
        echo -e "${GREEN}✓ Different libneo versions were fetched${NC}"
    else
        echo -e "${RED}✗ Could not verify different libneo versions${NC}"
    fi
else
    echo -e "${RED}✗ Could not find libneo in build directories${NC}"
fi


echo -e "\n${BLUE}Cleaning up test builds...${NC}"
rm -rf $TEST_DIR

echo -e "\n${GREEN}Test completed!${NC}"
