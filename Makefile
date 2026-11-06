BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja
CONFIG ?= Release

# CMake's automatic compiler search can select an unusable `flang` executable
# before gfortran (for example, when an incomplete LLVM installation is on
# PATH).  NEO-2 is developed and tested with gfortran, so use it by default.
# Keep both CMake's and the shell's compiler overrides working for users who
# need a different Fortran compiler.
ifneq ($(filter command line environment,$(origin CMAKE_Fortran_COMPILER)),)
FORTRAN_COMPILER_ARG := -DCMAKE_Fortran_COMPILER=$(CMAKE_Fortran_COMPILER)
else ifneq ($(filter command line environment,$(origin FC)),)
FORTRAN_COMPILER_ARG := -DCMAKE_Fortran_COMPILER=$(FC)
else
DEFAULT_FORTRAN_COMPILER := $(shell command -v gfortran 2>/dev/null)
ifneq ($(strip $(DEFAULT_FORTRAN_COMPILER)),)
FORTRAN_COMPILER_ARG := -DCMAKE_Fortran_COMPILER=$(DEFAULT_FORTRAN_COMPILER)
endif
endif

ifneq ($(filter command line environment,$(origin LIBNEO_GIT_TAG)),)
$(error LIBNEO_GIT_TAG is deprecated; use LIBNEO_REF instead)
endif
ifneq ($(filter command line environment,$(origin LIBNEO_BRANCH)),)
$(error LIBNEO_BRANCH is deprecated; use LIBNEO_REF instead)
endif

# Ignore an ambient LIBNEO_REF so the shell can't change the libneo fetch.
ifeq ($(origin LIBNEO_REF),environment)
LIBNEO_REF :=
endif

.PHONY: all ninja test install clean coverage clean-coverage
all: ninja

$(BUILD_NINJA):
	cmake --preset default -DCMAKE_COLOR_DIAGNOSTICS=ON -DCMAKE_BUILD_TYPE=$(CONFIG) $(FORTRAN_COMPILER_ARG) $(if $(LIBNEO_REF),-DLIBNEO_REF=$(LIBNEO_REF))

ninja: $(BUILD_NINJA)
	cmake --build --preset default

test: ninja
	cd $(BUILD_DIR) && ctest --test-dir TEST --output-on-failure

doc: $(BUILD_NINJA)
	cmake --build --preset default --target doc

clean:
	rm -rf $(BUILD_DIR)

coverage: clean
	@echo "=== Generating code coverage with lcov ==="
	cmake --preset default -DCMAKE_BUILD_TYPE=Coverage
	cmake --build --preset default
	@echo "Running tests with coverage..."
	cd $(BUILD_DIR) && ctest --test-dir TEST --output-on-failure
	@echo "Capturing coverage data..."
	cd $(BUILD_DIR) && lcov --capture --directory . --output-file coverage.info \
		--rc branch_coverage=1 \
		--rc geninfo_unexecuted_blocks=1 \
		--ignore-errors inconsistent,mismatch,empty,unused
	@echo "Filtering coverage data..."
	cd $(BUILD_DIR) && lcov --remove coverage.info \
		'*/build/*' \
		'*/TEST/*' \
		'*/libneo/*' \
		'*/thirdparty/*' \
		'*/DOC/*' \
		'*/MULTI-SPEC-TOOLS/*' \
		'*/tools/*' \
		'/usr/*' \
		'/tmp/*' \
		--output-file coverage_filtered.info \
		--rc branch_coverage=1 \
		--ignore-errors unused,empty
	@echo "Generating HTML report..."
	cd $(BUILD_DIR) && genhtml coverage_filtered.info --output-directory coverage_html \
		--branch-coverage \
		--legend \
		--ignore-errors source || echo "HTML generation completed with warnings"
	@echo "=== Coverage Summary ==="
	@cd $(BUILD_DIR) && lcov --summary coverage_filtered.info || echo "No coverage data found"
	@echo "Coverage report generated in $(BUILD_DIR)/coverage_html/index.html"

clean-coverage:
	rm -rf $(BUILD_DIR)/coverage_html/
	rm -f $(BUILD_DIR)/coverage.info $(BUILD_DIR)/coverage_filtered.info $(BUILD_DIR)/coverage.xml
	find $(BUILD_DIR) -name "*.gcov" -o -name "*.gcda" -o -name "*.gcno" -delete 2>/dev/null || true
