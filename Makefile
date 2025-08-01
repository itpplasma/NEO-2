BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja
CONFIG ?= Release

.PHONY: all ninja test install clean coverage clean-coverage
all: ninja

$(BUILD_NINJA):
	cmake --preset default -DCMAKE_COLOR_DIAGNOSTICS=ON -DCMAKE_BUILD_TYPE=$(CONFIG) $(if $(LIBNEO_GIT_TAG),-DLIBNEO_GIT_TAG=$(LIBNEO_GIT_TAG))

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
	cmake --preset default -DCMAKE_BUILD_TYPE=Debug
	cmake --build --preset default
	@echo "Running tests with coverage..."
	cd $(BUILD_DIR) && ctest --test-dir TEST --output-on-failure
	@echo "Capturing coverage data..."
	cd $(BUILD_DIR) && lcov --capture --directory . --output-file coverage.info \
		--rc branch_coverage=1 \
		--ignore-errors inconsistent,mismatch
	@echo "Filtering coverage data..."
	cd $(BUILD_DIR) && lcov --remove coverage.info \
		'*/build/*' \
		'*/TEST/*' \
		'*/libneo/*' \
		'*/thirdparty/*' \
		'*/DOC/*' \
		'*/MULTI-SPEC-TOOLS/*' \
		'*/tools/*' \
		--output-file coverage_filtered.info \
		--ignore-errors unused
	@echo "Generating HTML report..."
	cd $(BUILD_DIR) && genhtml coverage_filtered.info --output-directory coverage_html \
		--branch-coverage \
		--legend
	@echo "=== Coverage Summary ==="
	@cd $(BUILD_DIR) && lcov --summary coverage_filtered.info
	@if command -v lcov_cobertura >/dev/null 2>&1; then \
		echo "Generating XML report for CI/CD..."; \
		cd $(BUILD_DIR) && lcov_cobertura coverage_filtered.info -o coverage.xml; \
	else \
		echo "Note: Install lcov_cobertura (pip install lcov-cobertura) to generate XML reports"; \
	fi
	@echo "Coverage report generated in $(BUILD_DIR)/coverage_html/index.html"

clean-coverage:
	rm -rf $(BUILD_DIR)/coverage_html/
	rm -f $(BUILD_DIR)/coverage.info $(BUILD_DIR)/coverage_filtered.info $(BUILD_DIR)/coverage.xml
	find $(BUILD_DIR) -name "*.gcov" -o -name "*.gcda" -o -name "*.gcno" -delete 2>/dev/null || true
