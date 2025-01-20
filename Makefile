BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja
CONFIG ?= Release

.PHONY: all ninja test install clean
all: ninja

$(BUILD_NINJA):
	cmake --preset default -DCMAKE_COLOR_DIAGNOSTICS=ON -DCMAKE_BUILD_TYPE=$(CONFIG)

ninja: $(BUILD_NINJA)
	cmake --build --preset default

test: ninja
	cd $(BUILD_DIR) && ctest

doc: $(BUILD_NINJA)
	cmake --build --preset default --target doc

clean:
	rm -rf $(BUILD_DIR)
