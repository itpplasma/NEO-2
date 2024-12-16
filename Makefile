BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all ninja test install clean
all: ninja

$(BUILD_NINJA):
	cmake --preset default

ninja: $(BUILD_NINJA)
	cmake --build --preset default

test: ninja
	cd $(BUILD_DIR) && ctest

clean:
	rm -rf $(BUILD_DIR)
