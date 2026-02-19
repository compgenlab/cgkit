BIN_DIR := bin
GO_SOURCES := $(shell rg --files -g "*.go")

.PHONY: build clean test

build: $(BIN_DIR)/cgltk.darwin_arm64 $(BIN_DIR)/cgltk.linux_arm64 $(BIN_DIR)/cgltk.linux_amd64

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/cgltk.darwin_arm64: $(GO_SOURCES) | $(BIN_DIR)
	GOOS=darwin GOARCH=arm64 go build -o $@ .

$(BIN_DIR)/cgltk.linux_arm64: $(GO_SOURCES) | $(BIN_DIR)
	GOOS=linux GOARCH=arm64 go build -o $@ .

$(BIN_DIR)/cgltk.linux_amd64: $(GO_SOURCES) | $(BIN_DIR)
	GOOS=linux GOARCH=amd64 go build -o $@ .

clean:
	rm -rf $(BIN_DIR)
	mkdir $(BIN_DIR)

test:
	GOCACHE=/tmp/go-build-cache go test ./...
