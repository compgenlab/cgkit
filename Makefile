BIN_DIR := bin
GO_SOURCES := $(shell find . -type f -name '*.go')
VERSION ?= $(shell git describe --tags --always --dirty 2>/dev/null || echo dev)
GIT_HASH ?= $(shell git rev-parse --short HEAD 2>/dev/null)
LDFLAGS := -X 'github.com/compgenlab/cgkit/internal/cmd.Version=$(VERSION)' \
           -X 'github.com/compgenlab/cgkit/internal/cmd.GitHash=$(GIT_HASH)'

# Local builds resolve the github.com/compgenlab/hts dependency via the
# sibling go.work workspace; release builds (no go.work present) use the
# pinned module version from go.mod.

.PHONY: build clean test bump-hts

build: $(BIN_DIR)/cgkit.darwin_arm64 $(BIN_DIR)/cgkit.darwin_amd64 $(BIN_DIR)/cgkit.linux_arm64 $(BIN_DIR)/cgkit.linux_amd64 $(BIN_DIR)/cgkit.windows_amd64.exe

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/cgkit.darwin_arm64: $(GO_SOURCES) | $(BIN_DIR)
	GOOS=darwin GOARCH=arm64 go build -buildvcs=false -ldflags "$(LDFLAGS)" -o $@ .

$(BIN_DIR)/cgkit.darwin_amd64: $(GO_SOURCES) | $(BIN_DIR)
	GOOS=darwin GOARCH=amd64 go build -buildvcs=false -ldflags "$(LDFLAGS)" -o $@ .

$(BIN_DIR)/cgkit.linux_arm64: $(GO_SOURCES) | $(BIN_DIR)
	GOOS=linux GOARCH=arm64 CGO_ENABLED=0 go build -buildvcs=false -ldflags "$(LDFLAGS)" -o $@ .

$(BIN_DIR)/cgkit.linux_amd64: $(GO_SOURCES) | $(BIN_DIR)
	GOOS=linux GOARCH=amd64 CGO_ENABLED=0 go build -buildvcs=false -ldflags "$(LDFLAGS)" -o $@ .

$(BIN_DIR)/cgkit.windows_amd64.exe: $(GO_SOURCES) | $(BIN_DIR)
	GOOS=windows GOARCH=amd64 CGO_ENABLED=0 go build -buildvcs=false -ldflags "$(LDFLAGS)" -o $@ .

clean:
	rm -rf $(BIN_DIR)
	mkdir $(BIN_DIR)

test:
	GOCACHE=/tmp/go-build-cache go test ./...

# Pin the committed go.mod to the latest released hts from GitHub. Run this
# (then commit go.mod/go.sum) when cutting a release, after the matching hts tag
# has been pushed. GOWORK=off so it updates the module pin, not the workspace.
bump-hts:
	GOWORK=off GOPRIVATE=github.com/compgenlab/* GOCACHE=/tmp/go-build-cache \
		go get github.com/compgenlab/hts@latest
	GOWORK=off GOCACHE=/tmp/go-build-cache go mod tidy
