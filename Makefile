BIN_DIR := bin
GO_SOURCES := $(shell find . -type f -name '*.go')
VERSION ?= $(shell git describe --tags --always --dirty 2>/dev/null || echo dev)
GIT_HASH ?= $(shell git rev-parse --short HEAD 2>/dev/null)
# -s -w drop the symbol table and DWARF debug info, which is 30% of the binary
# (37.9 MB -> 26.4 MB for linux_amd64, measured). Panic tracebacks are NOT
# affected: Go's runtime carries its own pclntab for those, so a user-reported
# crash still names functions and line numbers. What is lost is attaching a
# debugger to a shipped binary and resolving symbols in a profile -- both of
# which mean rebuilding from source anyway, since this is a distributed CLI
# rather than a service anyone attaches to in place. Build without them when
# you need that: make build LDFLAGS_STRIP=
LDFLAGS_STRIP ?= -s -w
LDFLAGS := $(LDFLAGS_STRIP) \
           -X 'github.com/compgenlab/cgkit/internal/cmd.Version=$(VERSION)' \
           -X 'github.com/compgenlab/cgkit/internal/cmd.GitHash=$(GIT_HASH)'

# Local builds resolve the github.com/compgenlab/cghts dependency via the
# sibling go.work workspace; release builds (no go.work present) use the
# pinned module version from go.mod.

RELEASE_BRANCH ?= main

.PHONY: build clean test bump-cghts release-check

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

# Pin the committed go.mod to the latest released cghts from GitHub. Run this
# (then commit go.mod/go.sum) when cutting a release, after the matching cghts
# tag has been pushed. GOWORK=off so it updates the module pin, not the workspace.
bump-cghts:
	GOWORK=off GOPRIVATE=github.com/compgenlab/* GOCACHE=/tmp/go-build-cache \
		go get github.com/compgenlab/cghts@latest
	GOWORK=off GOCACHE=/tmp/go-build-cache go mod tidy

# Gate to run BEFORE tagging a release; it tags nothing itself. Verifies the
# state being released is the state that was tested -- clean tree, on the release
# branch, in sync with origin, tag unused -- then vets, tests and cross-builds
# with GOWORK=off.
#
# GOWORK=off matters here specifically: ordinary targets resolve cghts through
# the sibling go.work workspace, but a release resolves the go.mod pin. Passing
# tests locally says nothing about whether the released artifact builds.
#
#   make release-check
#   make release-check TAG=v1.2.3
#   make release-check RELEASE_BRANCH=topic   # exercise it off main
release-check:
	@sh scripts/release-check.sh "$(RELEASE_BRANCH)" "$(TAG)"
	@echo "--- verifying against the go.mod pin (GOWORK=off), not the workspace ---"
	GOWORK=off GOCACHE=/tmp/go-build-cache go vet ./...
	@out=$$(gofmt -l . 2>/dev/null); if [ -n "$$out" ]; then echo "gofmt needed:"; echo "$$out"; exit 1; fi
	GOWORK=off GOCACHE=/tmp/go-build-cache go test ./...
	GOWORK=off $(MAKE) --no-print-directory build
	@echo
	@echo "release-check PASSED -- safe to tag."
