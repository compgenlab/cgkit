#!/bin/sh
# Gate to run before cutting a release. Tags nothing itself.
#
# Two things this protects against.
#
# First, the state being released has to be the state that was tested. The gap is
# subtler than a broken build: work written and verified, then left uncommitted,
# and then either excluded by `git commit --amend` (which takes only what is
# staged) or destroyed by `git checkout <file>`. The suite still passes, because
# it ran against the working tree. The release ships without the fix. That has
# happened here, in both repositories.
#
# Second, and specific to cgkit: local builds resolve cghts through a sibling
# go.work workspace, while a release resolves the pin in go.mod. Those can differ,
# so "the tests pass" locally says nothing about whether the released artifact
# builds. The Makefile deliberately does not set GOWORK=off for ordinary targets;
# this gate does, because that is what CI and `go install` see.
#
# Every check reports rather than exiting early, so one run tells you everything.
#
# Usage: release-check.sh [release-branch] [tag]
#
# Note the second argument is a TAG, not the Makefile's VERSION -- that is already
# defined from `git describe` for ldflags, so reusing it made this gate refuse every
# run against the current describe string.

set -u

branch_expected="${1:-main}"
version="${2:-}"
fail=0

say() { printf '%s\n' "$*"; }

# 1. A dirty tree means the tested state is not the committed state.
if [ -n "$(git status --porcelain)" ]; then
	say "REFUSING: working tree is not clean."
	say "  Uncommitted changes mean what you tested is not what you would release."
	git status --short | sed 's/^/    /'
	fail=1
fi

# 2. Releases come from one branch.
branch="$(git rev-parse --abbrev-ref HEAD)"
if [ "$branch" != "$branch_expected" ]; then
	say "REFUSING: on branch '$branch', not '$branch_expected'."
	say "  Pass RELEASE_BRANCH to override."
	fail=1
fi

# 3. Local and remote must agree.
git fetch -q origin "$branch" 2>/dev/null || true
if git rev-parse --verify --quiet "origin/$branch" >/dev/null 2>&1; then
	if [ "$(git rev-parse HEAD)" != "$(git rev-parse "origin/$branch")" ]; then
		say "REFUSING: HEAD differs from origin/$branch."
		git --no-pager log --oneline "origin/$branch..HEAD" 2>/dev/null | sed 's/^/    unpushed: /'
		git --no-pager log --oneline "HEAD..origin/$branch" 2>/dev/null | sed 's/^/    missing:  /'
		fail=1
	fi
else
	say "NOTE: origin/$branch not found; skipping the sync check."
fi

# 4. A go.work makes local builds resolve a sibling checkout instead of the pin,
#    so its presence means the verification below would not be measuring the
#    thing being released. The gate runs with GOWORK=off regardless, but say so:
#    an unreleased local cghts change is the usual reason the pin is stale.
if [ -f go.work ]; then
	say "NOTE: go.work is present. This gate uses GOWORK=off, so it checks the"
	say "      go.mod pin rather than your local cghts checkout -- which is the"
	say "      point, but means local changes to cghts are invisible here."
fi

# 5. Never re-cut an existing tag.
if [ -n "$version" ]; then
	if git rev-parse --verify --quiet "refs/tags/$version" >/dev/null 2>&1; then
		say "REFUSING: tag $version already exists locally."
		fail=1
	fi
	if git ls-remote --exit-code --tags origin "refs/tags/$version" >/dev/null 2>&1; then
		say "REFUSING: tag $version already exists on origin."
		fail=1
	fi
fi

if [ "$fail" -ne 0 ]; then
	say ""
	say "release-check FAILED -- do not release."
	exit 1
fi

say "state: clean, on $branch_expected, in sync with origin."
say "cghts pin: $(grep 'compgenlab/cghts' go.mod | tr -s ' \t' ' ' | sed 's/^ //')"
