#!/bin/bash
# Resolve the pull request behind a golden run. Its head branch names the
# plasma/data branch with its records, and its golden-change label selects
# report-only base and stable comparisons. A merge-queue ref carries only the
# PR number, so the branch and labels come from the API.
set -euo pipefail
ref=${GITHUB_REF_NAME:?}
branch=${GITHUB_HEAD_REF:-$ref}
pr=${PR_NUMBER:-}
if [ "${GITHUB_EVENT_NAME:?}" = merge_group ]; then
    if ! [[ $ref =~ /pr-([0-9]+)-[0-9a-f]{40}$ ]]; then
        echo "Cannot read the PR number from queue ref $ref" >&2
        exit 2
    fi
    pr=${BASH_REMATCH[1]}
fi
change=0
if [ -n "$pr" ]; then
    json=$(curl -fsS -H "Authorization: Bearer ${GH_TOKEN:?}" \
        "${GITHUB_API_URL:-https://api.github.com}/repos/${GITHUB_REPOSITORY:?}/pulls/$pr")
    branch=$(jq -er .head.ref <<< "$json")
    if jq -e 'any(.labels[]; .name == "golden-change")' <<< "$json" > /dev/null; then
        change=1
    fi
fi
printf 'SCI_DATA_BRANCH=%s\nSCI_GOLDEN_CHANGE=%s\n' "$branch" "$change" >> "${GITHUB_ENV:?}"
