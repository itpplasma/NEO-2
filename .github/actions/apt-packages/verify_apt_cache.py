#!/usr/bin/env python3
"""Check restored APT archives against a fresh, authenticated download plan."""

import hashlib
import re
import shlex
import shutil
import sys
from pathlib import Path


def read_plan(path):
    expected = {}
    for line in path.read_text().splitlines():
        if not line.startswith("'"):
            continue
        fields = shlex.split(line)
        if len(fields) != 4:
            raise ValueError("Invalid APT package URI line")
        _, name, size, digest = fields
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._%+:~\-]*\.deb", name):
            raise ValueError("Invalid APT archive filename")
        if not re.fullmatch(r"[0-9]+", size):
            raise ValueError("Invalid APT archive size")
        if not re.fullmatch(r"SHA256:[0-9a-f]{64}", digest):
            raise ValueError("Missing APT archive SHA256")
        record = (int(size), digest.removeprefix("SHA256:"))
        if name in expected and expected[name] != record:
            raise ValueError("Conflicting APT archive metadata")
        expected[name] = record
    return expected


def verify(plan, archives):
    expected = read_plan(plan)
    kept = discarded = 0
    for path in archives.iterdir():
        if path.is_symlink():
            path.unlink()
            discarded += 1
        elif path.name.endswith(".deb"):
            record = expected.get(path.name)
            valid = path.is_file() and record is not None
            if valid:
                size, digest = record
                try:
                    valid = path.stat().st_size == size
                    if valid:
                        with path.open("rb") as stream:
                            actual = hashlib.sha256()
                            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                                actual.update(chunk)
                        valid = actual.hexdigest() == digest
                except OSError:
                    valid = False
            if valid:
                kept += 1
            else:
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
                discarded += 1
    canonical = "".join(
        f"{name} {size} {digest}\n"
        for name, (size, digest) in sorted(expected.items())
    )
    key = hashlib.sha256(canonical.encode()).hexdigest()
    print(f"cache-key={key}")
    print(f"APT archives: kept {kept}, discarded {discarded}", file=sys.stderr)


if __name__ == "__main__":
    try:
        verify(Path(sys.argv[1]), Path(sys.argv[2]))
    except (OSError, ValueError) as error:
        sys.exit(f"APT archive verification failed: {error}")
