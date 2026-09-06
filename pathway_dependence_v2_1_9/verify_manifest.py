"""Validate the release manifest after a clean extraction."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "SHA256SUMS.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def included_files() -> set[str]:
    files: set[str] = set()
    for path in ROOT.rglob("*"):
        if not path.is_file() or path == MANIFEST:
            continue
        if "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        files.add(path.relative_to(ROOT).as_posix())
    return files


def main() -> None:
    if not MANIFEST.is_file():
        raise RuntimeError("SHA256SUMS.csv is missing")
    with MANIFEST.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    listed: set[str] = set()
    for row in rows:
        rel = row["path"]
        pure = Path(rel)
        if pure.is_absolute() or ".." in pure.parts or "\\" in rel:
            raise RuntimeError(f"Unsafe manifest path: {rel}")
        if rel in listed:
            raise RuntimeError(f"Duplicate manifest path: {rel}")
        listed.add(rel)
        path = ROOT / pure
        if not path.is_file():
            raise RuntimeError(f"Missing file: {rel}")
        if path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Size mismatch: {rel}")
        if sha256(path) != row["sha256"]:
            raise RuntimeError(f"SHA-256 mismatch: {rel}")
    actual = included_files()
    if listed != actual:
        missing = sorted(listed - actual)
        unlisted = sorted(actual - listed)
        raise RuntimeError(f"Manifest inventory mismatch; missing={missing}, unlisted={unlisted}")
    print(f"Manifest verification passed: {len(rows)} files")


if __name__ == "__main__":
    main()