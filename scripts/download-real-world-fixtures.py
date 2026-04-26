#!/usr/bin/env python3
"""Download small public read datasets used by ignored real-world integration tests.

Example:
    scripts/download-real-world-fixtures.py \
      --dest /husky/henriksson/for_claude/skesa/external
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import sys
import urllib.parse
import urllib.request
from pathlib import Path


DEFAULT_MANIFEST = Path(__file__).resolve().parent.parent / "tests/external/real_world_fixtures.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download public real-world fixtures referenced by ignored integration tests."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="Manifest JSON describing downloadable fixtures.",
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=Path("tests/external-data"),
        help="Directory where fixture files are stored.",
    )
    parser.add_argument(
        "--fixture",
        action="append",
        default=[],
        help="Fixture id to download. Repeat to fetch multiple fixtures. Defaults to all.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Redownload files even if the expected checksum already exists.",
    )
    return parser.parse_args()


def load_manifest(path: Path) -> dict:
    return json.loads(path.read_text())


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensure_download(url: str, dest: Path, expected_md5: str | None, force: bool) -> None:
    if dest.exists() and not force:
        if expected_md5 is None or md5_file(dest) == expected_md5:
            print(f"reuse  {dest}")
            return
        print(f"stale  {dest} (checksum mismatch), redownloading", file=sys.stderr)

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    if tmp.exists():
        tmp.unlink()

    print(f"fetch  {url}")
    with urllib.request.urlopen(url) as response, tmp.open("wb") as handle:
        shutil.copyfileobj(response, handle)

    if expected_md5 is not None:
        observed = md5_file(tmp)
        if observed != expected_md5:
            tmp.unlink(missing_ok=True)
            raise SystemExit(
                f"checksum mismatch for {dest.name}: expected {expected_md5}, observed {observed}"
            )

    tmp.replace(dest)
    print(f"saved  {dest}")


def write_fastq_subset(src: Path, dest: Path, records: int, force: bool) -> None:
    if dest.exists() and not force:
        print(f"reuse  {dest}")
        return

    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"subset {src} -> {dest} ({records} records)")
    opener = gzip.open if src.suffix == ".gz" else open
    with opener(src, "rt", encoding="utf-8") as reader, dest.open("wt", encoding="utf-8") as writer:
        for _ in range(records):
            lines = [reader.readline() for _ in range(4)]
            if not lines[0]:
                break
            if any(line == "" for line in lines[1:]):
                raise SystemExit(f"truncated FASTQ record while subsetting {src}")
            writer.writelines(lines)


def ena_fastq_files(accession: str) -> list[dict[str, str]]:
    params = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": "run_accession,fastq_ftp,fastq_md5,fastq_bytes",
            "format": "tsv",
            "download": "false",
        }
    )
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?{params}"
    with urllib.request.urlopen(url) as response:
        rows = response.read().decode("utf-8").strip().splitlines()

    if len(rows) < 2:
        raise SystemExit(f"ENA file report for {accession} did not return any FASTQ rows")

    header = rows[0].split("\t")
    values = rows[1].split("\t")
    record = dict(zip(header, values, strict=True))
    ftp_paths = [part for part in record.get("fastq_ftp", "").split(";") if part]
    md5s = [part for part in record.get("fastq_md5", "").split(";") if part]
    if not ftp_paths:
        raise SystemExit(f"ENA file report for {accession} did not include fastq_ftp paths")

    files = []
    for index, ftp_path in enumerate(ftp_paths):
        url = ftp_path if "://" in ftp_path else f"https://{ftp_path}"
        files.append(
            {
                "name": Path(urllib.parse.urlparse(url).path).name,
                "url": url,
                "md5": md5s[index] if index < len(md5s) else None,
            }
        )
    return files


def main() -> int:
    args = parse_args()
    manifest = load_manifest(args.manifest)
    fixtures = manifest["fixtures"]
    selected = set(args.fixture)
    if selected:
        fixtures = [fixture for fixture in fixtures if fixture["id"] in selected]
        missing = selected - {fixture["id"] for fixture in fixtures}
        if missing:
            raise SystemExit(f"unknown fixture ids: {', '.join(sorted(missing))}")

    args.dest.mkdir(parents=True, exist_ok=True)
    for fixture in fixtures:
        fixture_dir = args.dest / fixture["id"]
        file_infos = fixture.get("files") or ena_fastq_files(fixture["accession"])
        for file_info in file_infos:
            ensure_download(
                file_info["url"],
                fixture_dir / file_info["name"],
                file_info.get("md5"),
                args.force,
            )
        subset_records = fixture.get("subset_records_per_file")
        if subset_records is not None:
            for index, file_info in enumerate(file_infos, start=1):
                write_fastq_subset(
                    fixture_dir / file_info["name"],
                    fixture_dir / f"subset_{index}.fastq",
                    subset_records,
                    args.force,
                )
        metadata_path = fixture_dir / "fixture.json"
        metadata_path.write_text(json.dumps(fixture, indent=2, sort_keys=True) + "\n")
        print(f"meta   {metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
