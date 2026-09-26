"""Download the two GSE109828 processed temporal count matrices and metadata."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import time
from urllib.request import urlopen

import certifi
import ssl


def fetch(task):
    url, path = task
    partial = path.with_name(path.name + ".partial")
    if not path.exists():
        with urlopen(
            url, timeout=90, context=ssl.create_default_context(cafile=certifi.where())
        ) as response, partial.open("wb") as dest:
            size = int(response.headers.get("Content-Length", 0))
            written = 0
            last = time.monotonic()
            while True:
                block = response.read(4 * 1024 * 1024)
                if not block:
                    break
                dest.write(block)
                written += len(block)
                if time.monotonic() - last > 25:
                    print(
                        path.name,
                        round(written / 1e6),
                        "MB of",
                        round(size / 1e6),
                        flush=True,
                    )
                    last = time.monotonic()
        if size and written != size:
            raise ValueError("Incomplete download: " + str(path))
        partial.replace(path)
    h = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(8 * 1024 * 1024):
            h.update(block)
    print("Ready", path.name, path.stat().st_size, flush=True)
    return dict(
        file=path.name, url=url, bytes=path.stat().st_size, sha256=h.hexdigest()
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    tasks = []
    for gsm, exp in [("GSM2970930", "HSMM1"), ("GSM2970931", "HSMM2")]:
        for kind in ["indextable", "counts"]:
            name = f"{gsm}_sciATAC_{exp}_{kind}.txt.gz"
            url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2970nnn/{gsm}/suppl/{name}"
            tasks.append((url, a.out / name))
    tasks.append(
        (
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
            a.out / "hg38ToHg19.over.chain.gz",
        )
    )
    with ThreadPoolExecutor(max_workers=4) as executor:
        manifest = list(executor.map(fetch, tasks))
    (a.out / "download_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
