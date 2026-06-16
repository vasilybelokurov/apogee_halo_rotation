#!/usr/bin/env python3
"""Download external catalogues used by the APOGEE halo rotation analysis."""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.request import urlopen


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "external"

DATASETS = {
    "allstar": {
        "url": "https://data.sdss.org/sas/dr17/apogee/spectro/aspcap/dr17/synspec_rev1/allStarLite-dr17-synspec_rev1.fits",
        "path": DATA_DIR / "allStarLite-dr17-synspec_rev1.fits",
    },
    "astronn": {
        "url": "https://data.sdss.org/sas/dr17/apogee/vac/apogee-astronn/apogee_astroNN-DR17.fits",
        "path": DATA_DIR / "apogee_astroNN-DR17.fits",
    },
    "vasiliev_gc_members": {
        "url": "https://zenodo.org/records/4891252/files/clusters.zip?download=1",
        "path": DATA_DIR / "vasiliev_baumgardt_2021_clusters.zip",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--only",
        choices=sorted(DATASETS),
        nargs="+",
        default=sorted(DATASETS),
        help="Datasets to download. Defaults to all.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Replace files that already exist.")
    return parser.parse_args()


def download(url: str, path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        print(f"Exists: {path}")
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    print(f"Downloading {url}")
    print(f"       to {path}")
    with urlopen(url) as response, tmp.open("wb") as out:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(path)


def main() -> None:
    args = parse_args()
    for name in args.only:
        dataset = DATASETS[name]
        download(dataset["url"], dataset["path"], overwrite=args.overwrite)


if __name__ == "__main__":
    main()
