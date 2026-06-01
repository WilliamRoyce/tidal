"""Deepen-pass helper: download arXiv e-print source into ``literature/<id>/`` so the
verification agents can deep-read full text (READ-ONLY w.r.t. the manuscript).

Used for the Task-2 deepen pass: the keys whose specific claims were UNVERIFIABLE /
WEAK from the abstract alone, and which HAVE an arXiv eprint, get their full TeX
fetched here, re-tiering them to "A" (offline) for a deep re-verify.
"""

from __future__ import annotations

import gzip
import io
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

LIT = Path("literature")
UA = "tidal-bibaudit/1.0 (mailto:wr286@cam.ac.uk)"


def fetch_eprint(arxiv_id: str) -> str:
    safe = arxiv_id.replace("/", "_")  # old-style ids (gr-qc/0305049) contain '/'
    dest = LIT / safe
    if dest.exists() and any(dest.glob("*.tex")):
        return "cached"
    url = f"https://arxiv.org/e-print/{arxiv_id}"
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        raw = urllib.request.urlopen(req, timeout=60).read()
    except Exception as exc:
        return f"error: {exc}"
    dest.mkdir(parents=True, exist_ok=True)
    # arXiv source is usually a gzipped tar; sometimes a single gzipped .tex
    try:
        with tarfile.open(fileobj=io.BytesIO(raw)) as tf:
            n = 0
            for m in tf.getmembers():
                if m.isfile() and m.name.endswith((".tex", ".bbl")):
                    data = tf.extractfile(m).read()
                    (dest / Path(m.name).name).write_bytes(data)
                    n += 1
    except tarfile.ReadError:
        pass
    else:
        return f"tar: {n} tex/bbl files"
    try:
        text = gzip.decompress(raw)
        (dest / f"{safe}.tex").write_bytes(text)
    except Exception:
        (dest / f"{safe}.raw").write_bytes(raw)
        return "raw (unrecognised format)"
    else:
        return "single-gz tex"


def main(ids: list[str]) -> int:
    for i, aid in enumerate(ids, 1):
        res = fetch_eprint(aid.strip())
        print(f"[{i}/{len(ids)}] {aid}: {res}")
        if res != "cached":
            time.sleep(3.0)  # arXiv courtesy
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
