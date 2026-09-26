"""Retrieve the GSE240061 MINiML sample manifest from GEO."""

from pathlib import Path
import argparse
import hashlib
import io
import json
import tarfile
import xml.etree.ElementTree as ET

import pandas as pd
import requests

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=True)
url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240061/miniml/GSE240061_family.xml.tgz"
response = requests.get(url, timeout=90)
response.raise_for_status()
with tarfile.open(fileobj=io.BytesIO(response.content), mode="r:gz") as archive:
    names = [x for x in archive.getnames() if x.endswith("_family.xml")]
    assert len(names) == 1
    xml = archive.extractfile(names[0]).read()
tree = ET.fromstring(xml)
ns = {"g": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}
rows = []
for sample in tree.findall("g:Sample", ns):
    row = dict(
        gsm=sample.attrib["iid"], title=sample.findtext("g:Title", namespaces=ns)
    )
    row["source"] = sample.findtext("g:Channel/g:Source", namespaces=ns)
    row["organism"] = sample.findtext("g:Channel/g:Organism", namespaces=ns)
    for item in sample.findall("g:Channel/g:Characteristics", ns):
        row[item.attrib.get("tag", "characteristic")] = (item.text or "").strip()
    row["supplementary_files"] = ";".join(
        x.text for x in sample.findall("g:Supplementary-Data", ns) if x.text
    )
    rows.append(row)
pd.DataFrame(rows).to_csv(
    args.out / "GEO_sample_manifest.tsv", sep="\t", index=False, na_rep="NA"
)
(args.out / "GEO_manifest_source.json").write_text(
    json.dumps(
        dict(
            url=url,
            sha256=hashlib.sha256(response.content).hexdigest(),
            n_assay_accessions=len(rows),
        ),
        indent=2,
    )
    + "\n"
)
print(pd.DataFrame(rows).drop(columns="supplementary_files").to_string(index=False))
