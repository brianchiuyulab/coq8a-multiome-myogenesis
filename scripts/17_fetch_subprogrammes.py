"""Freeze external muscle subprogrammes for intersection with the 221 genes."""

import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

SETS = [
    'GOBP_MUSCLE_CELL_FATE_COMMITMENT',
    'GOBP_MYOBLAST_DIFFERENTIATION',
    'GOBP_MYOTUBE_DIFFERENTIATION',
    'GOBP_POSITIVE_REGULATION_OF_MUSCLE_CELL_DIFFERENTIATION',
    'GOBP_MUSCLE_CELL_DIFFERENTIATION',
    'GOBP_MYOBLAST_FUSION',
    'GOBP_MYOFIBRIL_ASSEMBLY',
    'GOBP_MUSCLE_CONTRACTION',
    'DELASERNA_MYOD_TARGETS_UP',
    'DELASERNA_TARGETS_OF_MYOD_AND_SMARCA4',
    'VANOEVELEN_MYOGENESIS_SIN3A_TARGETS',
    'MEF2C_TARGET_GENES',
    'MEF2D_TARGET_GENES',
]


def fetch(name):
    url='https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName='+name+'&fileType=json'
    with urlopen(url, timeout=30) as response:
        raw=response.read()
    data=json.loads(raw)[name]
    assert data['geneSymbols'] and len(data['geneSymbols'])==len(set(data['geneSymbols']))
    data['download_url']=url
    data['download_sha256']=hashlib.sha256(raw).hexdigest()
    data['retrieved_utc']=datetime.now(timezone.utc).isoformat()
    print(name,len(data['geneSymbols']),flush=True)
    return name,data


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    with ThreadPoolExecutor(max_workers=4) as pool:
        data=dict(pool.map(fetch,SETS))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
