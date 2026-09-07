"""Fresh million-compound rerun, explicitly pinned to the archived scientific settings.

No cached scores or prepared conformers are supplied to the screening engine.
The archived headline used a Stage-0 percentage denominator and 0.4/0.6 weights;
the later engine defaults differ. These choices are fixed explicitly here.
"""
import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "absolute_floor_1000"


def sha(path):
    h = hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(4*1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers",type=int,default=12)
    args=ap.parse_args()
    OUT.mkdir(exist_ok=True,parents=True)
    manifest_path=OUT/"rerun_provenance.json"
    if manifest_path.exists():
        raise SystemExit("Run folder already has provenance; preserve it and use a new folder for another run.")
    command=[sys.executable,"-u","run_optimized_1M_topological_hashed_screening.py",
        "--total-target","1000000","--shortlist-pct","5","--shortlist-floor","1000",
        "--shortlist-basis","stage0","--cascade-hotspot-weight","0.4",
        "--feature-cap-mode","fixed","--native-workers",str(args.workers),
        "--chemistry-gate-mode","warn_only",
        "--workers",str(args.workers),"--prescreen-workers",str(args.workers),
        "--results-dir",str(OUT),"--output-label","full_1M_floor1000"]
    files=[ROOT/"run_optimized_1M_topological_hashed_screening.py",
        ROOT/"evidence/parallel_native.py",
        Path(__file__),ROOT/"maps/pharmacophore_rigorous.json",
        ROOT/"GLP1_top_ligand_analysis/configs/study_top1000_full_1M_native_tuned_no_docking.yaml",
        ROOT/"GLP1_top_ligand_analysis/6X18_GLP1_GLP1R.pdb"]
    files+= sorted(p for p in (ROOT/"tmp/zinc").glob("*.smi.gz") if p.name.startswith(("H17","H18","H19","H20")))
    files+=sorted((ROOT/"GLP1_top_ligand_analysis/src/glp1r_state_preference").glob("*.py"))
    import rdkit,numpy,pandas
    manifest={"status":"running","started_utc":datetime.now(timezone.utc).isoformat(),
        "command":command,"cwd":str(ROOT),"python":sys.version,"platform":platform.platform(),
        "rdkit":rdkit.__version__,"numpy":numpy.__version__,"pandas":pandas.__version__,
        "floor_rule":"min(N_eligible, max(ceil(0.05*N_stage0), 1000))",
        "baseline":"same original input library, Stage-0 denominator and 0.4/0.6 weights; structural-alert exclusions now disabled",
        "chemistry_gate_mode":"warn_only", "native_pains_filter":False,
        "feature_caps":"fixed (2,2,4,6,6,6), matching archived scores",
        "execution":"ordered process parallelism for prescreen, Stage 3, and native branch; unchanged per-ligand routines",
        "cache_reuse":False,"inputs_sha256":{str(p.relative_to(ROOT)):sha(p) for p in files}}
    manifest_path.write_text(json.dumps(manifest,indent=2))
    started=time.perf_counter()
    with open(OUT/"rerun.log","w",encoding="utf-8") as log:
        result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    manifest.update(status="completed" if result.returncode==0 else "failed",
        returncode=result.returncode,elapsed_seconds=time.perf_counter()-started,
        completed_utc=datetime.now(timezone.utc).isoformat())
    manifest_path.write_text(json.dumps(manifest,indent=2))
    raise SystemExit(result.returncode)


if __name__=="__main__":
    main()
