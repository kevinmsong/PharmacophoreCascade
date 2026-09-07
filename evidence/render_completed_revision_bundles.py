"""Render each genuine fresh native diagnostic bundle after its run completes."""
import argparse
import hashlib
import json
import time
from pathlib import Path
from render_floor_diagnostics import renderer,render_bundle

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'evidence/outputs'
META=OUT/'alert_disabled_revision'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--watch',action='store_true');args=ap.parse_args()
    run=renderer()
    module_sha=hashlib.sha256((ROOT/'evidence/floor_diagnostic_figures.py').read_bytes()).hexdigest()
    names=['benchmark_'+x for x in ['glp1r_full','ghsr_full','ntsr1_full','mdm2_full','glp1r_automated','glp1r_7ki0']]
    names += [f'{s}_seed{i}' for s in ['glp1r','ghsr','ntsr1','mdm2'] for i in range(1,6)]
    names += ['production_'+s for s in ['glp1r','ghsr','ntsr1','mdm2']]
    done=set()
    while True:
        pending=0
        bundles=[]
        for name in names:
            path=META/(name+'.json')
            if not path.exists():pending+=1;continue
            status=json.loads(path.read_text())['status']
            if status=='failed':raise RuntimeError(f'Failed run: {name}')
            if status!='completed':pending+=1;continue
            if name.startswith('production_'):
                system=name.removeprefix('production_')
                bundles += list((OUT/'absolute_floor'/system).glob('*/screening_production_native_terminal_bundle'))
            else:
                target=OUT/name if name.startswith('benchmark_') else OUT/'decoy_replicates'/name
                bundles.append(target/'benchmark_external/benchmark_native_bundle')
        head=ROOT/'results/absolute_floor_1000'
        if json.loads((head/'rerun_provenance.json').read_text())['status']=='completed':
            bundles.append(head/'screening_full_1M_floor1000_native_terminal_bundle')
        else:pending+=1
        for bundle in bundles:
            if bundle in done:continue
            manifest=bundle/'reports/figure_rendering_manifest.json'
            if manifest.exists():
                value=json.loads(manifest.read_text())
                if value.get('raster_dpi')==600 and value.get('figure_module_sha256')==module_sha:
                    done.add(bundle);continue
            render_bundle(bundle.resolve(),run)
            done.add(bundle)
        (META/'native_figure_queue.json').write_text(json.dumps({'pending_experiments':pending,
            'rendered_or_verified_bundles':[str(b.relative_to(ROOT)) for b in sorted(done)],
            'status':'completed' if pending==0 else 'waiting'},indent=2))
        if pending==0 or not args.watch:break
        time.sleep(30)
    print(f'Verified current 600-dpi/vector diagnostics in {len(done)} bundles',flush=True)


if __name__=='__main__':main()
