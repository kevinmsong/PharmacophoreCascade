"""Check the exact paired test against independent known null distributions."""
import json
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import analyze_equivalence as analysis


def main():
    original_root=analysis.ROOT
    results=[]
    with tempfile.TemporaryDirectory(prefix='docking_stat_check_') as folder:
        analysis.ROOT=Path(folder)
        target=analysis.ROOT/'evidence/outputs/revision/docking_top10_summary.csv'
        target.parent.mkdir(parents=True)
        def run(delta):
            delta=np.asarray(delta,dtype=float)
            pd.DataFrame({'final_rank':range(1,len(delta)+1),'ligand_id':[f'L{i}' for i in range(len(delta))],
                'best_active':delta-10,'best_inactive':np.full(len(delta),-10.),'active_pref':delta}).to_csv(target,index=False)
            return analysis.run_docking_test()
        for delta in [[-1,-2,-3,-4,-5,-6,-7,-8,-9,-10],[-1,2,-3,4,-5,6,-7,8,-9,10]]:
            observed=run(delta);expected=wilcoxon(delta,method='exact')
            assert observed['wilcoxon_W']==expected.statistic
            assert abs(observed['p_two_sided']-expected.pvalue)<1e-14
            results.append({'case':'untied vs SciPy exact','passed':True})
        tied=run([-1,-1,-1,-1])
        # Four equal ranks: only all-negative and all-positive are as extreme.
        assert tied['wilcoxon_W']==0 and tied['p_two_sided']==2/16
        results.append({'case':'equal-rank exact null, 2/16','passed':True})
        zeros=run([0,0,0,0]);assert zeros['p_two_sided']==1 and zeros['n_zero_differences']==4
        missing=run([-1,-2,0,np.nan]);assert missing['n']==3 and missing['n_selected']==4
        assert missing['p_two_sided']==.5 and missing['n_zero_differences']==1
        results.extend([{'case':'all zero differences','passed':True},{'case':'missing pair and zero excluded','passed':True}])
    analysis.ROOT=original_root
    output=original_root/'evidence/outputs/alert_disabled_revision/docking_stat_validation.json'
    output.write_text(json.dumps(results,indent=2))
    print(f'Passed {len(results)} exact signed-rank checks')


if __name__=='__main__':main()
