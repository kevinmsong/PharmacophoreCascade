"""Incorporate completed, validated floor results into the shared manuscript sources."""
import json
import re
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SUB=ROOT/'ACS_Omega_resubmission'
DATA=ROOT/'evidence/outputs/absolute_floor'
HEAD=ROOT/'results/absolute_floor_1000'


def replace_once(text,old,new):
    if text.count(old)!=1: raise ValueError(f'Expected exactly one source passage: {old[:90]}')
    return text.replace(old,new,1)


def main():
    raise RuntimeError('Superseded by the user-directed alert-disabled revision. All benchmark claims must be regenerated before manuscript integration.')
    summary=json.loads((DATA/'analysis_summary.json').read_text())
    head=summary['headline']
    run=json.loads((HEAD/'screening_full_1M_floor1000_run_summary.json').read_text())
    rows=pd.DataFrame(summary['benchmarks']).set_index(['system','policy'])
    def value(sys,policy,col): return int(rows.loc[(sys,policy),col])
    def pair(sys):
        a=value(sys,'percentage','final_ranked_actives');b=value(sys,'floor1000','final_ranked_actives')
        n=value(sys,'percentage','n_actives')
        return a,b,n
    ga,gb,gn=pair('ghsr');na,nb,nn=pair('ntsr1')
    nc=head['rerun_native_counts'];counts=head['rerun_counts']
    hours=head['pipeline_wall_seconds']/3600
    # Primary scientific comparisons must remain identifiable.
    if not head['shortlist_ids_same_order_as_archive']:
        raise RuntimeError('Headline shortlist does not reproduce archive; do not write attribution to the floor')
    if not head['final_ids_same_order_as_archive'] or head['native_score_max_abs_delta_common_final']>1e-9:
        raise RuntimeError('Headline ranking changed: update cohort/docking/correlations before revising the manuscript')
    for name, check in head['full_cohort_checks'].items():
        if not check['same_order'] or max(check['max_abs_deltas'].values())>1e-9:
            raise RuntimeError(f'{name} differs from archive: audit downstream supplementary claims before integration')
    body=(SUB/'manuscript_body.tex').read_text(encoding='utf-8')
    abstract=rf'''\noindent \new{{Peptide-receptor interfaces can be difficult to screen at library scale because
binding determinants extend across broad surfaces. We evaluate a staged workflow
ending in a pharmacophore built from the bound peptide's receptor contacts and ask
which component supplies the enrichment. Across GLP-1R, GHSR, NTSR1, and MDM2--p53,
native peptide-contact scoring exceeds a conventional single-pass 3D pharmacophore
by 0.08--0.66 ROC-AUC; a native-only baseline matches or exceeds the full cascade on
three systems. We then test an actionable shortlist policy: retain the larger of
5\% or 1,000 molecules, capped by eligible candidates. A fresh million-compound
GLP-1R screen produces {counts['shortlist_size']:,} shortlisted molecules,
{nc['successful_best_ligands']:,} successfully native-scored ligands, and the same ordered
top-1,000 as the archived screen; the floor is inactive at that scale. In paired
production-constrained benchmarks, the floor increases final active retention from
{ga}/{gn} to {gb}/{gn} for GHSR and from {na}/{nn} to {nb}/{nn} for NTSR1. GLP-1R and
MDM2--p53 have identical shortlists under both policies. The floor prevents an
additional shortlist cut on small pools, but does not reverse upstream gate losses
or guarantee survival through native preparation and scaffold selection. Native
contacts supply the enrichment; a minimum shortlist makes the efficiency-retention
trade more controllable. Per-molecule rerun outputs, configurations, and a reference
implementation accompany these computational results, which do not establish
biological activity.}}'''
    body,n=re.subn(r'(?<=\\begin\{abstract\}\n).*?(?=\n\\end\{abstract\})',lambda m:abstract,body,count=1,flags=re.S)
    if n!=1: raise ValueError('Abstract not found')
    body=body.replace('Peptide-mediated receptor activation and protein-protein recognition are important\nbut difficult targets for small-molecule discovery.',
                      'Peptide-mediated receptor activation and protein-protein recognition are important,\nbut difficult targets, for small-molecule discovery.')
    body=replace_once(body,'that survive the cheap ones. Survivors are ranked by cascade score, the top fraction\nis shortlisted into Stage 3 for conformer-level geometric reranking, and the native',
        'that survive the cheap ones. Survivors are ranked by cascade score; the revised\nshortlist retains the top fraction with a minimum of 1,000 eligible molecules for\nStage-3 conformer-level geometric reranking, and the native')
    start=body.index(r'\subsection{Million-scale tractability}')
    end=body.index(r'\subsection{\new{The terminal native stage supplies the enrichment}}',start)
    body=body[:start]+rf'''\subsection{{Million-scale tractability}}
\new{{We reran the complete million-compound GLP-1R screen from the original H17--H20
ZINC tranche files, with a 1,000-molecule shortlist floor and no reuse of archived
scores or conformers.\cite{{Irwin2020}} The run reproduces {counts['property_pass']:,}
Stage-0 survivors, {counts['hotspot_pass']:,} Stage-1/2 candidates, and the ordered
{counts['shortlist_size']:,}-molecule shortlist. At this scale the 5\% term exceeds
1,000, so adding the floor changes neither shortlist membership nor its order.
Fresh Stage-3 scoring yields {counts['final_hits']:,} successful ligand scores.
The native branch assembles {nc['candidate_pool_size']:,} candidates, selects
{nc['selected_ligands']:,} with an eight-per-scaffold cap,\cite{{Bemis1996}}
prepares {nc['prepared_microstates']:,} states, and successfully scores
{nc['scored_microstates']:,} states from {nc['successful_best_ligands']:,} ligands.
The ordered final top-1,000 agrees with the archived ranking (Tables~S1 and S12).
This is a demonstrated nonbinding-floor result, not a prediction from shortlist
size alone.}}

\new{{The fresh run took {head['pipeline_wall_seconds']:,.1f}~s ({hours:.2f}~h) end to end
using 12 process workers for the prescreen, Stage 3, and native calculations on the
same workstation. The total includes native preparation and scoring once.
The earlier 13.34~h statement added native-branch time to a pipeline total that
already included it; that double counting is corrected here and in Table~S1.
We report the measured wall time and numbers scored, without treating an
extrapolated full-library native cost as a measured speedup.}}

''' + body[end:]
    start=body.index(r'\subsection{\new{The trade between efficiency and active retention}}')
    end=body.index(r'\subsection{',start+15)
    efficiency=rf'''\subsection{{\new{{The trade between efficiency and active retention}}}}
\new{{The original fraction sweep identified a potentially avoidable loss: after
the Stage-1 gate, the GHSR and NTSR1 production libraries contain only 267 and 384
candidates. A further 5\% cut admits only 14 and 20. We tested the proposed repair
with fresh, paired production-constrained runs rather than equating shortlist
admission with successful native scoring (\figref{{fig:efficiency}};
\tabref{{tab:efficiency}}). Both arms use identical gates, scores, conformer settings,
native-pool quotas, and scaffold caps; they differ only in the shortlist minimum.}}

\new{{With a 1,000-molecule floor, all 267 GHSR and 384 NTSR1 gate survivors are
admitted to Stage 3. GHSR final active retention changes from {ga}/{gn}
({100*ga/gn:.0f}\%) to {gb}/{gn} ({100*gb/gn:.0f}\%); NTSR1 changes from {na}/{nn}
({100*na/nn:.0f}\%) to {nb}/{nn} ({100*nb/nn:.0f}\%). These are final-ranked counts after
3D scoring, native preparation, scaffold selection, and native scoring, not
shortlist-based upper bounds. The GLP-1R and MDM2--p53 5\% shortlists already
contain {value('glp1r','percentage','shortlist_n'):,} and {value('mdm2','percentage','shortlist_n'):,}
molecules, so the floor changes neither input nor output for those pairs.
Figure~S5 and Table~S13 account for every subsequent stage.}}

\new{{The floor is therefore a minimum on the percentage shortlist, not a fixed
1,000-molecule cap and not a guarantee that all gate-surviving actives reach the
final ranking. It prevents a second severe cut when the eligible pool is small;
upstream exclusions and downstream preparation or diversity constraints remain.
The native-branch wall times are measured and reported in Table~S14. The earlier
claims of universally complete retention and negligible or seconds-long cost are
replaced by these observed counts and timings.}}

\begin{{figure}}[tbp]
\centering
\includegraphics[width=\textwidth]{{fig6_efficiency.pdf}}
\caption{{\new{{\textbf{{The trade between efficiency and active retention under a
measured shortlist floor.}} (a) Molecules actually admitted to Stage 3 by each
policy. (b) In-domain actives in the final native ranking, with counts shown
directly. Orange circles denote 5\% only; blue squares denote 5\% with a minimum
of 1,000 eligible molecules. Coincident results are offset vertically for
visibility. Identical large-pool shortlists share one native execution; changed
small-pool shortlists undergo separate native preparation and scoring.}}}}
\label{{fig:efficiency}}
\end{{figure}}

\input{{efficiency_table.tex}}

'''
    body=body[:start]+efficiency+body[end:]
    body=body.replace('The dominant finding is that \\emph{every} Stage-0 active loss\nacross all four systems is a structural-alert rejection rather than a\nproperty-envelope failure:',
        'Most Stage-0 active losses\nacross the four systems are structural-alert rejections rather than\nproperty-envelope failures:')
    body=body.replace('Rerunning each system under the production settings (Methods), final active',
        'In the archived percentage-only production runs (Methods), final active')
    body=body.replace('of the 102-feature receptor-side model for Stages 1 and 2, shortlisted the top 5.0\\%\nof Stage-0 passers for Stage 3, and generated 16 conformers per shortlisted ligand',
        'of the 102-feature receptor-side model for Stages 1 and 2, shortlisted the top 5.0\\%\nof Stage-0 passers with a minimum of 1,000 eligible molecules for Stage 3, and\ngenerated 16 conformers per shortlisted ligand')
    marker=r'\subsection{\new{Production-constrained screening and orthogonal docking}}'
    methods=rf'''\new{{For the floor experiment, the original input identities, SMILES, and labels
were fixed, but all gate evaluations and molecular scores were recomputed.
The shortlist size was $K=\min(N,\max(\lceil0.05B\rceil,1000))$, where $N$ is the
eligible Stage-1/2 pool. To preserve each archived experiment, $B$ is the Stage-0
survivor count for the headline screen and the Stage-1/2 candidate count for the
production benchmarks. The headline uses hotspot/pair weights 0.4/0.6; the
archived production benchmarks use 0.25/0.75. Both use fixed per-type feature
caps (anion, cation, donor, acceptor, aromatic, hydrophobe: 2, 2, 4, 6, 6, 6).
These settings are pinned explicitly rather than inherited from later engine
defaults. Gate decisions and score columns are checked against archived outputs
before expensive reranking. Within each policy pair, Stages 0--2 and the Stage-3
union are freshly computed once; identical shortlists share downstream outputs,
whereas different shortlists undergo separate native preparation and scoring.
Full feature records, rather than zero-filled placeholders, are carried into
reranking and pool selection. Native preparation and scoring may be distributed
over ordered process chunks; serial/parallel validation gave identical prepared
SDFs and scores, mappings, and coordinates agreeing to $10^{{-12}}$ on a
51-microstate check. Per-molecule survival and execution provenance are archived.}}

'''
    body=replace_once(body,marker,marker+'\n'+methods)
    body=body.replace('shortlist fraction from 0.5\\% to 100\\% at fixed gate settings, recording molecules\nreaching native scoring and in-domain actives retained at each setting.',
        'shortlist fraction from 0.5\\% to 100\\% at fixed gate settings, recording\nshortlist membership; these values are admission counts, not measured native\nscoring outcomes. The new paired floor experiment measures the latter directly.')
    start=body.index(r'\new{The staged architecture should accordingly be judged as an efficiency mechanism.')
    end=body.index('\n\n',start)
    body=body[:start]+rf'''\new{{The staged architecture should accordingly be judged as an efficiency mechanism.
The fresh million-compound run supplies {nc['successful_best_ligands']:,} successful native ligand scores
from {nc['selected_ligands']:,} selected candidates in {hours:.2f}~h. At that scale the
minimum does not bind and the archived top-1,000 ranking is reproduced. The floor
does matter when a selective upstream gate leaves only a few hundred molecules:
the paired experiments improve final active retention for GHSR and NTSR1 by
avoiding another severe shortlist cut. That improvement is conditional on the
eligible pool and does not undo upstream losses or downstream selection limits.}}'''+body[end:]
    body=body.replace('Read as guidance rather than advocacy, these results suggest a different default configuration than the one we ran.',
        'These results support a minimum-protected percentage shortlist, now tested directly.')
    body=body.replace('Read as guidance rather than advocacy, these results suggest a different default\nconfiguration than the one we ran.',
        'These results support a minimum-protected percentage shortlist, now tested directly.')
    body=body.replace('Set the shortlist as an absolute floor rather than a percentage, so a selective\nupstream gate cannot compound into a near-empty pool.',
        'Retain the percentage shortlist with a minimum of 1,000 eligible candidates,\ncapped by the available pool, so a selective upstream gate is not followed by\nanother unnecessary cut to a near-empty shortlist.')
    body=body.replace('13.3~h','the measured rerun time')
    # Remove obsolete future-work wording wherever it might survive.
    if 'We did not re-run the headline screen' in body: raise RuntimeError('Obsolete claim survived revision')
    si=update_si(run,head)
    response=update_response(ga,gb,gn,na,nb,nn,counts['shortlist_size'])
    # Resolve every source anchor before replacing any manuscript file.
    for name, text in [('manuscript_body.tex',body),('supporting_information.tex',si),
                       ('response_to_reviewers.tex',response)]:
        (SUB/name).write_text(text,encoding='utf-8')
    print('Updated shared manuscript, abstract, Methods, Discussion, SI, and response sources')


def update_si(run,head):
    path=SUB/'supporting_information.tex';text=path.read_text(encoding='utf-8')
    text=text.replace('Figure S1 through Figure S4','Figure S1 through Figure S5').replace('Table S1 through Table S11','Table S1 through Table S14')
    text=replace_once(text, r'\item Figure S4. Structures of the top 20 native-ranked ligands.',
        r'\item Figure S4. Structures of the top 20 native-ranked ligands.'+'\n'+
        r'\item Figure S5. Active survival in fresh paired absolute-floor experiments.')
    text=replace_once(text, r'\item Table S11. Stage-by-stage active survival under production constraints.',
        r'\item Table S11. Archived percentage-only active survival under production constraints.'+'\n'+
        r'\item Table S12. Fresh million-compound absolute-floor rerun versus the archive.'+'\n'+
        r'\item Table S13. Stage-by-stage counts in paired absolute-floor experiments.'+'\n'+
        r'\item Table S14. Measured native-branch wall times for both policies.')
    text=text.replace('Workers & 15', 'Workers & 12 process workers for prescreen, Stage 3, and native branch')
    text=text.replace('Shortlist fraction & Top 5.0\\% of Stage-0-pass molecules',
        'Shortlist rule & Top 5.0\\% of Stage-0-pass molecules, minimum 1,000, capped by eligible candidates')
    t=run['timings'];nt=run['native_rerank']['timings'];c=run['counts'];nc=run['native_rerank']['counts']
    replacements={
        'Stage-0 property-pass count':f"{c['property_pass']:,}", 'Stage-1 hotspot-pass count':f"{c['hotspot_pass']:,}",
        'Stage-3 shortlist size':f"{c['shortlist_size']:,}", 'Stage-3 final hits written':f"{c['final_hits']:,}",
        'Native candidate-pool size':f"{nc['candidate_pool_size']:,}", 'Native selected ligands':f"{nc['selected_ligands']:,}",
        'Prepared states':f"{nc['prepared_microstates']:,}", 'Successfully scored states':f"{nc['scored_microstates']:,}",
        'Successfully native-scored ligands':f"{nc['successful_best_ligands']:,}",
        'Stage-0 time':f"{t['stage0_property_gate_sec']:,.1f}~s (summed worker time)",
        'Stage-1 time':f"{t['stage1_hotspot_scoring_sec']:,.1f}~s (summed worker time)",
        'Stage-2 time':f"{t['stage2_pair_hash_scoring_sec']:,.1f}~s (summed worker time)",
        'Stage-3 time':f"{t['stage3_rerank_sec']:,.1f}~s (wall time)",
        'Native branch selection time':f"{nt['selection_sec']:,.1f}~s",
        'Native branch preparation time':f"{nt['library_prepare_sec']:,.1f}~s",
        'Native branch scoring time':f"{nt['scoring_sec']:,.1f}~s",
        'Native branch wall time':f"{nt['total_sec']:,.1f}~s (included in total)"}
    for field,value in replacements.items():
        pattern=re.escape(field)+r' & [^\n]+'
        text,n=re.subn(pattern,lambda m:field+' & '+value+r' \\',text,count=1)
        if n!=1: raise ValueError(f'SI row missing: {field}')
    text=re.sub(r'Stage-0 through Stage-3 wall time & [^\n]+\n','',text,count=1)
    text=re.sub(r'\\textbf\{Combined total wall time\} & [^\n]+',
        lambda m:rf"\textbf{{Total pipeline wall time}} & \textbf{{{t['total_pipeline_sec']:,.1f}~s ({t['total_pipeline_sec']/3600:.2f}~h)}} \\",text,count=1)
    text=re.sub(r'Estimated cost of native-scoring all 1,000,000 & [^\n]+',
        lambda m:r'Full-library native-only runtime & Not measured; no extrapolated speedup asserted \\',text,count=1)
    text=text.replace('background and screened with the million-compound-run settings, so that active',
        'background and screened with the archived percentage-only production settings, so that active')
    extra=r'''
\clearpage
\section*{Measured Absolute-Floor Reruns}

The complete million-compound screen was recomputed from its original ZINC inputs
with a minimum shortlist size of 1,000. The original fixed feature caps,
0.4/0.6 hotspot/pair weights, Stage-0 denominator, conformer settings, and native
selection quotas were pinned explicitly. No archived scores or conformers were
used in this execution. Table~S12 compares the resulting counts with the archived
headline run. The floor does not bind at this scale.

\input{floor_headline_table.tex}

\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{figA5_floor_survival.pdf}
\caption{Stage-by-stage active retention in fresh paired production experiments.
Orange circles and dashed lines denote 5\% only; blue squares and solid lines
denote 5\% plus a 1,000-molecule floor capped by the eligible pool. The floor
affects shortlist membership, whereas the subsequent stages measure actual
preparation, scoring, and selection outcomes. Identical curves may overlap.
``Selected'' and ``Scored'' refer to the native branch.}
\label{fig:floor-survival}
\end{figure}

The paired benchmark inputs preserve the archived labels and 30,000-molecule
backgrounds. Their 0.25/0.75 hotspot/pair weights and candidate-count denominator
are retained in both arms. Fresh gate outcomes and score columns are checked
against the archived production evaluations. Full Stage-0/2 records are carried
into Stage 3 and pool selection; the earlier benchmark runner used zero-filled
placeholders for several of those fields. Consequently the paired current
executions, rather than an old-to-new difference, define the floor effect.
Shortlisting, successful native scoring, and final ranking are reported separately.
The archived attrition tables remain identified as percentage-only results.

\input{floor_survival_table.tex}
\input{floor_timing_table.tex}

The machine-readable release includes input identities and labels; all fresh
Stage-0/2 evaluations; Stage-3 inputs and scores; per-policy native selections,
scores, and final ranks; per-molecule stage survival; run logs; and SHA-256
provenance manifests. Shared computation between identical policy inputs is
declared explicitly. The earlier total runtime double-counted the native branch;
Table~S1 now gives one end-to-end total for the fresh rerun.

'''
    text=replace_once(text,r'\end{document}',extra+r'\end{document}')
    return text


def update_response(ga,gb,gn,na,nb,nn,k):
    path=SUB/'response_to_reviewers.tex';text=path.read_text(encoding='utf-8')
    start=text.index('We swept the Stage-3 shortlist fraction')
    end=text.index('Separately, and bearing directly',start)
    paragraph=rf'''We have now performed the full absolute-floor test, including a fresh rerun of
the headline million-compound screen. The tested rule retains the larger of 5\%
and 1,000 molecules, capped by eligible candidates. Archived score weights,
feature caps, and percentage denominators are pinned explicitly, and molecular
score/identity checks guard against later engine-default changes.

In paired production-constrained runs, GHSR final active retention changes from
{ga}/{gn} to {gb}/{gn}, and NTSR1 from {na}/{nn} to {nb}/{nn}. The floor expands their
shortlists from 14 to 267 and from 20 to 384 molecules. GLP-1R and MDM2--p53 have
identical shortlists under the two policies. These are measured final-ranked
counts after preparation, 3D scoring, native scoring, and scaffold selection,
not counts inferred from a shortlist sweep. In the million-compound rerun,
the floor is inactive: the {k:,}-molecule shortlist and ordered final top-1,000
are reproduced.

The revised abstract, Results, Methods, Discussion, Figure~5, and Table~5 now
report this completed experiment. Figure~S5 and Tables~S12--S14 supply stage counts,
headline reproduction, and measured native-branch wall times. The graphical
abstract has been redesigned around the measured floor results. We removed the
unsupported statements that all four systems would retain every gate-surviving
active and that the extra scoring would take seconds or have negligible cost.
The floor prevents a further severe cut on small pools; it does not reverse
upstream gate losses or guarantee survival of native/scaffold constraints.

We also corrected an independent timing-accounting error: the previously stated
13.34~h added native-branch time to a pipeline total that already included it.
The revised manuscript reports a single measured end-to-end time for the fresh
rerun and does not present an extrapolated full-library cost as a measured speedup.

'''
    text=text[:start]+paragraph+text[end:]
    return text


if __name__=='__main__': main()
