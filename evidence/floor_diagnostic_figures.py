"""Compact publication-style views of the unchanged native diagnostic statistics."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy import stats
import figstyle as fs

BLUE='#0072B2'
ORANGE='#E69F00'
RED='#D55E00'


def canvas(title, label, size=(6.5,4.2)):
    fs.use_target('acs');fs.apply()
    fig, ax=plt.subplots(figsize=size)
    fig.suptitle(title, fontsize=11, weight='bold', y=.98)
    ax.set_title(label, fontsize=8, pad=12)
    return fig,ax


def save(fig,path):
    fig.tight_layout(rect=[0,0,1,.94])
    fig.savefig(path,dpi=600)
    fig.savefig(path.with_suffix('.pdf'))
    plt.close(fig)


def install(namespace):
    def scatter(df, metrics_df, output_path, bootstrap_iterations, analysis_label):
        x=df.weighted_coverage_pct.to_numpy(float)
        y=df.native_weighted_coverage_pct.to_numpy(float)
        metrics=metrics_df.set_index('metric')
        r=metrics.loc['pearson_r','estimate'];rho=metrics.loc['spearman_rho','estimate']
        fig,ax=canvas('Stage-3 and native-contact coverage',
            f'{analysis_label}  |  Pearson r = {r:.3f}; Spearman ρ = {rho:.3f}')
        grouped=df.groupby(['weighted_coverage_pct','native_weighted_coverage_pct']).size().reset_index(name='n')
        # Duplicate observations are encoded by area; a small legend makes this explicit.
        scale=12
        ax.scatter(grouped.weighted_coverage_pct,grouped.native_weighted_coverage_pct,
                   s=scale*grouped.n,color=BLUE,alpha=.5,edgecolors='none')
        fit=namespace['fit_linear_model'](x,y)
        grid=np.linspace(x.min(),x.max(),160)
        low,_,high=namespace['bootstrap_regression_band'](x,y,grid,
            iterations=bootstrap_iterations,random_seed=namespace['RANDOM_SEED'])
        ax.fill_between(grid,low,high,color=ORANGE,alpha=.2,label='95% bootstrap band')
        ax.plot(grid,fit['slope']*grid+fit['intercept'],color=RED,lw=1.2,ls='--',label='Linear fit')
        ax.set(xlabel='Stage-3 weighted coverage (%)',ylabel='Native-contact coverage (%)')
        ax.legend(fontsize=7,loc='best')
        fig.text(.5,.005,'Point area is proportional to the number of coincident ligands.',ha='center',fontsize=7)
        save(fig,output_path)

    def ranks(rank_df,rank_summary,output_path,analysis_label):
        fig,ax=canvas('Change in coverage ranking',
            f"{analysis_label}  |  Median |shift| = {rank_summary['median_abs_rank_shift']:.1f}; "
            f"top-{rank_summary['top_k']} overlap = {rank_summary['top_k_overlap_count']}")
        n=len(rank_df);span=max(1,float(rank_df.rank_shift.abs().max()))
        cmap=LinearSegmentedColormap.from_list('floor_signed',[ORANGE,'#F5F5F5',BLUE])
        points=ax.scatter(rank_df.screen_rank,rank_df.native_rank,c=rank_df.rank_shift,
            cmap=cmap,norm=TwoSlopeNorm(0,-span,span),s=14,edgecolor='#777777',linewidth=.15)
        ax.plot([1,n],[1,n],color='#666666',ls='--',lw=.8)
        ax.set(xlim=(.5,n+.5),ylim=(n+.5,.5),
               xlabel='Stage-3 coverage rank (1 = highest)',ylabel='Native coverage rank (1 = highest)')
        ax.set_aspect('equal',adjustable='box')
        fig.colorbar(points,ax=ax,pad=.025).set_label('Stage-3 rank − native rank',fontsize=8)
        fig.text(.5,.005,'Coverage ranks break score ties by ligand ID; full production ranks use additional criteria.',ha='center',fontsize=6.5)
        save(fig,output_path)

    def quartiles(df,quartile_summary,output_path,analysis_label):
        fig,ax=canvas('Native coverage across Stage-3 quartiles',analysis_label)
        groups=[df.loc[df.screen_quartile==r.quartile_index,'native_weighted_coverage_pct'].to_numpy()
                for r in quartile_summary.itertuples()]
        boxes=ax.boxplot(groups,patch_artist=True,widths=.5,showfliers=False,
            medianprops={'color':'#222222','linewidth':1.3},
            boxprops={'facecolor':'#D9EAF3','edgecolor':BLUE})
        rng=np.random.default_rng(20260310)
        for i,values in enumerate(groups,1):
            ax.scatter(i+rng.uniform(-.16,.16,len(values)),values,s=5,color=BLUE,alpha=.3,edgecolors='none')
        median=df.native_weighted_coverage_pct.median()
        ax.axhline(median,color=RED,ls='--',lw=1,label=f'Overall median: {median:.2f}%')
        ax.set_xticks(range(1,len(groups)+1),[f'Q{i}\n(n = {len(v):,})' for i,v in enumerate(groups,1)])
        ax.set(xlabel='Stage-3 weighted-coverage quartile (low to high)',ylabel='Native-contact coverage (%)')
        ax.legend(fontsize=7)
        save(fig,output_path)

    def residuals(df,output_path,analysis_label,residual_summary):
        fs.use_target('acs');fs.apply()
        fit=namespace['fit_linear_model'](df.weighted_coverage_pct.to_numpy(float),df.native_weighted_coverage_pct.to_numpy(float))
        fitted=np.asarray(fit['fitted']);residual=np.asarray(fit['residuals'])
        fig,axes=plt.subplots(1,2,figsize=(6.5,3.2))
        fig.suptitle('Linear-fit residual diagnostics',fontsize=11,weight='bold',y=.98)
        axes[0].scatter(fitted,residual,s=7,color=BLUE,alpha=.45,edgecolors='none')
        axes[0].axhline(0,color='#555555',ls='--',lw=.8)
        axes[0].set(xlabel='Fitted native coverage (%)',ylabel='Residual (percentage points)')
        axes[0].set_title(f"RMSE = {residual_summary['rmse']:.2f}; MAE = {residual_summary['mae']:.2f}",fontsize=8)
        (theoretical,observed),(slope,intercept,_)=stats.probplot(residual,dist='norm')
        axes[1].scatter(theoretical,observed,s=7,color=ORANGE,edgecolors='none')
        axes[1].plot(theoretical,slope*theoretical+intercept,color='#555555',ls='--',lw=.8)
        axes[1].set(xlabel='Theoretical normal quantile',ylabel='Ordered residual')
        axes[1].set_title(f'{len(df):,} native-scored ligands',fontsize=8)
        for ax,letter in zip(axes,'ab'): fs.panel_label(ax,f'({letter})',dx=-.18)
        save(fig,output_path)

    namespace.update(build_bubble_scatter=scatter,build_rank_shift_plot=ranks,
                     build_quartile_plot=quartiles,build_residual_diagnostics=residuals)
