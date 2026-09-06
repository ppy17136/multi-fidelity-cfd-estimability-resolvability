"""Figures from frozen data; no simulation or benchmark times are regenerated."""
from pathlib import Path
import csv, itertools, json, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
DATA=HERE
OUT=HERE/'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':12,'axes.titlesize':12,'axes.labelsize':12,'legend.fontsize':10,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
COLORS=['#0072B2','#D55E00','#009E73','#CC79A7']

def rows(name):
    with (DATA/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def save(fig,name):
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(OUT/(name+'.png'),bbox_inches='tight',dpi=220)
    plt.close(fig)

def example():
    z=np.array([[5,1,4],[-1,-4,0],[-1,2,5]])
    parts=[[[0,1,2]],[[0],[1,2]],[[0,1],[2]],[[0],[1],[2]]]
    result=[]
    for p in parts:
        vals=[]
        for labels in itertools.product(range(3),repeat=len(p)):
            vals.append(int(sum(z[i,l] for block,l in zip(p,labels) for i in block)))
        result.append({'partition':p,'values':sorted(set(vals)),'lower':min(vals),'upper':max(vals),'width':max(vals)-min(vals)})
    assert [(x['lower'],x['upper']) for x in result]==[(-1,9),(-1,10),(-4,9),(-4,10)]
    assert max(result[1]['width'],result[2]['width'])==13
    assert max(result[1]['upper'],result[2]['upper'])-min(result[1]['lower'],result[2]['lower'])==14
    fig,axs=plt.subplots(1,3,figsize=(10,4.2),gridspec_kw={'width_ratios':[1,1.45,1.1]})
    ax=axs[0];ax.set_title('(a) Fixed contributions and graph',loc='left');ax.axis('off')
    table=ax.table(cellText=z,colLabels=['Path 1','Path 2','Path 3'],rowLabels=['Cell 1','Cell 2','Cell 3'],loc='upper center',bbox=[.08,.48,.9,.43]);table.auto_set_font_size(False);table.set_fontsize(10)
    for j in range(2):
        ax.plot([.18+j*.31,.49+j*.31],[.23,.23],color=COLORS[0],lw=2,transform=ax.transAxes)
        ax.text(.335+j*.31,.3,'0.5',ha='center',transform=ax.transAxes)
    for j in range(3):
        ax.scatter(.18+j*.31,.23,s=430,color='white',edgecolors=COLORS[0],zorder=3,transform=ax.transAxes)
        ax.text(.18+j*.31,.23,str(j+1),va='center',ha='center',transform=ax.transAxes)
    ax=axs[1];ax.set_title('(b) Ranges under one partition',loc='left')
    for j,r in enumerate(result):
        y=3-j;ax.hlines(y,r['lower'],r['upper'],color=COLORS[j],lw=3)
        ax.plot(r['values'],[y]*len(r['values']),'o',color=COLORS[j],ms=4)
        ax.text(10.65,y,f"width {r['width']}",va='center',fontsize=9)
    ax.axvline(-2,ls='--',color='.35',lw=1.3,label=r'$\tau=-2$')
    ax.set(yticks=range(4),yticklabels=['1 | 2 | 3','12 | 3','1 | 23','123'],xlabel='Scalar contrast',xlim=(-5,14.2),ylim=(-.5,3.7));ax.legend(loc='upper left')
    ax=axs[2];ax.set_title('(c) Distinct certificate targets',loc='left')
    ax.step([0,.5,1],[10,13,14],where='post',lw=2,color=COLORS[0],label='Common-partition width')
    ax.step([0,.5,1],[10,14,14],where='post',lw=2,ls='--',color=COLORS[1],label='Union width')
    ax.set(xlabel='Relaxation budget',ylabel='Width',xticks=[0,.5,1],ylim=(9.5,15.9))
    ax.text(.04,15.5,r'Decision loss: $\eta_{-2}=0.5$',fontsize=9)
    ax.text(.04,14.85,r'Full common width: $\rho_u=1$',fontsize=9)
    ax.legend(loc='lower right',fontsize=8)
    fig.tight_layout(w_pad=2.2);save(fig,'expanded_worked_example')
    (HERE/'worked_example_verified.json').write_text(json.dumps({'contributions':z.tolist(),'edge_costs':[.5,.5],'threshold':-2,'eta':.5,'rho':1,'partitions':result},indent=2))

def performance():
    f=rows('32_FOREST_DP_SCALING_BENCHMARK.csv');g=rows('30_GRAPH_COUPLING_MILP_BENCHMARK.csv')
    fig,axs=plt.subplots(1,3,figsize=(9,4.2))
    for j,top in enumerate(['chain','random_tree']):
        rr=[r for r in f if r['topology']==top];ns=sorted({int(r['cells']) for r in rr})
        for ax,key in zip(axs[:2],['runtime_s','peak_tracemalloc_mib']):
            a=[np.array([float(r[key]) for r in rr if int(r['cells'])==n]) for n in ns]
            ax.loglog(ns,[np.median(x) for x in a],'o-',color=COLORS[j],label=top.replace('_',' '),lw=1.7)
            ax.fill_between(ns,[x.min() for x in a],[x.max() for x in a],alpha=.12,color=COLORS[j])
    axs[0].set(title='(a) Forest runtime',xlabel='Vertices',ylabel='Seconds (median; range shaded)');axs[0].legend()
    axs[1].set(title='(b) Traced Python allocation',xlabel='Vertices',ylabel='Peak allocation (MiB)')
    times=np.sort([float(r['runtime_s']) for r in g]);nodes=np.array([float(r['branch_nodes']) for r in g]);assert len(g)==80
    axs[2].semilogx(times,np.arange(1,len(times)+1)/len(times),color=COLORS[0],lw=2)
    axs[2].set(title='(c) General-graph test pool',xlabel='MILP runtime (seconds)',ylabel='Empirical cumulative fraction',ylim=(0,1.06))
    axs[2].text(.44,.1,f'{len(g)} instances\n{int((nodes<=1).sum())} at 0-1 nodes',transform=axs[2].transAxes,fontsize=9)
    for ax in axs:ax.grid(alpha=.22)
    fig.tight_layout(w_pad=2);save(fig,'expanded_algorithm_evidence')

def sensitivity():
    rr=rows('37_REAL_CFD_COMMON_PARTITION_PROFILES.csv');summary=rows('42_REAL_CFD_COMMON_PARTITION_SUMMARY.csv')
    keys=[(r['case_id'],r['support']) for r in summary]
    models=['balanced','mesh_dominant','closure_dominant']
    costs=np.zeros((8,3));thresholds=np.geomspace(.1,40,150);front=np.full((8,len(thresholds)),np.nan)
    for i,key in enumerate(keys):
        for j,model in enumerate(models):
            records=[r for r in rr if (r['case_id'],r['support'])==key and r['cost_model']==model]
            costs[i,j]=min(float(r['relaxation_budget']) for r in records if float(r['diameter_relative_gap'])<=1e-7)
            if model=='balanced':
                base=next(float(r['minimum_single_partition_ratio']) for r in records if float(r['relaxation_budget'])==0)
                for k,t in enumerate(thresholds):
                    if base<t:continue
                    hit=[float(r['relaxation_budget']) for r in records if float(r['minimum_single_partition_ratio'])<t]
                    if hit:front[i,k]=min(hit)
    assert np.allclose(costs,np.tile([.5,.2,.5],(8,1)))
    fig,axs=plt.subplots(1,2,figsize=(9,4.4),gridspec_kw={'width_ratios':[1,1.45]})
    ax=axs[0];lab=[f'{chr(71)+str(i//2+1)}-{("fine" if i%2==0 else "shared")}' for i in range(8)]
    for j,model in enumerate(models):ax.plot(np.arange(8)+(j-1)*.12,costs[:,j],'o',color=COLORS[j],label=model.replace('_',' '))
    ax.set(xticks=range(8),xticklabels=lab,ylim=(0,.68),ylabel=r'$\rho_D(10^{-7})$',title='(a) Edge-cost sensitivity');ax.tick_params(axis='x',rotation=40);ax.legend(loc='upper right',fontsize=8)
    ax=axs[1]
    for i,r in enumerate(summary):
        lo=float(r['R_balanced_cost_0p5']);hi=float(r['R_zero_cost'])
        records=[x for x in rr if (x['case_id'],x['support'])==keys[i] and x['cost_model']=='balanced']
        assert all(float(x['minimum_single_partition_ratio'])>=hi-1e-12 for x in records if float(x['relaxation_budget'])<.5)
        ax.hlines(7-i,lo,hi,color=COLORS[i//2],lw=2)
        ax.plot(lo,7-i,'o',mfc='white',mec=COLORS[i//2],ms=6)
        ax.plot(hi,7-i,'o',color=COLORS[i//2],ms=6)
    ax.axvspan(1.6442596248017918,2.2831676147449826,alpha=.12,color='gray')
    ax.set(xscale='log',xlabel=r'Classification threshold $\tau$',yticks=range(8),yticklabels=lab[::-1],ylim=(-.7,7.7),title='(b) Threshold ranges with minimum cost 0.5')
    ax.grid(axis='x',alpha=.2)
    fig.tight_layout(w_pad=3);save(fig,'expanded_cost_sensitivity')
    (HERE/'figure_input_hashes.json').write_text(json.dumps({n:hashlib.sha256((DATA/n).read_bytes()).hexdigest() for n in ['32_FOREST_DP_SCALING_BENCHMARK.csv','30_GRAPH_COUPLING_MILP_BENCHMARK.csv','37_REAL_CFD_COMMON_PARTITION_PROFILES.csv','42_REAL_CFD_COMMON_PARTITION_SUMMARY.csv']},indent=2))

if __name__=='__main__':
    example();performance();sensitivity()
    print('Three additional figures built; example and sensitivity assertions passed.')
