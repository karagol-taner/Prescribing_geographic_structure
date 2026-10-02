import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Supplementary Figures S3 to S7 (phylogeographic signal, concentration, strips, space-time, national trends)."""
import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],'font.size':10,
  'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':11,'axes.titleweight':'bold',
  'figure.dpi':200,'savefig.dpi':200,'axes.labelsize':10})
CT='#2c3e50'; GENUINE='#1a7a5a'; POINT='#c0392b'; GREY='#95a5a6'; BLUE='#2f6db3'; ACC='#e08a1e'
OUT='../figures'; os.makedirs(OUT, exist_ok=True)
f1=lambda s: float(s.iloc[0])
R=pd.read_csv('mantel_af.csv')                 # per-item Mantel results: Drug,Mantel_r,P_Value,items,top_share
man84=pd.read_csv('mantel_84.csv')                # mantel_84.py
T=pd.read_csv('temporal_af.csv')
A=np.load('A_af.npy')
drugs=list(R['Drug']); idx={c.strip().lower():i for i,c in enumerate(drugs)}
def col(name):
    i=idx.get(name.strip().lower()); return A[:,i] if i is not None else None
def topshare(name):
    v=col(name)
    if v is None: return np.nan
    s=v.sum(); return (v.max()/s) if s>0 else np.nan
dk=['catheter','bag','ostomy','stoma','appliance','sheath','irrigation','faecal','drainage','adhesive','belt','tubing','protector','filler','plug','valve','filter','shield','plate','deodorant','lubricant','discharge solidifying','sensor','suspensory','collar','pouch','wipe','glove','dressing','bottle','teat','cup','syringe','needle','notification','urinal']
R['dev']=R['Drug'].str.lower().apply(lambda x:any(k in x for k in dk) and not any(e in x for e in ['baguette','eye tear']))
sigmask=R.P_Value<0.05
nstrong=int((sigmask&(R.Mantel_r>0.6)).sum()); nmod=int((sigmask&(R.Mantel_r>=0.3)&(R.Mantel_r<=0.6)).sum()); nweak=int((sigmask&(R.Mantel_r<0.3)).sum()); nsig=int(sigmask.sum())

# ---- Figure S3: phylogeographic signal ----
fig,ax=plt.subplots(1,2,figsize=(9,3.6))
sig=R[sigmask]['Mantel_r']; ns=R[~sigmask]['Mantel_r']
ax[0].hist(ns,bins=40,color=GREY,alpha=.8,label=f'P≥0.05 (n={len(ns):,})',edgecolor='white',linewidth=.3)
ax[0].hist(sig,bins=40,color=GENUINE,alpha=.85,label=f'P<0.05 (n={len(sig):,})',edgecolor='white',linewidth=.3)
ax[0].axvline(0,color=CT,ls='--',lw=1); ax[0].set_xlabel('Phylogeographic signal (Mantel r)'); ax[0].set_ylabel('Number of items')
ax[0].set_title(f'A  Phylogeographic signal, {len(R):,} items',loc='left'); ax[0].legend(frameon=False,fontsize=8.5)
lab=['Weak\n(r<0.3)','Moderate\n(0.3-0.6)','Strong\n(r>0.6)']; val=[nweak,nmod,nstrong]; colr=[GREY,ACC,POINT]
b=ax[1].bar(lab,val,color=colr,edgecolor='white'); ax[1].set_ylabel('Number of significant items')
ax[1].set_title(f'B  Strength of the {nsig} signals with P<0.05',loc='left'); ax[1].set_ylim(0,max(val)*1.22)
for r,v in zip(b,val): ax[1].text(r.get_x()+r.get_width()/2,v+2,f'{v}\n({100*v/nsig:.0f}%)',ha='center',fontsize=8.5)
plt.tight_layout(); plt.savefig(f'{OUT}/figS3_signal_distribution.png',bbox_inches='tight',facecolor='white'); plt.close()

# ---- Figure S6: 84 profiles versus 42 areas; dispersion by year ----
m=man84.rename(columns={'Mantel_r':'r84','Drug_Name':'Drug'})[['Drug','r84']].copy()
m['k']=m.Drug.str.strip().str.lower(); R['k']=R.Drug.str.strip().str.lower()
m=m.merge(R.rename(columns={'Mantel_r':'r42'})[['k','r42','P_Value','dev']],on='k',how='inner')
modern={'Dapagliflozin':'Dapagliflozin','Tirzepatide':'Tirzepatide','Liraglutide':'Liraglutide','Semaglutide':'Semaglutide'}
est={'Chlortalidone':'Chlortalidone','Ketoprofen':'Ketoprofen','Hydrocortisone butyrate':'Hydrocortisone butyrate','Lisinopril':'Lisinopril'}
fig,ax=plt.subplots(1,2,figsize=(11,5))
# panel A
ax[0].scatter(m['r84'],m['r42'],s=8,color=GREY,alpha=.35,edgecolor='none')
lim=[-0.25,0.75]; ax[0].plot(lim,lim,color=CT,ls='--',lw=1,label='no change')
for d in modern:
    row=m[m.Drug.str.lower()==d.lower()]
    if len(row): ax[0].scatter(row['r84'],row['r42'],s=55,color=POINT,zorder=5,edgecolor='white',lw=.6)
for d in est:
    row=m[m.Drug.str.lower()==d.lower()]
    if len(row): ax[0].scatter(row['r84'],row['r42'],s=55,color=BLUE,zorder=5,edgecolor='white',lw=.6)
for d,lb in {**modern,**est}.items():
    row=m[m.Drug.str.lower()==d.lower()]
    if len(row): ax[0].annotate(lb,(f1(row['r84']),f1(row['r42'])),textcoords='offset points',xytext=(7,-2),fontsize=8,color=(POINT if d in modern else BLUE),
                                path_effects=[pe.withStroke(linewidth=2.5,foreground='white')])
ax[0].scatter([],[],color=POINT,label='Fast-changing modern agents')
ax[0].scatter([],[],color=BLUE,label='Established interchangeable agents')
ax[0].set_xlabel('Mantel r, 84 area-period profiles (place and time mixed)'); ax[0].set_ylabel('Mantel r, 42 areas (time held constant)')
ax[0].set_title('A  Merging commissioning eras removes temporal signal',loc='left'); ax[0].set_xlim(lim); ax[0].set_ylim(lim); ax[0].legend(frameon=False,fontsize=7.5,loc='upper left')
# panel B: between-area dispersion (CV) trajectory over years, relative to 2021
import pickle
shares=np.load('shares_af.npy'); yrs=np.load('years_af.npy')
areas_td, drugs_td = pickle.load(open('td_index.pkl','rb'))
tdidx={d.strip().lower():i for i,d in enumerate(drugs_td)}
def cv_traj(name):
    i=tdidx.get(name.strip().lower())
    if i is None: return None
    out=[]
    for k in range(len(yrs)):
        v=shares[k,:,i]; mu=v.mean(); out.append(v.std()/mu if mu>0 else np.nan)
    return np.array(out)
for name,c in [('Dapagliflozin',POINT),('Semaglutide',ACC),('Lisinopril',BLUE),('Ketoprofen',GENUINE)]:
    tr=cv_traj(name)
    if tr is None or not np.isfinite(tr[0]) or tr[0]==0: continue
    ax[1].plot(yrs, tr/tr[0], '-o', color=c, ms=4, lw=1.9, label=name)
ax[1].axhline(1.0,color=CT,ls=':',lw=.9)
ax[1].set_xlabel('Year'); ax[1].set_ylabel('Between-area dispersion (CV, relative to 2021)')
ax[1].set_title('B  Between-area dispersion falls as new agents spread',loc='left')
ax[1].set_xticks(yrs); ax[1].set_ylim(0.55,1.15); ax[1].legend(frameon=False,fontsize=8,loc='lower left')
plt.tight_layout(); plt.savefig(f'{OUT}/figS6_space_time.png',bbox_inches='tight',facecolor='white'); plt.close()

# ---- Figure S4: concentration ----
sg=R[sigmask].copy(); sg['top']=sg['Drug'].apply(topshare)
fig,ax=plt.subplots(figsize=(6.6,5))
d1=sg[~sg.dev]; d2=sg[sg.dev]
ax.scatter(100*d2['top'],d2['Mantel_r'],s=26,color=GREY,alpha=.7,label='Dispensing appliance (supply artefact)',edgecolor='none')
ax.scatter(100*d1['top'],d1['Mantel_r'],s=30,color=GENUINE,alpha=.8,label='Medicine / other',edgecolor='none')
ax.axvline(40,color=CT,ls=':',lw=1); ax.text(41.5,0.535,'single-area\nconcentration >40%\n(point source)',fontsize=7.5,color=CT,va='top')
OFF={'Chlortalidone':(-5,3,'right'),'Hydrocortisone butyrate':(5,-9,'left'),'Ketoprofen':(-5,-9,'right'),'Padimate O':(-5,3,'right')}   # clear of other points
for d in ['Thyrotropin alfa','Cladribine (Immunomodulating)','Padimate O','Sucrase','Chlortalidone','Ketoprofen','Hydrocortisone butyrate','Lisinopril','Dapagliflozin']:
    row=sg[sg.Drug.str.lower()==d.lower()]
    if len(row):
        x=100*f1(row['top']); y=f1(row['Mantel_r']); nm=d.replace(' (Immunomodulating)','')
        dx,dy,ha=OFF.get(nm,(5,3,'left'))
        ax.annotate(nm,(x,y),textcoords='offset points',xytext=(dx,dy),ha=ha,fontsize=7.5,color=(POINT if x>40 else BLUE),
                    path_effects=[pe.withStroke(linewidth=2.5,foreground='white')])
ax.set_xlabel('Volume-adjusted share of use in the single highest area (%)'); ax.set_ylabel('Phylogeographic signal (Mantel r, 42 areas)')
ax.set_title('Volume-adjusted concentration separates point-source\nartefacts from distributed patterns')
ax.legend(frameon=False,fontsize=8,loc='lower right'); plt.tight_layout()
plt.savefig(f'{OUT}/figS4_concentration.png',bbox_inches='tight',facecolor='white'); plt.close()

# ---- Figure S5: three regimes ----
fig,ax=plt.subplots(1,3,figsize=(9.2,3.4))
for k,(d,ttl,c) in enumerate([('Thyrotropin alfa','Thyrotropin alfa\n(point source)',POINT),('Chlortalidone','Chlortalidone\n(distributed, long tail)',GENUINE),('Dapagliflozin','Dapagliflozin\n(narrow range)',BLUE)]):
    v=col(d); v=np.sort(v/v.sum())[::-1]*100
    ax[k].bar(range(len(v)),v,color=c,edgecolor='white',linewidth=.2)
    ax[k].set_title(ttl,fontsize=10); ax[k].set_xlabel('42 areas (ranked)')
    if k==0: ax[k].set_ylabel('Volume-adjusted share of use (%)')
    ax[k].axhline(100/42,color=CT,ls='--',lw=.8)
plt.tight_layout(); plt.savefig(f'{OUT}/figS5_strips.png',bbox_inches='tight',facecolor='white'); plt.close()

# ---- Figure S7: national trends (from temporal_data.csv) ----
Tt=pd.read_csv('temporal_data.csv')
gg=Tt.groupby('Period').agg(items=('Total_Items','sum'),cost=('Total_Cost','sum'),chem=('Unique_Chemicals','mean'),cv=('Total_Items',lambda s:s.std()/s.mean())).reset_index()
per=gg['Period'].values; x=range(len(per))
fig,ax=plt.subplots(2,2,figsize=(9,6))
ax[0,0].plot(x,gg['items']/1e9,'-o',color=BLUE,ms=3); ax[0,0].set_title('A  Dispensed items (billions)',loc='left')
ax[0,1].plot(x,gg['cost']/1e9,'-o',color=GENUINE,ms=3); ax[0,1].set_title('B  Net ingredient cost (GBP billions)',loc='left')
ax[1,0].plot(x,gg['chem'],'-o',color=ACC,ms=3); ax[1,0].set_title('C  Mean unique chemicals per area',loc='left')
ax[1,0].annotate('coding\ndiscontinuity\n(2023)',(8,gg['chem'].iloc[8]),xytext=(2.5,1150),fontsize=7.5,arrowprops=dict(arrowstyle='->',color=CT))
ax[1,1].plot(x,gg['cv'],'-o',color=POINT,ms=3); ax[1,1].set_title('D  Between-area CV of volume',loc='left')
for a in ax.flat:
    a.set_xticks(range(0,len(per),max(1,len(per)//5))); a.set_xticklabels([per[i] for i in range(0,len(per),max(1,len(per)//5))],rotation=45,ha='right',fontsize=7.5)
plt.tight_layout(); plt.savefig(f'{OUT}/figS7_temporal.png',bbox_inches='tight',facecolor='white'); plt.close()
print(f'strong={nstrong} mod={nmod} weak={nweak} sig={nsig}')
