"""Nonvisual preparation for MATLAB; use every verified-run input epoch."""
from pathlib import Path
import csv,json,math,re
from statistics import NormalDist
import numpy as np

REPO=Path(__file__).resolve().parent.parent
ROOT=REPO/'tmp/gls_fde_verified_20261015'
OUT=REPO/'matlab/figures/verified_paper/source_data'
OUT.mkdir(parents=True,exist_ok=True)
def keyed(path):
    return {' '.join(p[:2]):p for l in path.read_text(errors='replace').splitlines() if l.strip() and not l.startswith('%') for p in [l.split()]}
summary=[];audit=[]
for site in ['JC01','JC02','JC04']:
    text=(REPO/'matlab/result'/site/'plot_enu_error.m').read_text(encoding='utf-8-sig')
    ref=np.array(list(map(float,re.search(r'truthEcef\s*=\s*\[([^]]+)\]',text)[1].split(','))))
    x,y,z=ref;lon=math.atan2(y,x);rho=math.hypot(x,y);e2=6.6943799901413165e-3
    lat=math.atan2(z,rho*(1-e2))
    for _ in range(12):
        n=6378137/math.sqrt(1-e2*math.sin(lat)**2);lat=math.atan2(z+e2*n*math.sin(lat),rho)
    sl,cl,sp,cp=math.sin(lon),math.cos(lon),math.sin(lat),math.cos(lat)
    E=np.array([[-sl,cl,0],[-sp*cl,-sp*sl,cp],[cp*cl,cp*sl,sp]])
    for method in ['pure_fde1','igg3','ransac']:
        folder=ROOT/'verified'/method
        pos=keyed(folder/f'{site}.pos');diag=keyed(folder/f'{site}.glsdiag')
        epochs=list(csv.DictReader((ROOT/f'{site}_{method}_verified_epochs.csv').open()))
        assert len(epochs)==2864
        cols=['time_gpst','Q','HPE_m','VPE_m','HPL_m','VPL_m','Ksigma_H_m','Ksigma_V_m','sigma_H_m','sigma_V_m','Kmd']
        vals=[];neg=0
        for row in epochs:
            a={k:row[k] for k in cols if k in row}
            t=row['time_gpst'].replace('-','/')+'.000'
            if t in pos:
                p=pos[t];std=np.array(list(map(float,p[7:13])))
                Q=np.diag(std[:3]**2)
                Q[0,1]=Q[1,0]=math.copysign(std[3]**2,std[3])
                Q[1,2]=Q[2,1]=math.copysign(std[4]**2,std[4])
                Q[0,2]=Q[2,0]=math.copysign(std[5]**2,std[5])
                qe=E@Q@E.T
                vv=np.array([qe[0,0]+qe[1,1],qe[2,2]])
                neg+=int(np.any(vv<0))
                sig=np.sqrt(np.maximum(vv,0))
                k=-NormalDist().inv_cdf(1e-3/(6*(int(diag[t][3])+1)))
                a.update(Ksigma_H_m=k*sig[0],Ksigma_V_m=k*sig[1],sigma_H_m=sig[0],sigma_V_m=sig[1],Kmd=k)
            vals.append(a)
        with (OUT/f'{site}_{method}.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(vals)
        array=np.array([[float(r[k]) if r.get(k,'')!='' else np.nan for k in cols[1:]] for r in vals])
        for group,q in [('all',None),('fixed',1),('float',2)]:
            output=np.isfinite(array[:,0]) if q is None else array[:,0]==q
            for direction,pe,pl,ks in [('H',1,3,5),('V',2,4,6)]:
                paired=output&np.isfinite(array[:,[pe,pl,ks]]).all(axis=1)
                n=int(paired.sum());noutput=int(output.sum())
                record=dict(site=site,method=method,group=group,direction=direction,input=2864,
                    output=noutput,paired=n,availability_pct=100*np.sum(output&np.isfinite(array[:,3:5]).all(axis=1))/2864,
                    fixed_pct=100*np.sum(array[:,0]==1)/np.sum(np.isfinite(array[:,0])),
                    PL_coverage_pct=100*np.mean(array[paired,pe]<=array[paired,pl]) if n else None,
                    Ksigma_coverage_pct=100*np.mean(array[paired,pe]<=array[paired,ks]) if n else None,
                    PL_median_m=float(np.median(array[paired,pl])) if n else None,
                    Ksigma_median_m=float(np.median(array[paired,ks])) if n else None,
                    PE_RMS_m=float(np.sqrt(np.mean(array[output,pe]**2))) if noutput else None,
                    zero_PE=int(np.sum(output&(array[:,pe]==0))),zero_Ksigma=int(np.sum(output&(array[:,ks]==0))))
                summary.append(record)
        audit.append(dict(site=site,method=method,input=2864,output=len(pos),negative_HV_variance_after_rounded_covariance_rotation=neg))
with (OUT/'metrics.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
(OUT/'audit.json').write_text(json.dumps(audit,indent=2))
print(f'Prepared {len(audit)} runs, {len(summary)} group/direction summaries; covariance is already posterior-scaled.')
for r in summary:
    if r['group']=='all':print(r['site'],r['method'],r['direction'],round(r['PL_coverage_pct'],2),round(r['Ksigma_coverage_pct'],2))
