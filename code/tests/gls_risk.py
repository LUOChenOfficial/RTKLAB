"""Offline research prototype: conditional mode mixture, never production PL.

Inputs must partition the modeled cases and provide explicit probability masses.
Formal Gaussian covariances and observed separations are NOT calibrated
conditional-on-selection distributions. The caller must validate that model.
No fault probability is inferred from a RANSAC residual or consensus count.
"""
import argparse
import json
import math
from pathlib import Path


def tail_bound(limit, separation, sigma):
    """Triangle/union bound: P(|noise| > L-|separation|), capped at one."""
    if not all(math.isfinite(x) for x in (limit,separation,sigma)) or sigma<0:
        raise ValueError('Non-finite component or negative sigma')
    margin=limit-abs(separation)
    if margin<=0:return 1.0
    return math.erfc(margin/(math.sqrt(2)*sigma)) if sigma else 0.0


def solve_axis(components, budget, unmodeled_mass):
    """components: (probability, observed separation, conditional sigma).

    Probabilities partition modeled cases; unmodeled mass uses worst-case tail 1.
    This is an assumed-distribution risk calculation, not an integrity guarantee.
    """
    if not 0<budget<1 or not 0<=unmodeled_mass<=1:
        raise ValueError('Invalid risk budget or unmodeled mass')
    if not components or any(not math.isfinite(p) or p<0 for p,_,_ in components):
        raise ValueError('Missing components or invalid probability')
    if abs(sum(p for p,_,_ in components)+unmodeled_mass-1)>1e-10:
        raise ValueError('Mode probabilities plus unmodeled mass must sum to one')
    for _,sep,sigma in components:tail_bound(0,sep,sigma)
    if unmodeled_mass>=budget:return None
    def risk(limit):
        return unmodeled_mass+sum(p*tail_bound(limit,sep,sigma) for p,sep,sigma in components)
    if risk(0)<=budget:return 0.0
    lo,hi=0.0,1.0
    while risk(hi)>budget:
        hi*=2
        if hi>1e9:return None
    for _ in range(90):
        mid=(lo+hi)/2
        if risk(mid)>budget:lo=mid
        else:hi=mid
    return hi


def load_modes(path):
    epochs={}
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('%'):continue
        p=line.split();time=' '.join(p[:2]);key=':'.join(p[2:5])
        rows=epochs.setdefault(time,{})
        if key in rows:raise ValueError('Duplicate mode '+time+' '+key)
        rows[key]={'sep':list(map(float,p[5:8])),'sigma':list(map(float,p[8:11]))}
    return epochs


def run(mode_path, diag_path, model_path):
    modes=load_modes(mode_path)
    model=json.loads(Path(model_path).read_text(encoding='utf-8'))
    # Each epoch requires its own complete partition; no missing-mass default.
    output=[]
    for line in Path(diag_path).read_text().splitlines():
        if not line or line.startswith('%'):continue
        d=line.split();time=' '.join(d[:2]);rows=modes.get(time,{})
        result={'time':time,'status':'unavailable','hpl_m':None,'vpl_m':None}
        output.append(result)
        if int(d[5]) or int(d[6]) or not math.isfinite(float(d[10])):
            result['reason']='baseline missing modes, alert, or model failure';continue
        if len(rows)!=int(d[3])+1:
            result['reason']='incomplete mode export';continue
        weights=model['epochs'].get(time)
        if weights is None:
            result['reason']='no explicit residual-risk model for epoch';continue
        if set(weights['probabilities'])!=set(rows):
            raise ValueError('Probability keys must match every exported mode: '+time)
        components=[[(weights['probabilities'][key],r['sep'][j],r['sigma'][j]) for key,r in rows.items()] for j in range(3)]
        limits=[solve_axis(c,model['horizontal_budget']/2 if j<2 else model['vertical_budget'],weights['unmodeled_mass']) for j,c in enumerate(components)]
        if any(x is None for x in limits):
            result['reason']='unmodeled risk consumes budget';continue
        result.update(status='experimental_assumed_model',hpl_m=math.hypot(*limits[:2]),vpl_m=limits[2],axis_limits_m=limits)
    return {'model_label':model['label'],'integrity_validated':False,'epochs':output}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('modes');parser.add_argument('diagnostics');parser.add_argument('risk_model');parser.add_argument('output')
    args=parser.parse_args()
    Path(args.output).write_text(json.dumps(run(args.modes,args.diagnostics,args.risk_model),indent=2,allow_nan=False),encoding='utf-8')
