"""Vectorized version of the offline assumed-model risk solver."""
import math
import numpy as np
from scipy.special import erfc


def solve_batch(sep,sigma,weight,budget=5e-8):
    if sep.shape!=sigma.shape or weight.shape!=sep.shape[:2]:
        raise ValueError('Invalid component dimensions')
    if not (np.isfinite(sep).all() and np.isfinite(sigma).all() and np.isfinite(weight).all()):
        raise ValueError('Nonfinite input')
    if np.any(sigma<0) or np.any(weight<0) or not np.allclose(weight.sum(axis=1),1,rtol=0,atol=1e-12):
        raise ValueError('Invalid sigma or probability partition')
    if not 0<budget<1:raise ValueError('Invalid budget')
    budgets=np.array([budget/2,budget/2,budget])
    def risk(limit):
        margin=limit[:,None,:]-np.abs(sep)
        tail=erfc(np.maximum(margin,0)/(math.sqrt(2)*np.maximum(sigma,1e-300)))
        tail=np.where(sigma==0,0.0,tail)
        tail=np.where(margin<=0,1.0,tail)
        return np.sum(weight[:,:,None]*tail,axis=1)
    lo=np.zeros((sep.shape[0],3));hi=np.ones_like(lo)
    for _ in range(40):
        mask=risk(hi)>budgets
        if not mask.any():break
        hi=np.where(mask,hi*2,hi)
    else:raise ValueError('Unable to bracket risk limit')
    for _ in range(45):
        mid=(lo+hi)/2;mask=risk(mid)>budgets
        lo=np.where(mask,mid,lo);hi=np.where(mask,hi,mid)
    return np.hypot(hi[:,0],hi[:,1]),hi[:,2]
