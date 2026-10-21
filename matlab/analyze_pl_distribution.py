"""Descriptive PL distributions and same-epoch, same-state comparisons."""
from pathlib import Path
import csv
import numpy as np

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'figures/verified_paper/source_data'
OUT = ROOT / 'figures/verified_paper/pl_distribution'
OUT.mkdir(parents=True, exist_ok=True)
METHODS = ['pure_fde1', 'igg3', 'ransac']


def read(site, method):
    with (SOURCE / f'{site}_{method}.csv').open() as f:
        return {r['time_gpst']: {k: float(v) if v else np.nan
                                for k, v in r.items() if k != 'time_gpst'}
                for r in csv.DictReader(f)}


def stats(rows, direction):
    pl = np.array([r[f'{direction}PL_m'] for r in rows])
    pe = np.array([r[f'{direction}PE_m'] for r in rows])
    quantiles = np.quantile(pl, [.1, .25, .5, .75, .9, .95, .99])
    return dict(n=len(rows), **dict(zip(
        ['p10', 'p25', 'median', 'p75', 'p90', 'p95', 'p99'], quantiles)),
        maximum=pl.max(), mean=pl.mean(), coverage_pct=100*np.mean(pe <= pl),
        PE_RMS_m=np.sqrt(np.mean(pe**2)))


distributions, comparisons = [], []
for site in ['JC01', 'JC02', 'JC04']:
    runs = {m: read(site, m) for m in METHODS}
    for state, q in [('fixed', 1), ('float', 2)]:
        for direction in ['H', 'V']:
            valid = {m: {t: r for t, r in run.items()
                         if r['Q'] == q and np.isfinite(r[f'{direction}PL_m'])
                         and r[f'{direction}PL_m'] > 0
                         and np.isfinite(r[f'{direction}PE_m'])}
                     for m, run in runs.items()}
            for method in METHODS:
                result = dict(site=site, state=state, direction=direction,
                              method=method, **stats(list(valid[method].values()), direction))
                distributions.append(result)
                print('OWN', site, state, direction, method,
                      'N', result['n'], 'Q50/90/95/99',
                      '/'.join(f'{result[k]:.5f}' for k in ['median','p90','p95','p99']))
            # Pair each front end with the baseline: both must have the same Q.
            for method in METHODS[1:]:
                times = sorted(valid['pure_fde1'].keys() & valid[method].keys())
                if not times:
                    continue
                a = stats([valid['pure_fde1'][t] for t in times], direction)
                b = stats([valid[method][t] for t in times], direction)
                pa = np.array([valid['pure_fde1'][t][f'{direction}PL_m'] for t in times])
                pb = np.array([valid[method][t][f'{direction}PL_m'] for t in times])
                row = dict(site=site, state=state, direction=direction, method=method,
                           n=len(times), smaller_pct=100*np.mean(pb < pa),
                           median_epoch_ratio=np.median(pb/pa))
                for k in a:
                    if k != 'n':
                        row[f'baseline_{k}'] = a[k]
                        row[f'robust_{k}'] = b[k]
                comparisons.append(row)
                print('PAIRED', site, state, direction, method, 'N', len(times),
                      'Q50/95 baseline', f"{a['median']:.5f}/{a['p95']:.5f}",
                      'robust', f"{b['median']:.5f}/{b['p95']:.5f}",
                      'smaller%', f"{row['smaller_pct']:.2f}",
                      'coverage', f"{a['coverage_pct']:.2f}/{b['coverage_pct']:.2f}")

for name, rows in [('distributions.csv', distributions), ('paired_comparisons.csv', comparisons)]:
    with (OUT / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
