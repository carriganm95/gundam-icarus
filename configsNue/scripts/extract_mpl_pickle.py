"""Dump the plotted data of a pickled spineplot/medulla matplotlib Figure to JSON.

The pickles are written with numpy 2 / recent matplotlib, so they can't be loaded
from the LCG ROOT kernel. Run this with the spineplot pixi python instead:

    python extract_mpl_pickle.py plot_reco_leading_electron_energy.pkl out.json
"""
import json
import pickle
import sys

import numpy as np
from matplotlib.colors import to_hex


def extract_axis(ax):
    out = {
        'xlabel': ax.get_xlabel(), 'ylabel': ax.get_ylabel(),
        'xlim': list(ax.get_xlim()), 'ylim': list(ax.get_ylim()),
        'xscale': ax.get_xscale(), 'yscale': ax.get_yscale(),
        'texts': [{'text': t.get_text(), 'pos': list(t.get_position()),
                   'transform_axes': t.get_transform() == ax.transAxes} for t in ax.texts],
        'stack': [], 'errorbars': [], 'bands': [], 'fills': [], 'steps': [],
    }

    # Step histograms drawn with ax.stairs (e.g. the fractional-uncertainty breakdown plots)
    for p in ax.patches:
        if type(p).__name__ != 'StepPatch':
            continue
        d = p.get_data()
        out['steps'].append({
            'label': p.get_label(),
            'edges': np.asarray(d.edges, dtype=float).tolist(),
            'values': np.asarray(d.values, dtype=float).tolist(),
            'color': to_hex(p.get_edgecolor()),
            'linestyle': p.get_linestyle(),
            'linewidth': p.get_linewidth(),
        })

    # Stacked histogram components (one BarContainer per category)
    for c in ax.containers:
        if type(c).__name__ != 'BarContainer':
            continue
        rects = list(c)
        edges = [r.get_x() for r in rects] + [rects[-1].get_x() + rects[-1].get_width()]
        out['stack'].append({
            'label': rects[0].get_label(),
            'edges': edges,
            'heights': [r.get_height() for r in rects],
            'bottoms': [r.get_y() for r in rects],
            'color': to_hex(rects[0].get_facecolor()),
        })

    # Data points with error bars
    for c in ax.containers:
        if type(c).__name__ != 'ErrorbarContainer':
            continue
        line, _, barcols = c.lines
        segs = barcols[0].get_segments() if barcols else []
        y = np.asarray(line.get_ydata(), dtype=float)
        out['errorbars'].append({
            'label': c.get_label(),
            'x': np.asarray(line.get_xdata(), dtype=float).tolist(),
            'y': y.tolist(),
            'ylo': [float(y[i] - s[0][1]) for i, s in enumerate(segs)],
            'yhi': [float(s[1][1] - y[i]) for i, s in enumerate(segs)],
        })

    # Systematic band drawn as a PatchCollection of rectangles, and fill_between bands
    for col in ax.collections:
        name = type(col).__name__
        if name == 'PatchCollection':
            rects = [p.vertices for p in col.get_paths()]
            out['bands'].append({
                'label': col.get_label(),
                'edges': [float(v[0, 0]) for v in rects] + [float(rects[-1][1, 0])],
                'lo': [float(v[0, 1]) for v in rects],
                'hi': [float(v[2, 1]) for v in rects],
            })
        elif 'FillBetween' in name or name == 'PolyCollection':
            v = col.get_paths()[0].vertices
            out['fills'].append({'label': col.get_label(), 'vertices': v.tolist()})
    return out


if __name__ == '__main__':
    with open(sys.argv[1], 'rb') as f:
        fig = pickle.load(f)
    with open(sys.argv[2], 'w') as f:
        json.dump({'axes': [extract_axis(ax) for ax in fig.axes]}, f, indent=1,
                  default=lambda o: o.item() if hasattr(o, 'item') else o.tolist())
