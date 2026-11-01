"""
digitize_fig6.py -- semi-automatic digitisation of Ratnam et al. (arXiv v2) Fig. 6(a) and 6(b).

Fig. 6 is embedded as two raster images on page 9 of the arXiv PDF. Extract them first:
    pdfimages -png -f 9 -l 9 <paper.pdf> /tmp/fig6        ->  /tmp/fig6-000.png (a), /tmp/fig6-001.png (b)
then:  python digitize_fig6.py /tmp/fig6-000.png /tmp/fig6-001.png ../Results/gain/fig6_paper_digitised.csv

Calibration: x from the gamma gridlines (0.6 ... 0.9) and the axis box; y (log scale) by least squares on the
pixel rows of the labelled major/minor gridlines.  Markers: colour mask per series, the marker is the largest
row-cluster of that colour inside a 13-px column window centred on each gamma of the paper's grid.
An overlay PNG with every detected point is written next to the CSV for visual verification.
Accuracy ~ +-3 % in value (1-2 px); overlapping markers near gamma > 0.94 are less reliable.
This replaces manual WebPlotDigitizer clicking; spot-check a few points in WebPlotDigitizer if required.
"""
import sys, csv, os
import numpy as np
from PIL import Image, ImageDraw

GAM = [0.5, 0.6, 0.7, 0.8, 0.85, 0.87, 0.9, 0.92, 0.94, 0.96, 0.97, 0.98, 0.99]


def masks(a):
    a = a.astype(int); R, G, B = a[..., 0], a[..., 1], a[..., 2]
    return {"Algorithm 3": (R < 90) & (G < 90) & (B < 90), "Algorithm 1": (B > 180) & (R < 120) & (G < 120),
            "Baseline (6)": (R > 180) & (G < 110) & (B < 110), "Baseline (8)": (G > 170) & (R < 170) & (B < 120) & (G - R > 40)}


CFG = {"a": dict(box=(82, 506), c=lambda g: 82 + (g - 0.5) * 850, rows=[54, 196, 339], vals=[1000, 100, 10],
                 legend=(85, 300, 250, 382), gmax=0.99),
       "b": dict(box=(73, 506), c=lambda g: 72 + (g - 0.5) * 870, rows=[39, 79, 135, 231, 246, 262, 281, 302, 328, 359],
                 vals=[40, 30, 20, 10, 9, 8, 7, 6, 5, 4], legend=(380, 40, 506, 135), gmax=0.98)}


def digitise(path, C):
    a = np.array(Image.open(path).convert("RGB")); M = masks(a)
    A = np.stack([np.ones(len(C["vals"])), -np.log10(C["vals"])], 1)
    r0, s = np.linalg.lstsq(A, np.array(C["rows"], float), rcond=None)[0]
    yv = lambda r: 10 ** ((r0 - r) / s)
    lx0, ly0, lx1, ly1 = C["legend"]; x0, x1 = C["box"]
    img = Image.open(path).convert("RGB"); dr = ImageDraw.Draw(img); out = []
    for name, m in M.items():
        mm = np.zeros_like(m); mm[35:382, x0 + 3:x1 - 1] = m[35:382, x0 + 3:x1 - 1]; mm[ly0:ly1, lx0:lx1] = False
        for g in GAM:
            if g > C["gmax"]: continue
            col = int(round(C["c"](g)))
            if name == "Baseline (6)":                       # solid line without markers
                rows = np.where(mm[:, max(col - 1, 0):col + 2].any(1))[0]
                if len(rows) == 0: continue
                r = float(np.median(rows))
            else:
                w = mm[:, col - 6:col + 7]; rows = np.where(w.any(1))[0]
                if len(rows) == 0: continue
                groups = np.split(rows, np.where(np.diff(rows) > 3)[0] + 1)
                gr = groups[int(np.argmax([w[x].sum() for x in groups]))]
                r = float(np.average(gr, weights=w[gr].sum(1)))
            out.append((name, g, float(yv(r))))
            dr.ellipse([col - 3, r - 3, col + 3, r + 3], outline=(255, 0, 255), width=2)
    return out, img


if __name__ == "__main__":
    pa, pb, dst = sys.argv[1:4]
    rows = []
    for tag, p in (("a", pa), ("b", pb)):
        pts, img = digitise(p, CFG[tag])
        img.resize((img.width * 2, img.height * 2)).save(os.path.join(os.path.dirname(dst), f"fig6{tag}_digitisation_overlay.png"))
        rows += [dict(panel=f"6({tag})", dynamic_type="i" if tag == "a" else "ii", series=n, gamma=g,
                      snr_metric_chi2_over_1_minus_chi2=round(v, 3)) for n, g, v in pts]
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"wrote {len(rows)} points to {dst}")
