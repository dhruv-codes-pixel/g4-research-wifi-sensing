"""Loader for the Brescia 20 MHz Nexmon-CSI dataset (Tonini & Lo Cigno, Zenodo). MAT v7.3 files with one variable `capture`
(256 FFT rows x N frames, compound real/imag float64). Timestamps are NOT in the MAT file; they are read from the matching pcap
(libpcap record headers; the pcap packet count equals the MAT frame count for every capture)."""
import glob, json, os, struct
import h5py, numpy as np
import csi_methods as M

ROOT_DEFAULT = "data/20MHz"

def load_mat(path):
    """Return complex (N, 256) array (frames x FFT rows)."""
    with h5py.File(path, "r") as f:
        a = f["capture"][()]
    return (a["real"] + 1j * a["imag"]).T

def load_active(path):
    """Complex (N, 242) CSI on the 242 usable tones (k = HE20_KIDX)."""
    return load_mat(path)[:, M.HE20_ROWS]

def pcap_timestamps(path):
    """Seconds (float) of every packet record of a classic libpcap file."""
    with open(path, "rb") as f:
        g = f.read(24); magic = struct.unpack("<I", g[:4])[0]
        end = "<" if magic in (0xa1b2c3d4, 0xa1b23c4d) else ">"; scale = 1e-9 if magic == 0xa1b23c4d else 1e-6
        ts = []
        while True:
            h = f.read(16)
            if len(h) < 16: break
            s, u, il, ol = struct.unpack(end + "IIII", h); f.seek(il, 1); ts.append(s + u * scale)
    return np.array(ts)

def list_captures(root=ROOT_DEFAULT):
    """All captures with their config.json: list of dicts (class label, path, config...)."""
    out = []
    for mat in sorted(glob.glob(os.path.join(root, "*", "*", "capture*.mat"))):
        d = os.path.dirname(mat); cfg = json.load(open(os.path.join(d, "config.json")))
        people = cfg["people"]["num"]; date = os.path.basename(d)
        cls = {"0people": "E_day" if date == "2024-07-12" else "E_night", "1people": "S1", "4people": "M4"}[os.path.basename(os.path.dirname(d))]
        out.append(dict(cls=cls, people=people, date=date, capture=os.path.basename(mat)[:-4], mat=mat, pcap=mat[:-4] + ".pcap", code=cfg["code"], notes=cfg["notes"]))
    return out
