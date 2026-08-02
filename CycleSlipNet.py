#!/usr/bin/env python3
"""
CycleSlipNet
============
Three-model deep learning GNSS carrier-phase cycle-slip detection.
Models  : Transformer | CNN | LSTM
"""

import os, sys, re, argparse, warnings, itertools, time, gc
import numpy as np
import pandas as pd
import joblib
from collections import defaultdict
from scipy.stats import chi2
from sklearn.metrics import confusion_matrix

if sys.platform == "win32":
    import io
    if hasattr(sys.stdout, "buffer"):
        sys.stdout = io.TextIOWrapper(
            sys.stdout.buffer, encoding="utf-8", errors="replace", line_buffering=True)
    if hasattr(sys.stderr, "buffer"):
        sys.stderr = io.TextIOWrapper(
            sys.stderr.buffer, encoding="utf-8", errors="replace", line_buffering=True)

warnings.filterwarnings("ignore")


def _safe_print(*args, **kwargs):
    try:
        print(*args, **kwargs)
    except (OSError, UnicodeEncodeError):
        text = " ".join(str(a) for a in args)
        safe = text.encode(sys.stdout.encoding or "utf-8", errors="replace").decode(
            sys.stdout.encoding or "utf-8", errors="replace")
        try:
            sys.stdout.write(safe + "\n"); sys.stdout.flush()
        except Exception:
            pass


import tensorflow as tf
try:
    from tensorflow.keras.models import load_model
    import tensorflow.keras.saving as keras_saving
except ImportError:
    raise ImportError("pip install tensorflow")


# ==============================================================
# §0  Global configuration
# ==============================================================
DEFAULT_SLIP_PATH  = r"D:\data10.0\onsaSDEDclean\slip_normalized.out"
DEFAULT_RINEX_PATH = r"D:\data10.0\15s\onsa0450.24o"
DEFAULT_TF_DIR     = r"E:\train\out\Transformer"
DEFAULT_CNN_DIR    = r"E:\train\out\CNN"
DEFAULT_LSTM_DIR   = r"E:\train\out\LSTM"
DEFAULT_OUTPUT_DIR = r"E:\CycleSlipNet_output"
DEFAULT_SEED       = 42
DEFAULT_THRESHOLD  = 2.0
DEFAULT_D          = 1
DEFAULT_SCENARIO   = "all"

MONTE_CARLO_N = 15
BOOTSTRAP_B   = 1000

METRIC_BOUNDS = {
    "Precision": (0.0, 1.0),
    "Recall":    (0.0, 1.0),
    "F1":        (0.0, 1.0),
    "Accuracy":  (0.0, 1.0),
    "RSR":       (0.0, 1.0),
    "MRE":       (0.0, None),
}

DL_CONF_THRESHOLD = 1   # detection threshold in wavelength multiples

CARRIER_CANDS = {
    0: ["L2I", "L1C", "L1P", "L1X"],
    1: ["L7I", "L5Q", "L5X", "L5P"],
    2: ["L6I", "L8X", "L6X", "L3X"],
}
# Reference wavelengths (m); actual values are read per-row from slip.out
WAVELENGTH = {0: 0.192, 1: 0.248, 2: 0.255}

PAL = {"TF": "#2878b5", "CNN": "#c82423", "LSTM": "#d62728"}

DL_MODELS = {
    "TF":   {"full_name": "Transformer", "color": PAL["TF"],   "ls": "-"},
    "CNN":  {"full_name": "CNN",          "color": PAL["CNN"],  "ls": "--"},
    "LSTM": {"full_name": "LSTM",         "color": PAL["LSTM"], "ls": ":"},
}

SCENARIOS = {
    "sc_A": {"desc": "Scenario A: small slips (0-10%, 1-10 cyc)",
             "ratio_min": 0.00, "ratio_max": 0.10, "slip_min": 1,  "slip_max": 10},
    "sc_B": {"desc": "Scenario B: engineering baseline (0-10%, 1-50 cyc)",
             "ratio_min": 0.00, "ratio_max": 0.10, "slip_min": 1,  "slip_max": 50},
    "sc_C": {"desc": "Scenario C: moderate concurrency (10-20%, 1-50 cyc)",
             "ratio_min": 0.10, "ratio_max": 0.20, "slip_min": 1,  "slip_max": 50},
    "sc_D": {"desc": "Scenario D: extreme stress (20-30%, 1-100 cyc)",
             "ratio_min": 0.20, "ratio_max": 0.30, "slip_min": 1,  "slip_max": 100},
}


def _slip_bins(sc_def):
    mx = sc_def["slip_max"]
    edges_all  = [0.5, 1.5, 2.5, 5.5, 10.5, 20.5, 50.5, 9999.5]
    labels_all = ["1", "2", "3-5", "6-10", "11-20", "21-50", ">50"]
    keep_idx = None
    for i, right in enumerate(edges_all[1:]):
        if right >= mx:
            keep_idx = i; break
    if keep_idx is None:
        keep_idx = len(labels_all) - 1
    return edges_all[:keep_idx + 2], labels_all[:keep_idx + 1]


# ==============================================================
# §1  Transformer positional encoding
# ==============================================================
@keras_saving.register_keras_serializable()
class PositionalEncoding(tf.keras.layers.Layer):
    def __init__(self, sequence_length, d_model, **kwargs):
        super().__init__(**kwargs)
        self.sequence_length = sequence_length
        self.d_model         = d_model
        pos = np.arange(sequence_length)[:, np.newaxis]
        dt  = np.exp(np.arange(0, d_model, 2) * -(np.log(10000.) / d_model))
        pe  = np.zeros((sequence_length, d_model))
        pe[:, 0::2] = np.sin(pos * dt)
        pe[:, 1::2] = np.cos(pos * dt)[:, :d_model // 2]
        self.pos_encoding = tf.convert_to_tensor(pe, dtype=tf.float32)

    def call(self, inputs):
        return inputs + self.pos_encoding

    def get_config(self):
        c = super().get_config()
        c.update({"sequence_length": self.sequence_length, "d_model": self.d_model})
        return c

    @classmethod
    def from_config(cls, config):
        return cls(**config)


# ==============================================================
# §3  Parsing
# ==============================================================
_SLIP_RE = re.compile(
    r"Epoch:\s*(\d+),\s*PRN=(\d+),\s*freq_index=(\d+),"
    r"\s*wavelength=([\d.eE+\-]+),\s*TDCP cycle_slip=([\d.eE+\-]+),"
    r"\s*Simulated slip=([\d.eE+\-]+)(.*)"
)


def parse_slip_out(path, dt_sec=1.0):
    _safe_print(f"\n[Parse] slip.out: {path}  (dt={dt_sec:.3f}s)")
    recs = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = _SLIP_RE.search(line.rstrip("\n"))
            if not m:
                continue
            tdcp_raw = float(m.group(5))
            recs.append({
                "Epoch":      int(m.group(1)),
                "PRN":        int(m.group(2)),
                "freq_index": int(m.group(3)),
                "wavelength": float(m.group(4)),
                "TDCP_m":     tdcp_raw,
                "TDCP_norm":  tdcp_raw / dt_sec,
                "sim_slip":   float(m.group(6)),
                "tail":       m.group(7),
            })
    if not recs:
        raise ValueError(f"No records parsed: {path}")
    df = pd.DataFrame(recs)
    df = df.sort_values(["PRN", "freq_index", "Epoch"]).reset_index(drop=True)
    prns = sorted(df["PRN"].unique())
    _safe_print(f"  {len(df)} records  {df['Epoch'].nunique()} epochs  "
                f"{len(prns)} satellites  PRN={prns}")
    for fi in sorted(df["freq_index"].unique()):
        _safe_print(f"  freq_index={fi}: {len(df[df['freq_index'] == fi])} rows")
    return df


def _parse_rinex_time(line):
    parts = line.split()
    if len(parts) < 7:
        return None
    try:
        return int(parts[4]) * 3600. + int(parts[5]) * 60. + float(parts[6])
    except Exception:
        return None


def parse_rinex3_obs(path):
    _safe_print(f"\n[Parse] RINEX: {path}")
    hdr, sys_obs, elist, odict = [], {}, [], {}
    in_hdr = True; eidx = -1; dt_sec = None
    with open(path, "r", encoding="ascii", errors="ignore") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        raw = lines[i].rstrip("\n")
        if in_hdr:
            hdr.append(raw)
            if "SYS / # / OBS TYPES" in raw:
                sc = raw[0].strip()
                if sc:
                    n = int(raw[3:6]); ts = raw[7:60].split()
                    sys_obs.setdefault(sc, []).extend(ts)
                    while len(sys_obs[sc]) < n:
                        i += 1; c = lines[i].rstrip("\n"); hdr.append(c)
                        sys_obs[sc].extend(c[7:60].split())
            if "INTERVAL" in raw:
                try:
                    dt_sec = float(raw[:10].strip())
                except Exception:
                    pass
            if "END OF HEADER" in raw:
                in_hdr = False
            i += 1; continue
        if raw.startswith(">"):
            parts = raw.split()
            flag  = int(parts[7]) if len(parts) > 7 else 0
            eidx += 1; elist.append((eidx, raw, flag, 0)); odict[eidx] = {}
            i += 1; continue
        if eidx >= 0 and len(raw) >= 3:
            sid = raw[0:3].strip(); sc = sid[0] if sid else ""
            if sc in sys_obs:
                ts = sys_obs[sc]; obs_str = raw[3:]; vals = {}
                for k, ot in enumerate(ts):
                    seg = obs_str[k * 16:k * 16 + 14].strip()
                    try:
                        v = float(seg); vals[ot] = v if v != 0. else None
                    except Exception:
                        vals[ot] = None
                odict[eidx][sid] = vals
        i += 1
    if dt_sec is None or dt_sec <= 0:
        et = [t for _, el, _, _ in elist[:200] if (t := _parse_rinex_time(el))]
        if len(et) >= 2:
            diffs = [abs(et[j+1]-et[j]) for j in range(len(et)-1)
                     if 0 < abs(et[j+1]-et[j]) <= 3600]
            if diffs:
                dt_sec = float(np.median(diffs))
        if not dt_sec or dt_sec <= 0:
            dt_sec = 1.0
    _safe_print(f"  {len(elist)} epochs  systems={list(sys_obs.keys())}  dt={dt_sec:.3f}s")
    return odict, elist, sys_obs, hdr, dt_sec


def resolve_obs(sys_obs):
    all_t = set(t for ts in sys_obs.values() for t in ts)
    cm = {}
    for fi in range(3):
        for c in CARRIER_CANDS[fi]:
            if c in all_t:
                cm[fi] = c; break
    _safe_print("  Carrier obs codes: " + " | ".join(
        f"F{fi}:{cm.get(fi,'---')}" for fi in range(3)))
    for fi in range(3):
        if fi not in cm:
            _safe_print(f"  WARNING: freq_index={fi} has no matching obs type in RINEX "
                        f"(candidates: {CARRIER_CANDS[fi]}). "
                        f"Injection still covers this freq via slip.out.")
    return cm


def _prn2sat(prn):
    if 300 < prn < 400: return f"C{prn-300:02d}"
    if 100 < prn < 200: return f"G{prn-100:02d}"
    return f"C{prn:02d}"


# ==============================================================
# §4  Cycle-slip injection
# ==============================================================
def _epoch_inject_single(df_raw, odict, elist, cm, sc,
                         master_seed, mc_run, dt_sec=1.0):
    """
    Inject simulated cycle slips into df_raw.

    epoch_pool candidates are drawn directly from df_raw — not filtered
    through rinex_avail — so every (PRN, freq_index, Epoch) triple that
    exists in slip.out is eligible for injection regardless of whether
    the corresponding obs code appears in the RINEX file.
    """
    slip_epoch_min = int(df_raw["Epoch"].min())

    def ep2idx(ep):
        return ep - slip_epoch_min

    df_sim = df_raw.copy()
    if "sim_slip" not in df_sim.columns:
        df_sim["sim_slip"] = 0.0

    # Map (PRN, freq_index, Epoch) -> DataFrame row index
    slip_index = {
        (int(r.PRN), int(r.freq_index), int(r.Epoch)): idx
        for idx, r in df_sim.iterrows()
    }

    # RINEX phase observations (kept in sync for downstream use)
    oi = {ei: {s: dict(v) for s, v in sv.items()} for ei, sv in odict.items()}

    # Build epoch_pool directly from df_raw (no RINEX filter)
    epoch_pool = defaultdict(list)
    for (prn, fi, ep_abs) in slip_index.keys():
        epoch_pool[ep_abs].append((prn, fi))

    fi_count = defaultdict(int)
    for ep_abs, cands in epoch_pool.items():
        for (prn, fi) in cands:
            fi_count[fi] += 1
    _safe_print(f"  [MC={mc_run:02d}] epoch_pool freq coverage: "
                + ", ".join(f"F{fi}:{cnt}" for fi, cnt in sorted(fi_count.items())))

    if not epoch_pool:
        _safe_print("  [Inject] Warning: epoch_pool is empty")
        return df_sim, oi, {}, pd.DataFrame()

    slip_dict = {}; inject_log_rows = []

    for ep_abs in sorted(epoch_pool.keys()):
        candidates = epoch_pool[ep_abs]
        N_obs_e    = len(candidates)
        if N_obs_e == 0:
            continue
        ep_seed  = (master_seed * 100003 + ep_abs * 7 + mc_run * 31 + 20240101)
        rng      = np.random.default_rng(ep_seed & 0xFFFFFFFFFFFFFFFF)
        r_e      = float(rng.uniform(sc["ratio_min"], sc["ratio_max"]))
        N_slip_e = max(0, min(int(round(r_e * N_obs_e)), N_obs_e))
        inject_log_rows.append({
            "epoch":         ep_abs,
            "N_obs_e":       N_obs_e,
            "sampled_ratio": round(r_e, 6),
            "N_slip_e":      N_slip_e,
            "actual_ratio":  round(N_slip_e / N_obs_e, 6) if N_obs_e > 0 else 0.,
        })
        if N_slip_e == 0:
            continue

        idx_perm       = rng.permutation(len(candidates))
        selected       = [candidates[idx_perm[k]] for k in range(N_slip_e)]
        slip_cands_arr = np.arange(sc["slip_min"], sc["slip_max"] + 1)
        weights        = np.ones(len(slip_cands_arr), dtype=float)
        if sc["slip_min"] == 1:
            weights[0] = 0.1          # downweight single-cycle slips
        weights /= weights.sum()

        for (prn, fi) in selected:
            slip_abs = int(rng.choice(slip_cands_arr, p=weights))
            slip_val = slip_abs * (1 if rng.random() > 0.5 else -1)
            key      = (prn, fi, ep_abs)
            row_idx  = slip_index.get(key)
            if row_idx is None:
                continue
            wl = float(df_sim.at[row_idx, "wavelength"])
            df_sim.at[row_idx, "TDCP_norm"] += slip_val * wl / dt_sec
            df_sim.at[row_idx, "TDCP_m"]    += slip_val * wl
            df_sim.at[row_idx, "sim_slip"]  += slip_val
            slip_dict[key] = slip_val

            sat_id = _prn2sat(prn)
            ei     = ep2idx(ep_abs)
            ct     = cm.get(fi)
            if ct and 0 <= ei < len(elist):
                sv = oi.get(ei, {}).get(sat_id)
                if sv and sv.get(ct) is not None:
                    sv[ct] += slip_val

    inject_log = pd.DataFrame(inject_log_rows)
    sizes = [abs(v) for v in slip_dict.values()]
    if sizes:
        _safe_print(f"  [MC={mc_run:02d}] Injected {len(slip_dict)} slips  "
                    f"magnitude [{min(sizes)}-{max(sizes)}] cyc")
    else:
        _safe_print(f"  [MC={mc_run:02d}] Injected 0 slips")

    for fi in sorted(df_sim["freq_index"].unique()):
        n_slip = int((df_sim[df_sim["freq_index"] == fi]["sim_slip"] != 0).sum())
        _safe_print(f"  [MC={mc_run:02d}] freq_index={fi}: {n_slip} slips injected")

    return df_sim, oi, slip_dict, inject_log


def save_inject_log(il, sc_name, mc_run, out_dir):
    if il.empty: return
    il.to_csv(os.path.join(out_dir, f"inject_log_{sc_name}_mc{mc_run:02d}.csv"), index=False)


def summarize_inject_logs(ll, sc_name, out_dir):
    if not ll: return
    al = pd.concat(ll, ignore_index=True)
    if al.empty or "actual_ratio" not in al.columns: return
    pd.DataFrame([{
        "scenario":   sc_name,
        "mc_runs":    len(ll),
        "mean_ratio": round(float(al["actual_ratio"].mean()), 6),
    }]).to_csv(os.path.join(out_dir, f"inject_summary_{sc_name}.csv"), index=False)


# ==============================================================
# §5  DL model loading and inference
# ==============================================================
def load_dl_model(model_dir, model_short_name="TF"):
    candidates = [
        "transformer_model.keras", "cnn_model.keras", "lstm_model.keras",
        f"{model_short_name.lower()}_model.keras", "model.keras", "best_model.keras",
    ]
    model_path = None
    for fname in candidates:
        p = os.path.join(model_dir, fname)
        if os.path.exists(p):
            model_path = p; break
    if model_path is None:
        for fname in os.listdir(model_dir):
            if fname.endswith(".keras"):
                model_path = os.path.join(model_dir, fname); break
    if model_path is None:
        raise FileNotFoundError(f"[{model_short_name}] No .keras found in: {model_dir}")
    m  = load_model(model_path, custom_objects={"PositionalEncoding": PositionalEncoding})
    ls = joblib.load(os.path.join(model_dir, "label_scaler.pkl"))
    ws = joblib.load(os.path.join(model_dir, "window_size.pkl"))
    try:
        d_model = m.input_shape[-1]
    except Exception:
        d_model = 1
    _safe_print(f"  [{model_short_name}] {os.path.basename(model_path)} "
                f"window={ws} d_model={d_model}")
    return m, ls, ws


def _best_integer_estimate(sf_raw, max_search=2):
    base = round(sf_raw); frac = abs(sf_raw - base)
    if frac <= 0.35:
        return float(base)
    best_int = base; best_dist = frac
    for c in range(int(base) - max_search, int(base) + max_search + 1):
        d = abs(sf_raw - c)
        if d < best_dist:
            best_dist = d; best_int = c
    return float(best_int)


def _build_window(ts_norm, ws):
    """1-D feature window: raw TDCP_norm sequence, shape (1, ws, 1)."""
    tw = np.array(ts_norm[-ws:], dtype=float).flatten()
    if len(tw) < ws:
        tw = np.pad(tw, (ws - len(tw), 0), mode="edge")
    return tw.reshape(1, ws, 1)


def _flush_dl(m, bf, bm, ls, res, mname, dt_sec):
    X = np.vstack(bf)
    try:
        raw_out = m(X.astype(np.float32), training=False).numpy().flatten()
    except Exception:
        raw_out = m.predict(X, batch_size=len(X), verbose=0).flatten()

    pred_norm_arr = ls.inverse_transform(raw_out.reshape(-1, 1)).flatten()

    det_col   = f"Slip_Detected_{mname}"
    conf_col  = f"DL_Confidence_{mname}"
    size_col  = f"Slip_Size_{mname}_cycles"
    rep_col   = f"TDCP_repaired_{mname}_m"
    float_col = f"DL_SlipFloat_{mname}"

    for i, pred_norm in enumerate(pred_norm_arr):
        mt       = bm[i]
        row      = mt["row"]
        obs_norm = float(row["obs_norm"])
        wl       = float(row["wavelength"])

        delta_norm     = pred_norm - obs_norm
        thr_norm       = DL_CONF_THRESHOLD * wl / dt_sec
        slip_float_raw = -delta_norm * dt_sec / max(wl, 1e-9)
        conf           = abs(delta_norm * dt_sec) / max(wl, 1e-9)

        if abs(delta_norm) >= thr_norm:
            slip_int = _best_integer_estimate(slip_float_raw)
            if abs(slip_int) < 1.0:
                det = False; slip_int = 0.; repaired_norm = obs_norm
            else:
                det = True
                repaired_norm = obs_norm + slip_int * wl / dt_sec
        else:
            det = False; slip_int = 0.; repaired_norm = obs_norm; conf = 0.

        res["Epoch"].append(int(row["Epoch"]))
        res["PRN"].append(int(row["PRN"]))
        res["freq_index"].append(int(row["freq_index"]))
        res["wavelength"].append(wl)
        res["TDCP_obs_m"].append(float(row["TDCP_m"]))
        res["TDCP_norm_obs"].append(obs_norm)
        res["TDCP_norm_pred"].append(pred_norm)
        res["Delta_norm"].append(delta_norm)
        res["TDCP_pred_m"].append(pred_norm * dt_sec)
        res["Delta_TDCP_m"].append(delta_norm * dt_sec)
        res["Threshold_norm"].append(thr_norm)
        res["Threshold_m"].append(thr_norm * dt_sec)
        res[det_col].append(bool(det))
        res[conf_col].append(round(conf, 4))
        res[size_col].append(slip_int if det else 0.)
        res[float_col].append(round(slip_float_raw, 6))
        res[rep_col].append(repaired_norm * dt_sec)
        res["sim_slip"].append(float(row["sim_slip"]))


def run_dl_model(m, ls, ws, df_sim, mname, dt_sec=1., fi_list=None, batch=1024):
    t0 = time.time()
    df = df_sim.copy()
    if fi_list:
        df = df[df["freq_index"].isin(fi_list)].copy()

    det_col = f"Slip_Detected_{mname}"
    res = defaultdict(list)
    gw: dict = defaultdict(list)
    bf = []; bm = []

    for gk, grp in df.groupby(["PRN", "freq_index"]):
        grp = grp.sort_values("Epoch").reset_index(drop=True)
        gw[gk] = []
        prev_epoch = None

        for nt in grp.itertuples(index=False):
            obs_norm  = float(nt.TDCP_norm)
            obs_m     = float(nt.TDCP_m)
            cur_epoch = int(nt.Epoch)
            wl        = float(nt.wavelength)

            if prev_epoch is not None and (cur_epoch - prev_epoch) != 1:
                gw[gk] = []
            prev_epoch = cur_epoch
            gw[gk].append(obs_norm)
            if len(gw[gk]) > ws + 10:
                gw[gk].pop(0)

            if len(gw[gk]) < ws:
                res["Epoch"].append(cur_epoch)
                res["PRN"].append(int(nt.PRN))
                res["freq_index"].append(int(nt.freq_index))
                res["wavelength"].append(wl)
                res["TDCP_obs_m"].append(obs_m)
                res["TDCP_norm_obs"].append(obs_norm)
                res["TDCP_norm_pred"].append(float("nan"))
                res["Delta_norm"].append(float("nan"))
                res["TDCP_pred_m"].append(float("nan"))
                res["Delta_TDCP_m"].append(float("nan"))
                res["Threshold_norm"].append(float("nan"))
                res["Threshold_m"].append(float("nan"))
                res[det_col].append(False)
                res[f"DL_Confidence_{mname}"].append(0.)
                res[f"Slip_Size_{mname}_cycles"].append(0.)
                res[f"DL_SlipFloat_{mname}"].append(0.)
                res[f"TDCP_repaired_{mname}_m"].append(obs_m)
                res["sim_slip"].append(float(nt.sim_slip))
                continue

            bf.append(_build_window(gw[gk], ws))
            bm.append({"row": {
                "Epoch":      cur_epoch,
                "PRN":        int(nt.PRN),
                "freq_index": int(nt.freq_index),
                "wavelength": wl,
                "TDCP_m":     obs_m,
                "obs_norm":   obs_norm,
                "sim_slip":   float(nt.sim_slip),
            }})

            if len(bf) >= batch:
                _flush_dl(m, bf, bm, ls, res, mname, dt_sec)
                bf.clear(); bm.clear()

    if bf:
        _flush_dl(m, bf, bm, ls, res, mname, dt_sec)
        bf.clear(); bm.clear()

    gc.collect()
    dr = (pd.DataFrame(dict(res))
          .sort_values(["PRN", "freq_index", "Epoch"])
          .reset_index(drop=True))
    elapsed = time.time() - t0
    dr.attrs["elapsed"]   = elapsed
    dr.attrs["n_samples"] = len(dr)
    _safe_print(f"  [{mname}] freq={fi_list}  "
                f"Detected:{int(dr[det_col].sum())}  {elapsed:.2f}s")
    return dr


# ==============================================================
# §7  Metrics and McNemar test
# ==============================================================
def _repair_metrics(true_s, det_s, est_s):
    ta = np.asarray(true_s); da = np.asarray(det_s).astype(bool); ea = np.asarray(est_s)
    mask = (ta != 0) & da; n = int(mask.sum())
    if n == 0:
        return {"RSR": 0., "MRE": 0., "n_repaired": 0}
    err = np.abs(ea[mask] - ta[mask])
    return {"RSR": round(float((err == 0).mean()), 4),
            "MRE": round(float(err.mean()), 4),
            "n_repaired": n}


def calc_metrics(true_s, det_s, est_s, tag=""):
    tl     = (np.asarray(true_s) != 0).astype(int)
    pl     = np.asarray(det_s).astype(int)
    cm_mat = confusion_matrix(tl, pl, labels=[0, 1])
    tn, fp, fn, tp = cm_mat.ravel(); tot = tn + fp + fn + tp
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.
    rec  = tp / (tp + fn) if (tp + fn) > 0 else 0.
    f1   = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.
    acc  = (tp + tn) / tot if tot > 0 else 0.
    rm   = _repair_metrics(true_s, det_s, est_s)
    _safe_print(f"  [{tag}]")
    _safe_print(f"    P={prec:.4f}  R={rec:.4f}  F1={f1:.4f}  Acc={acc:.4f}")
    _safe_print(f"    RSR={rm['RSR']:.4f}  MRE={rm['MRE']:.4f}  repaired={rm['n_repaired']}")
    _safe_print(f"    TP={tp}  FP={fp}  FN={fn}  TN={tn}  total={tot}  "
                f"true_slips={int(tl.sum())}")
    return {
        "tag": tag,
        "TN": int(tn), "FP": int(fp), "FN": int(fn), "TP": int(tp),
        "Accuracy":   round(acc,  4),
        "Precision":  round(prec, 4),
        "Recall":     round(rec,  4),
        "F1":         round(f1,   4),
        "RSR":        rm["RSR"],
        "MRE":        rm["MRE"],
        "n_repaired": rm["n_repaired"],
        "total":      tot,
        "true_slips": int(tl.sum()),
        "_cm":        cm_mat,
    }


def mcnemar_test(y_true, pred1, pred2, label1="M1", label2="M2"):
    tl = (np.asarray(y_true).ravel() != 0).astype(bool)
    p1 = np.asarray(pred1).ravel().astype(bool)
    p2 = np.asarray(pred2).ravel().astype(bool)
    b  = int(((p1 == tl) & (p2 != tl)).sum())
    c  = int(((p1 != tl) & (p2 == tl)).sum()); n = b + c
    if n == 0:
        stat, p_val = 0., 1.
    else:
        stat  = max(0., (abs(b - c) - 1.) ** 2 / n)
        p_val = float(chi2.sf(stat, df=1))
    sig    = ("***" if p_val < 0.001 else "**" if p_val < 0.01
              else "*" if p_val < 0.05 else "ns")
    better = (f"{label1}>>{label2}" if b > c
              else f"{label2}>>{label1}" if c > b else "equal")
    _safe_print(f"    McNemar [{label1} vs {label2}]: "
                f"b={b} c={c} chi2={stat:.4f} p={p_val:.6f} {sig} [{better}]")
    return {
        "label1": label1, "label2": label2, "b": b, "c": c,
        "chi2": round(stat, 4), "p_value": round(p_val, 6),
        "significance": sig, "better": better,
    }


# ==============================================================
# §7b  Error record export
# ==============================================================
def save_error_records(df, det_col, size_col, method_name, out_dir, tag):
    if df.empty or "sim_slip" not in df.columns: return
    os.makedirs(out_dir, exist_ok=True)
    det_arr  = df[det_col].astype(bool)
    true_arr = df["sim_slip"] != 0
    fp_df    = df[det_arr & ~true_arr].copy();  fp_df["error_type"] = "FP"
    fn_df    = df[~det_arr & true_arr].copy();  fn_df["error_type"] = "FN"
    misdet   = pd.concat([fp_df, fn_df], ignore_index=True)
    if not misdet.empty:
        misdet.to_csv(
            os.path.join(out_dir, f"misdetect_{method_name}_{tag}.csv"), index=False)
        _safe_print(f"    [{method_name}] Misdetect: FP={len(fp_df)} FN={len(fn_df)}")
    tp_df = df[det_arr & true_arr].copy()
    wrong = tp_df[tp_df[size_col] != tp_df["sim_slip"]].copy()
    if not wrong.empty:
        wrong["repair_error"] = wrong[size_col] - wrong["sim_slip"]
        wrong["error_type"]   = "wrong_repair"
        wrong.to_csv(
            os.path.join(out_dir, f"wrong_repair_{method_name}_{tag}.csv"), index=False)
        _safe_print(f"    [{method_name}] Wrong repairs: {len(wrong)}")


# ==============================================================
# §7c  Metric printing and CSV export
# ==============================================================
def print_and_save_metrics(df, det_col, size_col, method_name,
                           out_dir, tag, filter_notnan_col=None):
    dv = df.copy()
    if filter_notnan_col and filter_notnan_col in dv.columns:
        dv = dv[~dv[filter_notnan_col].isna()]
    if len(dv) == 0:
        _safe_print(f"  [{method_name}] No valid samples, skipped"); return None
    n_true = int((dv["sim_slip"] != 0).sum())
    if n_true == 0:
        _safe_print(f"  [{method_name}] Warning: true_slips=0  total={len(dv)}")
        os.makedirs(out_dir, exist_ok=True)
        pd.DataFrame([{
            "tag": f"{method_name}_{tag}", "warning": "no_true_slips",
            "n_total": len(dv),
            "n_detected": int(dv[det_col].astype(bool).sum()),
        }]).to_csv(
            os.path.join(out_dir, f"metrics_{method_name}_{tag}_diag.csv"), index=False)
        return None
    mt      = calc_metrics(dv["sim_slip"], dv[det_col], dv[size_col],
                           f"{method_name}_{tag}")
    mt_copy = {k: v for k, v in mt.items() if k != "_cm"}
    os.makedirs(out_dir, exist_ok=True)
    pd.DataFrame([mt_copy]).to_csv(
        os.path.join(out_dir, f"metrics_{method_name}_{tag}.csv"), index=False)
    _safe_print(f"    -> metrics_{method_name}_{tag}.csv")
    save_error_records(dv, det_col, size_col, method_name, out_dir, tag)
    return mt


# ==============================================================
# §7d  Bootstrap 95% CI
# ==============================================================
def _compute_bootstrap_ci95(values, metric_name="", B=None):
    if B is None:
        B = BOOTSTRAP_B
    arr = np.array(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    n   = len(arr)
    if n == 0:
        return 0., 0., 0., 0.
    mean = float(np.mean(arr))
    std  = float(np.std(arr, ddof=1)) if n > 1 else 0.
    if n == 1:
        return mean, std, mean, mean
    rng        = np.random.default_rng(seed=12345)
    boot_means = np.empty(B, dtype=float)
    for i in range(B):
        sample        = rng.choice(arr, size=n, replace=True)
        boot_means[i] = sample.mean()
    ci_lo = float(np.percentile(boot_means, 2.5))
    ci_hi = float(np.percentile(boot_means, 97.5))
    bounds    = METRIC_BOUNDS.get(metric_name, (None, None))
    lo_bound, hi_bound = bounds
    if lo_bound is not None:
        ci_lo = max(ci_lo, lo_bound); ci_hi = max(ci_hi, lo_bound)
    if hi_bound is not None:
        ci_hi = min(ci_hi, hi_bound); ci_lo = min(ci_lo, hi_bound)
    return mean, std, ci_lo, ci_hi


# ==============================================================
# §7e  MC summary with Bootstrap CI
# ==============================================================
def summarize_mc_metrics_with_ci(mc_metrics_list, out_dir, tag_prefix=""):
    if not mc_metrics_list:
        return pd.DataFrame()
    df          = pd.DataFrame(mc_metrics_list)
    metric_cols = ["Precision", "Recall", "F1", "Accuracy", "RSR", "MRE"]
    rows        = []
    group_cols  = []
    if "scenario" in df.columns: group_cols.append("scenario")
    if "model"    in df.columns: group_cols.append("model")
    if "freq"     in df.columns: group_cols.append("freq")
    if not group_cols:            group_cols = ["tag"]

    for gkeys, grp in df.groupby(group_cols, sort=True):
        if not isinstance(gkeys, tuple):
            gkeys = (gkeys,)
        base = dict(zip(group_cols, gkeys))
        base["n_runs"] = len(grp)
        for mc in metric_cols:
            if mc not in grp.columns: continue
            vals = grp[mc].dropna().tolist()
            mean, std, ci_lo, ci_hi = _compute_bootstrap_ci95(vals, metric_name=mc)
            rows.append({
                **base,
                "metric":      mc,
                "mean":        round(mean,  4),
                "std":         round(std,   4),
                "ci95_lo":     round(ci_lo, 4),
                "ci95_hi":     round(ci_hi, 4),
                "ci_method":   "bootstrap_percentile",
                "bootstrap_B": BOOTSTRAP_B,
            })

    result = pd.DataFrame(rows)
    if not result.empty:
        fname = (f"metrics_with_CI_{tag_prefix}.csv" if tag_prefix
                 else "metrics_with_CI.csv")
        result.to_csv(os.path.join(out_dir, fname), index=False)
        _safe_print(f"    -> {fname}  ({len(result)} rows, "
                    f"Bootstrap {BOOTSTRAP_B}-resample 95% CI, "
                    f"grouped by {group_cols})")
    return result


# ==============================================================
# §8  Slip-size stratified statistics
# ==============================================================
def compute_slip_size_stats(df, det_col, size_col, method_name,
                            sc_name, ft, sc_def, filter_notnan_col=None):
    dv = df.copy()
    if filter_notnan_col and filter_notnan_col in dv.columns:
        dv = dv[~dv[filter_notnan_col].isna()]
    dv_slip = dv[dv["sim_slip"] != 0].copy()
    if len(dv_slip) == 0:
        return pd.DataFrame()
    dv_slip["abs_slip"] = dv_slip["sim_slip"].abs().astype(int)
    edges, labels = _slip_bins(sc_def)
    dv_slip["size_grp"] = pd.cut(
        dv_slip["abs_slip"], bins=edges, labels=labels,
        right=True, include_lowest=True)
    rows = []
    for gl in labels:
        grp = dv_slip[dv_slip["size_grp"] == gl]; total = len(grp)
        if total == 0:
            rows.append({
                "scenario": sc_name, "method": method_name, "freq_tag": ft,
                "slip_size_group": str(gl), "total_simulated": 0,
                "detected": 0, "missed": 0, "detection_rate": 0.,
                "repair_correct": 0, "repair_wrong": 0, "RSR": 0.}); continue
        det     = int(grp[det_col].astype(bool).sum())
        dm_mask = grp[det_col].astype(bool)
        if det > 0 and size_col in grp.columns:
            est  = grp.loc[dm_mask, size_col].values
            true = grp.loc[dm_mask, "sim_slip"].values
            rok  = int((est == true).sum()); rer = det - rok; rsr = rok / det
        else:
            rok = rer = 0; rsr = 0.
        rows.append({
            "scenario": sc_name, "method": method_name, "freq_tag": ft,
            "slip_size_group": str(gl), "total_simulated": total,
            "detected": det, "missed": total - det,
            "detection_rate": round(det / total, 4),
            "repair_correct": rok, "repair_wrong": rer,
            "RSR": round(rsr, 4)})
    ts = sum(r["total_simulated"] for r in rows)
    ds = sum(r["detected"]        for r in rows)
    ro = sum(r["repair_correct"]  for r in rows)
    rows.append({
        "scenario": sc_name, "method": method_name, "freq_tag": ft,
        "slip_size_group": "TOTAL", "total_simulated": ts,
        "detected": ds, "missed": ts - ds,
        "detection_rate": round(ds / ts, 4) if ts > 0 else 0.,
        "repair_correct": ro, "repair_wrong": ds - ro,
        "RSR": round(ro / ds, 4) if ds > 0 else 0.})
    return pd.DataFrame(rows)


def enum_combos(freqs):
    combos = []
    for r in range(1, len(freqs) + 1):
        for c in itertools.combinations(freqs, r):
            combos.append(list(c))
    return combos


def ftag(fi_list):
    return "F" + "F".join(str(f) for f in sorted(fi_list))


def mkdir(path):
    os.makedirs(path, exist_ok=True); return path


# ==============================================================
# §18  True slip count statistics
# ==============================================================
def report_true_slip_counts(df_sim, slip_dict, sc_name, out_dir):
    mkdir(out_dir)
    n_total    = len(df_sim)
    n_slips    = int((df_sim["sim_slip"] != 0).sum())
    slip_sizes = [abs(v) for v in slip_dict.values()] if slip_dict else []
    report = {
        "scenario":       sc_name,
        "total_obs":      n_total,
        "n_cycle_slips":  n_slips,
        "n_injected":     len(slip_dict) if slip_dict else 0,
        "slip_ratio_pct": round(100. * n_slips / max(n_total, 1), 4),
        "min_slip_size":  int(min(slip_sizes))             if slip_sizes else 0,
        "max_slip_size":  int(max(slip_sizes))             if slip_sizes else 0,
        "mean_slip_size": round(float(np.mean(slip_sizes)), 2) if slip_sizes else 0.,
    }
    pd.DataFrame([report]).to_csv(
        os.path.join(out_dir, f"true_slip_counts_{sc_name}.csv"), index=False)
    _safe_print(f"  [Slip count] {sc_name}: {n_slips}/{n_total} "
                f"({report['slip_ratio_pct']:.2f}%)")
    return report


# ==============================================================
# §12  Single-scenario run
# ==============================================================
def run_scenario(sc_name, sc, df_raw, odict, elist, cm_map,
                 dl_model_dict, seed, dt_sec, out_dir):
    sc_root = mkdir(os.path.join(out_dir, sc_name))
    _safe_print(f"\n{'='*66}")
    _safe_print(f"  Scenario: {sc_name}  {sc['desc']}")
    _safe_print(f"  Bootstrap CI: B={BOOTSTRAP_B} resamples, seed=12345")
    _safe_print(f"{'='*66}")

    all_freqs = sorted(df_raw["freq_index"].unique().tolist())
    combos    = enum_combos(all_freqs)
    _safe_print(f"  Freq indices in data: {all_freqs}  Combos: {len(combos)}")

    mc_all_metrics = []; mc_inject_logs = []
    inject_log_dir = mkdir(os.path.join(sc_root, "inject_logs"))
    all_slip_stats = []; all_mcnemar = []

    for mc in range(MONTE_CARLO_N):
        _safe_print(f"\n  -- MC={mc:02d} --")
        df_sim, oi, slip_dict, inject_log = _epoch_inject_single(
            df_raw, odict if odict else {}, elist if elist else [],
            cm_map, sc, seed, mc, dt_sec)
        mc_inject_logs.append(inject_log)
        save_inject_log(inject_log, sc_name, mc, inject_log_dir)

        report_true_slip_counts(
            df_sim, slip_dict, f"{sc_name}_mc{mc:02d}",
            mkdir(os.path.join(sc_root, "slip_counts")))

        combo_model_preds = defaultdict(dict)

        for mname, (m, ls, ws) in dl_model_dict.items():
            _safe_print(f"\n  {'─'*60}")
            _safe_print(f"  Model: {mname}  MC={mc:02d}")
            _safe_print(f"  {'─'*60}")

            for fi_list in combos:
                ft  = ftag(fi_list)
                tag = f"{sc_name}_{ft}_mc{mc:02d}"

                _safe_print(f"\n  [DL]  {mname}  {ft}")
                df_dl = run_dl_model(m, ls, ws, df_sim, mname, dt_sec, fi_list)

                dv      = df_dl[~df_dl["TDCP_pred_m"].isna()].copy()
                det_col = f"Slip_Detected_{mname}"
                siz_col = f"Slip_Size_{mname}_cycles"
                dl_dir  = mkdir(os.path.join(sc_root, "dl_only"))

                mt = print_and_save_metrics(
                    dv, det_col, siz_col, mname, dl_dir, tag,
                    filter_notnan_col="TDCP_pred_m")
                if mt:
                    mt.pop("_cm", None)
                    mt["mc_run"]   = mc
                    mt["model"]    = mname
                    mt["freq"]     = ft
                    mt["scenario"] = sc_name
                    mc_all_metrics.append(mt)

                    ss = compute_slip_size_stats(
                        dv, det_col, siz_col, mname, sc_name, ft, sc,
                        filter_notnan_col="TDCP_pred_m")
                    if not ss.empty:
                        all_slip_stats.append(ss)

                if len(dv) > 0 and dv["sim_slip"].abs().sum() > 0:
                    combo_model_preds[ft][mname] = (
                        dv["sim_slip"].values, dv[det_col].values)

        for ft, model_preds in combo_model_preds.items():
            models = sorted(model_preds.keys())
            for idx_a, ma in enumerate(models):
                for idx_b, mb in enumerate(models):
                    if idx_b <= idx_a: continue
                    true_s, pred_a = model_preds[ma]
                    _,      pred_b = model_preds[mb]
                    if len(true_s) != len(pred_a): continue
                    try:
                        mn = mcnemar_test(true_s, pred_a, pred_b, ma, mb)
                        mn.update({"scenario": sc_name, "freq": ft, "mc": mc})
                        all_mcnemar.append(mn)
                    except Exception:
                        pass

    summary_dir = mkdir(os.path.join(sc_root, "summary"))
    summarize_inject_logs(mc_inject_logs, sc_name, summary_dir)

    if mc_all_metrics:
        mc_df = pd.DataFrame(mc_all_metrics)
        mc_df.to_csv(os.path.join(summary_dir, "ALL_models_metrics.csv"), index=False)
        summarize_mc_metrics_with_ci(mc_all_metrics, summary_dir, tag_prefix=sc_name)

        _safe_print(f"\n  {'='*62}")
        _safe_print(f"  [{sc_name}] Cross-model summary (mean +/- Bootstrap 95% CI)")
        _safe_print(f"  {'─'*62}")
        _safe_print(f"  {'Model':<14} {'P':>8} {'R':>8} {'F1':>8} {'RSR':>8} {'MRE':>8}")
        for mname in sorted(mc_df["model"].unique()
                            if "model" in mc_df.columns else []):
            sub = mc_df[mc_df["model"] == mname]
            p_m, _, p_lo, p_hi = _compute_bootstrap_ci95(
                sub["Precision"].tolist(), "Precision")
            r_m, _, r_lo, r_hi = _compute_bootstrap_ci95(
                sub["Recall"].tolist(), "Recall")
            f_m, _, f_lo, f_hi = _compute_bootstrap_ci95(
                sub["F1"].tolist(), "F1")
            rs_m = sub["RSR"].mean(); mr_m = sub["MRE"].mean()
            _safe_print(f"  {DL_MODELS[mname]['full_name']:<14} "
                        f"{p_m:>7.4f} {r_m:>8.4f} {f_m:>8.4f} "
                        f"{rs_m:>8.4f} {mr_m:>8.4f}")
            _safe_print(f"  {'':14} "
                        f"[{p_lo:.4f},{p_hi:.4f}] "
                        f"[{r_lo:.4f},{r_hi:.4f}] "
                        f"[{f_lo:.4f},{f_hi:.4f}]")
        _safe_print(f"  {'='*62}")

    if all_mcnemar:
        pd.DataFrame(all_mcnemar).to_csv(
            os.path.join(summary_dir, "mcnemar_cross_model.csv"), index=False)
    if all_slip_stats:
        pd.concat(all_slip_stats, ignore_index=True).to_csv(
            os.path.join(summary_dir, "slip_size_stats_all.csv"), index=False)

    _safe_print(f"\n  Scenario {sc_name} done -> {sc_root}")
    return mc_all_metrics


# ==============================================================
# §13  Main pipeline
# ==============================================================
def run(slip_path, rinex_path, tf_dir, cnn_dir, out_dir,
        seed=DEFAULT_SEED, fd=DEFAULT_D, scenario=DEFAULT_SCENARIO,
        lstm_dir=DEFAULT_LSTM_DIR):
    mkdir(out_dir)
    _safe_print("=" * 66)
    _safe_print("  CycleSlipNet: Three-model DL GNSS cycle-slip detection")
    _safe_print(f"  CI: Bootstrap percentile (B={BOOTSTRAP_B}, seed=12345)")
    _safe_print(f"  Detection threshold: {DL_CONF_THRESHOLD} x wavelength")
    _safe_print("=" * 66)

    has_rinex = bool(rinex_path and os.path.isfile(rinex_path))
    odict = elist = None
    cm_map = {}; dt_sec = fd

    if has_rinex:
        odict, elist, sys_obs, hdr, dt_from_rinex = parse_rinex3_obs(rinex_path)
        cm_map  = resolve_obs(sys_obs)
        dt_sec  = dt_sec or dt_from_rinex
        _safe_print(f"  [Sampling interval] {dt_sec:.3f} s")
    else:
        _safe_print("  RINEX not found -> using slip.out data only")
        dt_sec = dt_sec or 1.0

    df_raw = parse_slip_out(slip_path, dt_sec=dt_sec)

    _safe_print("  -- Loading DL models --")
    dl_model_dict = {}
    for mname, mdir in [("TF", tf_dir), ("CNN", cnn_dir), ("LSTM", lstm_dir)]:
        if not mdir or not os.path.isdir(mdir): continue
        try:
            dl_model_dict[mname] = load_dl_model(mdir, mname)
        except Exception as e:
            _safe_print(f"  [{mname}] Load failed: {e}")
    if not dl_model_dict:
        raise RuntimeError("No DL models loaded — check model paths")

    sc_list     = list(SCENARIOS.keys()) if scenario == "all" else [scenario]
    all_metrics = []

    for sc_n in sc_list:
        if sc_n not in SCENARIOS: continue
        ms = run_scenario(
            sc_n, SCENARIOS[sc_n], df_raw, odict, elist, cm_map,
            dl_model_dict, seed, dt_sec, out_dir)
        all_metrics.extend(ms)

    global_dir = mkdir(os.path.join(out_dir, "ALL_summary"))
    if all_metrics:
        gdf = pd.DataFrame(all_metrics)
        gdf.to_csv(
            os.path.join(global_dir, "ALL_scenarios_all_models.csv"), index=False)
        summarize_mc_metrics_with_ci(
            all_metrics, global_dir, tag_prefix="ALL_scenarios")
        per_sc_dir = mkdir(os.path.join(global_dir, "per_scenario"))
        if "scenario" in gdf.columns:
            for sc_n in gdf["scenario"].unique():
                sc_sub = gdf[gdf["scenario"] == sc_n].to_dict("records")
                summarize_mc_metrics_with_ci(sc_sub, per_sc_dir, tag_prefix=sc_n)

    _safe_print(f"\n[Done] -> {out_dir}")
    return all_metrics


# ==============================================================
# §14  CLI entry point
# ==============================================================
if __name__ == "__main__":
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="CycleSlipNet: Three-model DL GNSS cycle-slip detection")
    p.add_argument("--slip",           default=DEFAULT_SLIP_PATH)
    p.add_argument("--rinex",          default=DEFAULT_RINEX_PATH)
    p.add_argument("--tf",             default=DEFAULT_TF_DIR)
    p.add_argument("--cnn",            default=DEFAULT_CNN_DIR)
    p.add_argument("--lstm",           default=DEFAULT_LSTM_DIR)
    p.add_argument("--output",         default=DEFAULT_OUTPUT_DIR)
    p.add_argument("--seed",           type=int,   default=DEFAULT_SEED)
    p.add_argument("--interval",       type=float, default=None,
                   help="Sampling interval in seconds (auto-detected from RINEX if omitted)")
    p.add_argument("--scenario",       default=DEFAULT_SCENARIO,
                   help="all / sc_A / sc_B / sc_C / sc_D")
    p.add_argument("--mc_n",           type=int,   default=MONTE_CARLO_N,
                   help="Monte Carlo runs per scenario (recommend >= 30)")
    p.add_argument("--bootstrap_b",    type=int,   default=BOOTSTRAP_B,
                   help="Bootstrap resamples B (recommend >= 1000)")
    p.add_argument("--conf_threshold", type=float, default=DL_CONF_THRESHOLD,
                   help="Detection threshold in wavelength multiples")
    a = p.parse_args()
    MONTE_CARLO_N     = a.mc_n
    BOOTSTRAP_B       = a.bootstrap_b
    DL_CONF_THRESHOLD = a.conf_threshold
    run(a.slip, a.rinex, a.tf, a.cnn, a.output,
        seed=a.seed,
        fd=a.interval if a.interval else DEFAULT_D,
        scenario=a.scenario, lstm_dir=a.lstm)