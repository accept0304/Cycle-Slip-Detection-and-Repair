# Cycle-Slip Detection and Repair

Three-model deep learning framework for GNSS carrier-phase cycle-slip detection and repair.

**Models:** Transformer · CNN · LSTM  
**Input:** 1-D sliding window of raw TDCP\_norm values  
**Evaluation:** Monte Carlo injection · Bootstrap percentile 95% CI · McNemar test

---

## Installation

```bash
pip install -r requirements.txt
```

Requires Python ≥ 3.9, TensorFlow ≥ 2.12.

---

## Usage

```bash
python cycle_slip_net.py \
  --slip     /path/to/slip_normalized.out \
  --rinex    /path/to/station.24o \
  --tf       /path/to/transformer_dir \
  --cnn      /path/to/cnn_dir \
  --lstm     /path/to/lstm_dir \
  --output   /path/to/output \
  --scenario all
```

Each model directory must contain a `.keras` model file, `label_scaler.pkl`, and `window_size.pkl`.

---

## Evaluation Scenarios

| ID | Injection ratio | Slip magnitude |
|---|---|---|
| sc_A | 0–10 % | 1–10 cycles |
| sc_B | 0–10 % | 1–50 cycles |
| sc_C | 10–20 % | 1–50 cycles |
| sc_D | 20–30 % | 1–100 cycles |

---

## Output Structure

```
output/
├── sc_A/
│   ├── dl_only/          # per-run P/R/F1/RSR/MRE CSVs
│   ├── inject_logs/      # per-epoch injection records
│   ├── slip_counts/      # injected slip statistics
│   └── summary/
│       ├── ALL_models_metrics.csv
│       ├── metrics_with_CI_sc_A.csv   # Bootstrap 95% CI
│       ├── mcnemar_cross_model.csv
│       └── slip_size_stats_all.csv
├── sc_B/  sc_C/  sc_D/
└── ALL_summary/
    ├── ALL_scenarios_all_models.csv
    ├── metrics_with_CI_ALL_scenarios.csv
    └── per_scenario/
```

**Reported metrics:** Precision · Recall · F1 · Accuracy · RSR (Repair Success Rate) · MRE (Mean Repair Error in cycles)
