V5r competitor comparison package
================================

The original V5r detector functions were split directly from:
run_v5r_validation-2.m

Main comparison:
  run_competitor_comparison.m

Run:
  run_competitor_comparison('C:\path\to\MITBIH_60s_all_48_records')

Files:
  run_competitor_comparison.m
  run_v5r_validation.m                 (original validation split)
  TMA_PLUS_MANUAL_V5r.m                (exact original V5r core)
  TMA_PLUS_MANUAL_V5r_local.m
  tma_edge_safe.m
  tma_edge_safe_local.m
  chain_merge.m
  chain_merge_local.m
  apply_min_distance.m
  apply_min_distance_local.m
  event_shape.m
  event_shape_local.m
  estimate_W_fast.m
  estimate_W_fast_local.m
  local_peaks.m
  local_peaks_local.m
  load_record.m
  load_record_same_format.m
  match_beats.m
  match_beats_same.m
  pan_tompkins_sedghamiz_style.m
  cwt_ricker_detector.m

Important:
The V5r implementation is taken directly from the uploaded original validation
file. Pan-Tompkins and CWT-Ricker are toolbox-free comparison implementations;
they should be described as implementations/styles rather than claimed verbatim
reproductions of a particular published source unless independently verified.
