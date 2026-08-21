# Claude Development Guide

This file contains instructions and context for Claude Code to help with development and maintenance of the QSPIO PDAC repository.

## Repository Overview

QSPIO PDAC is a quantitative systems pharmacology (QSP) modeling framework specifically designed for pancreatic ductal adenocarcinoma (PDAC) immunotherapy simulation and analysis. This repository extends the core QSPIO framework with PDAC-specific implementations and parameters. The codebase is MATLAB-based with SimBiology models.

### Where things are since the 2026-08-19 split

**This repository is the model.** The corpus, the priors, the elicitation and the
virtual-population campaign moved to `pdac-paper`, and the SBI/NPE workflows that
preceded the campaign were removed rather than moved -- they depended on the
deprecated `qsp_inference` and are preserved in the frozen `jeliason/pdac-build`
mirror.

So paths named in the mechanism notes below -- `calibration_targets/` (now
`corpus/`), `notes/calibration/`, `parameters/pdac_priors.csv`,
`submodel_targets/` -- refer to **the pdac-paper checkout**, not to this one.

`notes/` was removed here on the same day, so citations of `notes/model/...`
resolve in the frozen **`jeliason/pdac-build`** mirror. The design record they
point at is preserved there in full: the mechanism reviews, the ECM and TLA
investigations, the GPT-5 audit transcripts and the LaTeX audit documents.

All of these are cited as provenance for a modelling decision, and the
citations are left as they were written rather than rewritten in bulk. The sync checks that
compare the two sides (`check-priors-sync`, `check-cpp-sync`,
`check-observables-sync`) are run from pdac-paper with `MODEL_ROOT` pointing
here.

## Key Commands

### Testing
```bash
# Run Python tests (use .venv)
.venv/bin/pytest tests/

# Run slow tests (includes MATLAB simulations)
.venv/bin/pytest tests/ -m slow
```

### Setup
```matlab
% Add paths and configure environment
startup
```

```bash
# Python environment (IMPORTANT: always use .venv for Python commands)
source .venv/bin/activate  # Activate virtual environment
# OR directly use: .venv/bin/python <script.py>
```

## Project Structure

- `core/` - Core MATLAB functionality (shared framework)
  - `modules/` - SimBiology model modules for different therapies (including GVAX)
  - `visualization/` - the D2 schematic source and its ELK renderer
  - `postprocessing/` - `howManyClones` and `simulation_config`, the two the modules and tests call
  - `model_compilation/` - Model compilation utilities
  - `utils/` - Helper functions
- `scripts/` - model assembly, export and diagnostics
  - `immune_oncology_model_PDAC.m` - the main model file
  - `export_sbml.m` - what `make export-sbml` builds from
  - `evolve_to_diagnosis.m`, `set_healthy_populations.m` - the diagnosis-state solve
  - `trace_*.m`, `diag_*.m`, `smoke_*.m`, `toy_*.py` - mechanism traces and probes
  - `run_median_simulation.m` - reads the priors from pdac-paper via `paper_root`
  - `hpc_status.py` - cluster/local git drift report
- `cpp/` - the generated C++ simulator (`cpp/qsp/ode`) and its CMake build (`cpp/sim`)
- `resources/cpp/param_all.xml` - the parameter XML the simulator reads
- `parameters/` - `parameters_PDAC.m` (the model's own defaults), posterior medians
- `scenarios/` - Treatment scenario YAML configurations
- `notes/model/` - mechanism audits, structural reviews, the design record
- `tests/` - the model's own tests: MATLAB<->C++ parity and
  evolve-to-diagnosis. Corpus-facing qualitative and reduction checks live in
  pdac-paper; the obsolete classifier-coverage suite was retired at the split.
- `notes/` - Working notes and planning documents
- `presentations/` - Marp markdown slides (build with `npx @marp-team/marp-cli --no-stdin`)
- `tests/` - Python test suite (pytest)
- `startup.m` - MATLAB environment setup

## Common Development Tasks

### Iterative model debugging — read the runbook first

Before adjusting parameters, calibration targets, submodel targets, or kicking
off SBI verification runs, read **`notes/workflows/iterative_model_debugging.md`**.
It covers the recurring touchpoints (worktree setup with `cpp/sim/build` symlink
and `.compare_cache/` copy, C++ rebuild via `make refresh-cpp` /
`cmake --build cpp/sim/build`, single-theta `qsp_sim` traces with hand-edited
XMLs, `submodel_priors.yaml` refresh via
`scripts/regen_submodel_priors.py [--invalidate PARAM ...]`, the
`python -m workflows.sbi_runner --with-submodel-sync ...` invocation, and how
to read joint reachability + pairwise tension from the prior-PC log).

### PDAC-Specific Model Development
1. Modify PDAC parameters in `parameters/`
2. Update PDAC model configuration in `scripts/immune_oncology_model_PDAC.m`
3. Test changes with `scripts/single_patient_validation.m`

### Adding New Therapy Modules
1. Create new module in `core/modules/`
2. Update `core/modules/gvax_module.m` for PDAC-specific implementations
3. Follow existing module naming conventions
4. Add appropriate tests in `scripts/test_*.m`
5. Document mechanism in `notes/model/pdac_changelog.md`

### Adding New Mechanistic Modules
When adding new biological mechanisms (e.g., chemokine pathways, cell recruitment):
1. Create a self-contained module in `core/modules/`
2. Add parameters to `parameters/parameters_PDAC.m` with `.derivedFrom = {}`
3. Add module call in `scripts/immune_oncology_model_PDAC.m`
4. Add parameters to priors CSV if they should be varied in SBI
5. Document in `notes/model/pdac_changelog.md` with literature references

### Visualization and Analysis
- RECIST response evaluation tools
- Survival analysis capabilities

## Development Guidelines

### Code Style
- Follow existing MATLAB conventions
- Use descriptive function and variable names
- Maintain modular architecture for project extension
- Design for reusability across cancer types

### Testing
- Run `pytest tests/` before committing changes
- Test module integration with core framework
- Verify compatibility with existing modules

### Module Design
- Create self-contained therapy modules
- Design for PDAC-specific tumor microenvironment
- Document module interfaces and dependencies
- Ensure compatibility with existing PDAC parameters

### Module Assembly Caveats
The macrophage module (`macrophage_module.m`) modifies reactions from earlier modules at assembly time:
- **Appends** `*H_IL10` to Treg-mediated CD8 exhaustion reaction rate
- Creates IL-12, IL-10, IFNg species and Hill functions used by APC_module's maturation formula

The APC module (`APC_module.m`) splits dendritic cells into two subtypes: **cDC1** (BATF3+, cross-presenting, primes CD8 via MHC-I) and **cDC2** (IRF4+, MHC-II presenting, primes CD4 Th/Treg). Species: `cDC1`, `cDC2`, `mcDC1`, `mcDC2` in both V_T and V_LN (8 total). Homeostatic densities are **independently parameterized** per subtype and compartment (`APC0_cDC1_T`, `APC0_cDC2_T`, `APC0_cDC1_LN`, `APC0_cDC2_LN`) — the old `f_cDC1 × APC0_T` parameterization was replaced because the locked ratio prevented SBI from matching both densities simultaneously (LOO D² dominated by the cDC1/cDC2 pair). Hill function routing: `H_mAPC` uses `V_LN.mcDC1` (CD8 priming), `H_APC_Treg` and `H_APC_Th` use `V_LN.mcDC2` (Treg/Th priming). Both CD4 Hill functions use `V_LN.nCD4` (shared naive CD4+ pool) and `n_CD4_clones` (self-antigen clone diversity) in their denominators. GVAX-matured APCs route to `V_LN.mcDC2` (GM-CSF recruits cDC2/moDCs). Aggregate alias rules (`APC_total_T`, `mAPC_total_T`, `mAPC_total_LN`) sum subtypes for backward compatibility. IL-12 secretion uses `V_T.mcDC1` (not total mAPC). Antigen uptake uses `APC_total_T` (all immature DCs).

APC maturation logic (DAMP+IL-12 OR, IL-10+TGFb inhibition) is self-contained in `APC_module.m` with subtype-specific rates (`k_APC_mat_cDC1`, `k_APC_mat_cDC2`). The macrophage module provides the cytokine Hill functions (`H_IL12`, `H_IL10`) that appear in both cDC1 and cDC2 maturation rate expressions, but does not modify them.

**Th exhaustion routes** (`Th_module.m`): conventional Th cells transition to `Th_exh` via three parallel mechanisms mirroring CD8 in `Teff_module.m`:
- chronic-antigen baseline: `k_Th_basal_exh × Th × C_total/(C_total+K_C_exh+cell)` (shares `K_C_exh` with CD8 baseline route)
- TGFβ-driven: `k_TGFb_Th_exh × Th × H_TGFb_Th_exh` where `H_TGFb_Th_exh = V_T.TGFb / (V_T.TGFb + TGFb_50_Th_exh)` is a cell-specific TGFβ Hill (defined in `cytokine_tcell_module.m` next to `H_TGFb_Teff`)
- antigen-clearance funnel (post-tumor-eradication): `k_cell_clear × Th × K_C2_recruit/(C_total² + K_C2_recruit)`
Regular Th death is `k_CD4_death × Th → null` (NOT routed to Th_exh — the prior `V_T.Th -> V_T.Th_exh` at `k_CD4_death` was a death-relabel that locked Th_exh fraction to `k_CD4_death/(k_CD4_death+k_cell_clear)`, which made matching Shibuya 2014's PD-1+ CD4 fraction impossible without crashing CD4 density past the `k_CD4_death > k_Th_pro_tumor` bifurcation).
Th_exh decay uses `k_cell_clear × Th_exh → null` (in `simbio_init.m`). CD8_exh decay is NOT via `k_cell_clear` — it uses the purpose-built, submodel-anchored `k_CD8_exh_death × CD8_exh → null` reaction in `Teff_module.m` (`k_CD8_exh_death_PDAC_deriv001.yaml`). A redundant `k_cell_clear` clearance of CD8_exh formerly also lived in `simbio_init.m`, double-counting exhausted-CD8 death at `(k_cell_clear + k_CD8_exh_death)` ≈ 0.25/day instead of the calibrated 0.15/day (suppressing the exhausted-CD8 fraction, gvax_nivo d21 median ~0.20 → ~0.25 single-decay); it was removed. Submodel targets in `calibration_targets/submodel_targets/k_Th_basal_exh_*.yaml` (Lawton 2024 chronic CD3/CD28) and `k_TGFb_Th_exh_TGFb_50_Th_exh_*.yaml` (Park 2016 TGFβ dose-response on PD-1+ CD4).

The APC module uses `H_pore_APC` (pore-size-based APC exclusion) defined in `fibroblast_module_enhanced.m`. Since APC_module is called before fibroblast_module_enhanced in model assembly, this is a forward reference resolved at simulation time by SimBiology.

**TLA / urelumab module** (`tla_module.m`): owns the entire urelumab (anti-CD137 / 4-1BB agonist) effect and the vaccine-induced intratumoral tertiary lymphoid aggregate (TLA) it acts in. Assemble AFTER `gvax_module` (needs `V_ID.mAPC`), `Teff_module` (needs `V_T.CD8`, `CD8_pex`), and `pk_module('aCD137')`. Tier-1 design — the TLA and its CD8 states are **tagged V_T sub-pools, not a new compartment**, seeded only by the GVAX depot, so baseline/Arm-A/Arm-B fits are untouched (all TLA species are zero without vaccine; `CD8_vic ≈ 0` without urelumab). Species: `CD8_TLA` (resting within-TLA), `CD8_TLA_act` (CD137+ activated, the Heumann Fig-5b observable), and `CD8_vic` (egressed neoplastic-vicinity effector-memory pool, Fig-5d). `TLA` (niche fraction, 0–1) and `CXCL13` are dimensionless/quasi-steady **parameters with rules**, not species (CXCL13 has a sub-mm diffusion length, so it is niche-local, not well-mixed in V_T — the feed-forward source is diluted in `V_niche = TLA·V_T`). Two urelumab antitumor levers: (1) the small `CD8_vic` pool kills via a **separate** `k_C_vic_eff` reaction that drops the `(1−H_pore_T)` pore throttle and Treg-saturation factor and carries the `phi_vic_kill` quality multiplier — quality not quantity, and throttle-bypassing because vicinity effectors are positioned at the neoplastic interface; (2) `k_CD137_reinvig` refluxes the large progenitor-exhausted `CD8_pex` (Tpex) reservoir back to functional effectors, gated by `H_CD137_T` and synergizing with the anti-PD-1 `k_reinvig` arm — the dominant lever. **Assembly-time modification:** like `macrophage_module`, the module rewrites an existing `Teff_module` reaction in place — the bulk-CD8 logistic denominator `V_T.CD8/(rho_CD8_max*V_T + cell)` becomes `CD8_T_tot/(...)` where `CD8_T_tot = V_T.CD8 + CD8_TLA + CD8_TLA_act + CD8_vic`, so niche/vicinity sub-pools still count against carrying capacity (the proliferation prefactor stays bulk `V_T.CD8`). Within-TLA CD8 deaths were removed (folded into the `CD8_T_tot` logistic); vicinity killing also appends `+k_C_vic_eff*V_T.C1` to the `R_Tcell` antigen-release accumulator. **CD137 activation is LATCHED, not occupancy-tracked** (2026-07-19): the `CD8_TLA → CD8_TLA_act` conversion (reaction 8) is driven by a durable exposure latch `L_CD137` (a parameter with a rate rule, `dL_CD137/dt = k_CD137_latch_on·H_CD137_T·(1−L_CD137) − k_CD137_latch_off·L_CD137`), **not** instantaneous `H_CD137_T`. Urelumab is a single day-7 dose that washes out (aCD137 t½ ~5.5 d) well before the d21 mIHC readout, so occupancy-gated activation decayed the within-TLA CD137+ fraction to ~0.41 (peak ~0.67 at d10) vs Heumann Arm C's ~0.83; the slow-decaying latch (`k_CD137_latch_off` ~0.015/day) keeps the activated state up AND keeps activating fresh niche CD8 (which otherwise accumulate in the resting `CD8_TLA` denominator — it grows ~1.7M→5.3M over the window — and dilute the fraction). Zero without urelumab (`H_CD137_T=0 → L_CD137=0 →` `f_act_base` floor), so Arms A/B are untouched. `H_CD137_T` still *directly* gates the `k_CD137_reinvig`/`k_CD137_prolif` Tpex arm (instantaneous is correct there — those are during-drug fluxes, not a latched marker). See `notes/model/tla_cxcl13_urelumab_investigation.md`.

**Vaccine-driven CD8 expansion (cd8_fc Paths A, B & C)** — mechanisms that make the GVAX-induced intratumoral CD8 fold-change treatment-differential. (1) **Path B — CD40L-licensed division destiny** (`cytokine_tcell_module.m`): the CD8 division-destiny rule's costim term `N_div_costim*H_CD28_APC` becomes `N_div_costim*H_CD28_APC*H_div_help`. `H_CD28_APC` (the B7–CD28 synapse Hill) sits saturated (~0.98) and treatment-invariant at the cold-tumor baseline, so on its own the costim divisions fire as hard at diagnosis as under vaccine. The licensing gate is **`H_div_help`** (a DEDICATED Hill in `Th_module.m`, `V_LN.Th/(V_LN.Th+Th_50_div)`), **not** the priming-throughput `H_CD4_help` — a trace showed `H_CD4_help` does double duty (it also gates cDC1→CD8 priming) and sits **saturated** (~0.999), so it cannot serve as a vaccine-cold costim gate. `H_div_help` reads `V_LN.Th` against a higher EC50 (`Th_50_div`≈1e6) so it is genuinely vaccine-cold (~0.61→0.92), making the costim divisions a treatment differential. `N_aTh` (CD4) keeps the unlicensed `H_CD28_APC` form. **Critical dependency:** `H_div_help` only swings the right way because the CD4 help arm now responds to vaccine — see Path A's CD4 arm; without it `V_LN.Th` *drops* under vaccine (0.76×) and the gate moves backwards. A caveat note in `checkpoint_module.m` flags the cleaner long-term form (make `H_CD28_APC`'s B7 itself licensing-dependent). (2) **Path A — epitope-spreading clonal recruitment** (`gvax_module.m`, assembly-time rewrite like `macrophage_module`/`tla_module`): applied to **both** the CD8 and CD4 commitment reactions. The naive-CD8→`aCD8` reaction (`H_CD4_help*n_CD8_clones` substring) has `n_CD8_clones` rewritten to `n_CD8_clones*(1 + f_epitope_spread_CD8*H_P1_spread)`; the naive-CD4→`aTh` reaction (`H_P1_II*V_LN.nCD4*n_CD4_clones` substring) has `n_CD4_clones` rewritten to `n_CD4_clones*(1 + f_epitope_spread_CD4*H_P1_spread)`. `H_P1_spread` is the P1_GVAX-concentration Hill (reuses `EC50_P1_mature`, =0 at diagnosis). The clone multiplier sits on the *activated*-pool production, not the naive-pool drain, so it boosts output without accelerating naive depletion. **Spreading is asymmetric — two SBI params:** `f_epitope_spread_CD4` (median 6) ≥ `f_epitope_spread_CD8` (median 4), because CD4 broadening uses the promiscuous MHC-II default pathway while CD8 broadening is rate-limited by cDC1 cross-presentation and is downstream of CD4 help (Thomas 2004 cross-priming). The `n_CD8_clones`/`n_CD4_clones` *parameters* are unchanged, so `antigen_rate.m` and `TCR_KPR_module.m` usages are untouched. The CD4 arm is load-bearing: it un-mutes the help response (`V_LN.aTh` 1.19×→~6×, flips `V_LN.Th` from dropping to rising ~7×), which is what makes Path B's `H_div_help` fire. (3) **Path C — in-situ intratumoral proliferation** (`tla_module.m`, 2026-07-18, added as a new reaction — NOT an assembly-time rewrite): a SEPARATE saturating CD8 proliferation source `null -> V_T.CD8` at `k_CD8_T_vac_pro*V_T.CD8*H_P1_spread/(1 + CD8_T_tot/(rho_CD8_max*V_T + cell))`, addressing the LN→tumor translation gap Paths A/B left open. Paths A/B lift LN priming/production, but that thin influx dilutes into the large, locally-proliferated tumor CD8 pool (the trace: `V_T.CD8` ~22× circulating `V_C.CD8`, so tumor CD8 is proliferation-dominated, not influx-limited) — nothing trafficking-side moves cd8_fc, only `k_CD8_T_pro` does (at −0.48, negative: it raises baseline and shrinks the fold). The constitutive `k_CD8_T_pro` (Blessin-anchored) is vaccine-blind, so it fires equally at diagnosis and d21; Path C is the vaccine-DIFFERENTIAL increment to that same in-situ proliferation, gated by `H_P1_spread` (the GVAX antigen Hill, identically 0 at diagnosis / no-vaccine arms) so it decouples cd8_fc from the baseline `cd8_density` that pins `k_CD8_T_pro`. Biology: GVAX builds intratumoral TLAs with Ki-67+ CD8 (Lutz 2014, Cancer Immunol Res 2:616 — same Zheng group as the cd8_fc data); the aggregate is the histological signature of the vaccine's proliferative program, but Path C proliferates the BULK `V_T.CD8` (the non-LA countable pool, `CD8_total_T = V_T.CD8+exh+pex+vic`), NOT the `CD8_TLA` sub-pool — the within-TLA pool is only ~`a_TLA` (~0.55%) of bulk, far too small a reservoir to source a bulk-scale fold at any Ki-67-anchored division rate. Shares the constitutive term's `CD8_T_tot` density brake so the two jointly saturate at carrying capacity. `k_CD8_T_vac_pro` is a PLACEHOLDER (CSV prior, median 0.5/day, σ 1.5), to be anchored to the Lutz Ki-67-in-TLA fraction once SBI confirms the pathway. All three paths are vaccine-gated (≈0 at diagnosis) so baseline / Arm-A / Arm-B fits are preserved. Before Path C, at submodel medians cd8_fc moved ~1.6×→2.43× with the residual gap in the LN→tumor translation (`V_C.CD8` ~4.4× diluting to `V_T.CD8` ~2.4×). See `notes/calibration/cd8_fc_vaccine_expansion_deficit.md`.

The fibroblast module (`fibroblast_module_enhanced.m`) manages three CAF subtypes: `iCAF` (inflammatory, secretes IL-6/CXCL12), `myCAF` (myofibroblastic, secretes collagen/TGFb/VEGF), and `apCAF` (antigen-presenting, sequesters P1 and induces Tregs). The old single `aPSC` species no longer exists. The module tracks `V_T.collagen` (mass in mg) and derives `phi_collagen = V_T.collagen / (V_T × rho_collagen)` — the collagen volume fraction (dimensionless, 0–1). `phi_collagen` is the primary input for all ECM-dependent Hill functions. The volume rule appends `+V_T.collagen/rho_collagen` to V_T, creating a self-limiting ceiling as collagen fills the available volume. Secretion saturation uses `(1 - phi_collagen)` which engages as collagen approaches physical packing limits. See `notes/model/ecm_formulation_review.md` for the full ECM module design rationale.

**IL-1 / iCAF mechanism:** `V_T.IL1` (IL-1alpha, nM) is secreted by cancer cells and drives iCAF formation via `H_IL1_eff = IL1 / (IL1_50 * (1 + TGFb/TGFb_50_IL1R1) + IL1)`. TGFb antagonizes IL-1 signaling by shifting the effective EC50 rightward (IL-1R1 downregulation), not by directly suppressing the iCAF program. This replaces the old `(1-H_TGFb)` proxy which was mechanistically incorrect. The fibroblast module uses three separate TGFb EC50s: `TGFb_50_CAF_act` (activation/interconversion), `TGFb_50_col_sec` (collagen secretion), `TGFb_50_IL1R1` (IL-1R1 repression) — no remaining generic `TGFb_50` references.

**CAF clearance and stiffness niche:** Activated CAF (iCAF/myCAF/apCAF) tumor-eradication clearance is modulated by `(1 - H_stiff_fib)` — stiff ECM suppresses clearance, allowing fibrotic scar persistence post-therapy. This is consistent with `K_stiff_niche × H_stiff_fib` in the carrying capacity.

**ECM-derived quantities** (all defined in `fibroblast_module_enhanced.m`):
- `d_pore = d_pore_ref × (phi_col_ref / phi_collagen)^0.5` — pore size between collagen fibers (μm)
- `stiffness = E_ref × (phi_collagen / phi_col_ref)^n_stiff` — tissue stiffness (dimensionless, kPa-scaled)
- `H_pore_T`, `H_pore_APC` — pore-size-based immune cell exclusion (replaces old H_ECM_T_mot, H_ECM_APC_mig)
- `H_stiff_Texh` — stiffness-mediated T cell exhaustion (replaces old H_ECM_T_exh)
- `H_stiff_fib` — stiffness-mediated fibroblast activation (YAP/TAZ mechanotransduction)

**CAF carrying capacity** has two terms: `f_stroma_max × C_total` (cancer paracrine) + `K_stiff_niche × H_stiff_fib` (stiffness-dependent niche). The stiffness term allows CAFs to persist in stiff ECM even without cancer cells (post-therapy fibrotic scar). Collagen degradation includes both constitutive (`k_col_deg`) and M2 macrophage-dependent (`k_MMP_M2 × V_T.Mac_M2 / V_T`) terms. M2-like TAMs are the dominant collagen degraders via mannose-receptor endocytosis (Madsen 2013) and TAM-mediated collagen turnover in PDAC (LaRue 2022).

The angiogenic factor species is named `VEGF` (renamed from the legacy `c_vas`). Rate constant parameter names (`k_vas_Csec`, `k_vas_deg`, `k_vas_Msec`, `k_vas_myCAF_sec`) retain the `vas` prefix for continuity.

The antigen release function (`antigen_rate.m`) searches for killing reactions by string-matching parameter names in reaction rates. If killing parameters are renamed, the search string must be updated or the killing-dependent antigen release term silently evaluates to zero.

**Scenario YAML requirements**: All scenario YAMLs that simulate from diagnosis must include `initialization_function: 'evolve_to_diagnosis'` in their `sim_config` section. Without this, `batch_worker.m` starts from default initial conditions (C1=1 cell) instead of diagnosis state, producing near-zero values. The `evolve_to_diagnosis` function evolves the model from healthy state to a target tumor diameter (extracted from the variant's `initial_tumour_diameter` parameter).

**Test statistics function signature**: All test statistics in `scenarios/test_statistics/*.csv` must use the 3-argument Pint-aware signature: `compute_test_statistic(time, species_dict, ureg)`. The `species_dict` values are Pint Quantities with units from `model_structure.json`. Functions returning dimensionless values should return `value * ureg.dimensionless`. The old 2-argument signature `(time, species_dict)` without Pint will silently produce NaN.

**Calibration/prediction target time-series reduction**: The reduction from an observable time-series to the scalar compared against the data is declared on the `observable:` block, and exactly one of the two forms is required — there is no implicit default, so omitting both is a hard error rather than a silent t=0.

- `readout_time: <float>` + `readout_time_unit: <str>` — linearly interpolate the series at that simulation time. This is the common case (e.g. a day-21 resection endpoint uses `readout_time: 21.0`, `readout_time_unit: day`). A diagnosis/baseline snapshot states `readout_time: 0.0` explicitly.
- `reduce_observable: |` — a `def reduce_observable(time, series) -> float` for peak/Cmax, AUC, final-value, or time-to-threshold measurements that have no single readout time.

The older `empirical_data.index_values` / `index_unit` / `index_type` fields (and the top-level equivalents on prediction targets) are gone. Cross-scenario targets do not use this reduction and keep their own per-arm `observable.inputs[*].observable_code` shape.

See `notes/model/model_additions_ranked.md` (Implementation Audit section) for full details.

See `notes/calibration/calibration_strategy.md` for calibration workflow, target coverage analysis, and priority targets to add. See `notes/calibration/fibroblast_param_audit.md` for fibroblast module parameter coverage. See `notes/calibration/calibration_sufficiency_analysis.md` for calibration norms benchmarking (vs published QSP-IO models), submodel target selection rationale, and audit of `to-review/` old-schema targets.

## PDAC Model Integration

### PDAC Module Types
- **Therapy modules**: GVAX, checkpoint inhibitors for PDAC
- **Visualization modules**: PDAC-specific plotting and analysis functions
- **Parameter modules**: PDAC-specific parameter distributions
- **Batch processing**: SSH cluster integration for large-scale analysis

## Useful File Patterns

### Therapy Modules
- `core/modules/*_module.m` - Therapy-specific modules (GVAX, checkpoint inhibitors, etc.)

### Configuration Files
- `parameters/*.m` - PDAC parameter files
- `pdac-paper`'s `parameters/pdac_priors.csv` - the prior master, kept in step with `parameters_PDAC.m` by `make check-priors-sync` (run from there)
- `scenarios/*.yaml` - Treatment scenario configurations
  - `baseline_no_treatment.yaml` - Diagnosis snapshot (stop_time=1, evolve_to_diagnosis)
  - `gvax_neoadjuvant_zheng2022.yaml` - Neoadjuvant GVAX monotherapy (stop_time=21, evolve_to_diagnosis)
  - `gvax_nivo_neoadjuvant_zheng2022.yaml` - Neoadjuvant GVAX+nivo (stop_time=21, evolve_to_diagnosis)
  - `clinical_progression.yaml` - 90-day progression for clinical outcome calibration
- `model_structure.json` - Species metadata (names, compartments, units) for Pint conversions in test statistics

### Core Framework Files
- `simbio_*.m` - Model compilation utilities
- `startup.m` - Environment configuration
- Functions in `core/utils/` - General utilities

## PDAC-Specific Workflows

### Paper-facing behavior checks

Calibration bounds, prior/design-space coverage, and model-reduction ablations
belong to the pdac-paper checkout because their meaning depends on its corpus
and priors. This repository keeps only tests of the assembled model itself.

# In bayesflow_workflow.py or sbi_workflow.py
MODEL_VERSION = 'baseline_gvax'  # Descriptive version name
MODEL_DESCRIPTION = 'PDAC baseline: 8 params, 12 observables, GVAX therapy'
NUM_TOTAL_SIMULATIONS = 5000  # Will accumulate/reuse from pool
```

**When to bump MODEL_VERSION:**
- ✅ Priors change (parameter ranges, distributions)
- ✅ Test statistics change (new/removed observables)
- ✅ Model equations change (affects simulation outputs)
- ✅ Therapy modules change (dosing, mechanisms)

**Don't bump for:**
- ❌ Visualization/plotting changes
- ❌ Analysis code changes
- ❌ Documentation updates
- ❌ Comments, formatting
- ❌ Bug fixes that don't affect outputs

#### Usage Examples

**Running workflows:**
```bash
# First run: generates 2000 new simulations
python workflows/sbi_workflow.py

# Second run with same version: reuses all 2000 simulations
python workflows/sbi_workflow.py

# Third run requesting more: reuses 2000, generates 3000 new
# (Edit NUM_TOTAL_SIMULATIONS = 5000 first)
python workflows/sbi_workflow.py
```

#### How It Works

1. **Configuration Hashing**: Pools are identified by hash of priors CSV content, test stats CSV content, model script, and model version
2. **Batch Accumulation**: Each simulation run adds a new batch to the pool
3. **Automatic Reuse**: When requesting N simulations, the system:
   - Checks pool for available simulations
   - If enough exist, loads from pool (no MATLAB execution)
   - If not enough, runs only the additional simulations needed
   - Adds new batch to pool for future reuse
4. **Random Sampling**: If pool has more simulations than requested, randomly samples subset

#### Benefits

- **Efficiency**: Reuse expensive simulations across runs
- **Flexibility**: Change NUM_TOTAL_SIMULATIONS or SEED without invalidating cache
- **Iteration**: Accumulate simulations over multiple runs
- **Transparency**: Manifest tracks all batches with metadata
- **Self-Documenting**: Pool names include version and config hash

#### Implementation Files

- Simulation pooling and QSP simulator provided by the `maple` package
- `workflows/sbi_workflow.py` - Single-scenario SBI workflow with pooling
- `workflows/multi_scenario_sbi_workflow.py` - Multi-scenario SBI workflow with pooling

**Posterior predictive caching**: `simulate_with_parameters()` also caches results in suffix pools keyed by a hash of the theta matrix (`posterior_predictive_{hash}`). Repeated calls with the same posterior draws hit the cache instantly, avoiding redundant MATLAB runs during iterative development.
