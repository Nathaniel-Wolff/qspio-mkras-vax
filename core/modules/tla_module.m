function model = tla_module(model, params)
% TLA_MODULE  Tertiary lymphoid aggregate (TLA) niche + urelumab/CD137 mechanism.
%
% Owns the entire urelumab (anti-CD137, 4-1BB agonist) effect and the
% vaccine-induced intratumoral TLA compartment it acts in. Replaces the retired
% LN division-destiny term (N_div_CD137*H_CD137 in cytokine_tcell_module).
%
% Motivation (Heumann 2023, Nat Commun 14:3650, NCT02451982 Arm C). The
% urelumab signal is NOT a bulk CD8 expansion — it is an activation-STATE
% conversion of CD8 resident inside vaccine-induced intratumoral TLAs:
%   - CD137+ fraction of within-TLA CD8: ~10% (Arm A, GVAX) / ~0% (Arm B,
%     GVAX+nivo) / ~75% (Arm C, +urelumab)            [Fig 5b]
%   - general within-TLA CD8 NUMBER is flat across arms [Suppl Fig 9]
%   - bulk cytotoxic (GZMB+) fraction unchanged          [Fig 5c]
%   - a small activated subset (2-4%) appears in the vicinity OUTSIDE the TLA
%     only in Arm C (egress)                             [Fig 5d]
%
% Design (Tier 1 — tagged V_T sub-pools, not a new compartment). The aggregate
% and the CD8 states inside it are sub-pools of V_T, seeded ONLY by the GVAX
% depot, so the calibrated baseline / Arm-A / Arm-B fit is untouched (all TLA
% species are zero without vaccine).
%
% State:
%   CD8_TLA       cell  V_T species; resident within-TLA CD8 (resting / CD137-)
%   CD8_TLA_act   cell  V_T species; CD137+ activated within-TLA CD8 (Heumann obs)
%   CD8_vic       cell  V_T species; egressed activated effector in the neoplastic
%                       vicinity (Fig 5d 2-4% outside-TLA subset); higher killing
%                       quality (phi_vic_kill) + durable (k_vic_decay) — the
%                       urelumab antitumor/survival effector pool
%   TLA           []    dimensionless PARAMETER with a rate rule (abstract niche
%                       fraction [0,1); not a physical pool -> not a species)
%   CXCL13        nM    PARAMETER (quasi-steady repeatedAssignment, not a species/
%                       ODE) — fast chemokine turnover; avoids a stiff fast mode
%
% Rules:
%   H_CD137_T = V_T.aCD137/(V_T.aCD137 + CD137_50)            (tumor urelumab Hill)
%   CXCL13    = (k_CXCL13_sec*V_ID.mAPC + k_CXCL13_ff*V_T.CD8_TLA)/(k_CXCL13_deg*V_T)
%   H_CXCL13  = CXCL13/(CXCL13 + CXCL13_50_TLA)
%   d(TLA)/dt = k_TLA_form*H_CXCL13*(1-TLA) - k_TLA_decay*TLA (rate rule)
%   CD8_T_tot = V_T.CD8 + V_T.CD8_TLA + V_T.CD8_TLA_act       (carrying-cap denom)
%
% Reaction network:
%   (CXCL13 source/feed-forward/clearance and TLA formation/decay are the
%    algebraic/rate rules above, not reactions.)
%   (6)  CD8 -> CD8_TLA                 k_TLA_seq * a_TLA * CD8       (residence in)
%   (7)  CD8_TLA -> CD8                 k_TLA_release * CD8_TLA       (residence out, resting)
%   (7b) CD8_TLA_act -> CD8             k_TLA_release * CD8_TLA_act   (residence out, activated; decouples # from act-fraction)
%   (8)  CD8_TLA -> CD8_TLA_act         k_CD137_act*(f_act_base+(1-f_act_base)*H_CD137_T)*CD8_TLA
%   (9)  CD8_TLA_act -> CD8_TLA         k_CD137_deact * CD8_TLA_act
%   (10) CD8_TLA_act -> CD8_vic         k_TLA_egress * CD8_TLA_act    (egress to vicinity effector pool)
%   (10b) CD8_vic -> CD8                k_vic_decay * CD8_vic         (vicinity effector durability/turnover)
%   (14) CD8_pex -> CD8                 k_CD137_reinvig * CD8_pex * H_CD137_T  (urelumab Tpex reinvigoration; dominant antitumor lever, synergizes with anti-PD-1 k_reinvig)
%   (deaths) within-TLA + vicinity turnover are folded into the CD8_T_tot logistic,
%            not separate k_CD8_death reactions (see note at bottom).
%   Killing coupling: CD8_vic enters Teff_module's k_C_Tcell_eff effector-density
%   factor with quality weight phi_vic_kill (vicinity effectors kill better per cell).
%
% with H_CXCL13 = CXCL13/(CXCL13 + CXCL13_50_TLA).
%
% Assemble AFTER gvax_module (needs V_ID.mAPC), Teff_module (V_T.CD8,
% k_CD8_death) and pk_module('aCD137') (V_T.aCD137).
%
% Observables (defined downstream, e.g. aggregate_observables_module):
%   within-TLA CD137+ fraction        = CD8_TLA_act / (CD8_TLA + CD8_TLA_act)   [Fig 5b]
%   CD8+CD137+ % of TLA ROI cells      = CD8_TLA_act / (TLA ROI cell count)      [Fig 3b]
%   CD137+ vicinity (egressed)         = egressed fraction of bulk V_T.CD8       [Fig 5d]
%   TLA prevalence                     = V_T.TLA above threshold                 [Li/Lutz]

V_T = model.Compartment(3);

%% ---- Species (V_T sub-pools) -------------------------------------------------
% (CXCL13 is NOT a species — it is a quasi-steady repeatedAssignment parameter
% below, to avoid a stiff fast chemokine mode at the GVAX dose event.)
CD8_TLA = addspecies(V_T, 'CD8_TLA', 0, 'InitialAmountUnits', 'cell', ...
    'Notes', 'Resting (CD137-) CD8 resident within the TLA');
CD8_TLA_act = addspecies(V_T, 'CD8_TLA_act', 0, 'InitialAmountUnits', 'cell', ...
    'Notes', 'CD137+ activated CD8 within the TLA (Heumann within-TLA CD8+CD137+ observable)');
CD8_vic = addspecies(V_T, 'CD8_vic', 0, 'InitialAmountUnits', 'cell', ...
    'Notes', 'Egressed activated effector CD8 in the neoplastic vicinity outside the TLA (Heumann Fig 5d 2-4% subset); higher per-cell killing quality (phi_vic_kill) and effector-memory-like durability (k_vic_decay) — the urelumab antitumor/survival effector pool');

%% ---- Parameters --------------------------------------------------------------
% CD137_50 is added by parameters block / shared; add here only if absent so the
% module is self-contained if assembled standalone.
addp = @(nm) add_param_if_absent(model, params, nm);
addp('CD137_50');
addp('k_CXCL13_sec');
addp('k_CXCL13_ff');
addp('k_CXCL13_deg');
addp('CXCL13_50_TLA');
addp('k_TLA_form');
addp('k_TLA_decay');
addp('k_TLA_seq');
addp('k_TLA_release');
addp('a_TLA');
addp('k_CD137_act');
addp('k_CD137_deact');
addp('f_act_base');
addp('k_TLA_egress');
addp('phi_vic_kill');
addp('k_vic_decay');
addp('k_CD137_reinvig');
addp('f_CD137_func');
addp('k_CD137_prolif');
addp('k_CD8_T_vac_pro');
addp('k_CD137_latch_on');
addp('k_CD137_latch_off');

%% ---- Rules -------------------------------------------------------------------
% Tumor-local urelumab engagement (was the LN division-destiny Hill; relocated to
% V_T here). =0 when urelumab undosed (V_T.aCD137=0) so all TLA conversion sits at
% the f_act_base floor in Arms A/B.
p = addparameter(model, 'H_CD137_T', 0, 'ValueUnits', 'dimensionless', ...
    'ConstantValue', false, 'Notes', 'Hill of tumor urelumab (V_T.aCD137) gating TLA CD8 CD137-activation');
addrule(model, 'H_CD137_T = V_T.aCD137/(V_T.aCD137 + CD137_50)', 'repeatedAssignment');

% Integrated-exposure LATCH for durable CD137/4-1BB activation (2026-07-19).
% Urelumab is a SINGLE day-7 dose; by the day-21 readout it has washed out to
% ~1/15 of peak (aCD137 t1/2 ~5.5 d), so the instantaneous occupancy H_CD137_T
% collapses (trace: 0.78 -> 0.19). Yet Heumann observes ~83% within-TLA CD8+CD137+
% at d21 (14 d after the dose) — a DURABLE latched activation program, not a live
% occupancy readout. Driving reaction 8 directly off H_CD137_T made the activated
% fraction peak ~0.67 at d10 then decay to ~0.41 at d21 (two causes: the drive
% collapses AND fresh CD8 keep entering the niche unactivated, growing the resting
% CD8_TLA denominator ~1.7M->5.3M). L_CD137 integrates the exposure and decays
% slowly (k_CD137_latch_off ~0.015/day, weeks-long memory), so activation persists
% after washout AND fresh niche cells keep getting activated. =0 without urelumab
% (H_CD137_T=0), so Arms A/B stay at the f_act_base floor. See
% notes/model/tla_cxcl13_urelumab_investigation.md.
p = addparameter(model, 'L_CD137', 0, 'ValueUnits', 'dimensionless', ...
    'ConstantValue', false, 'Notes', 'Integrated-exposure latch for durable CD137/4-1BB activation: charges from H_CD137_T during the urelumab window, decays slowly (k_CD137_latch_off) — drives the TLA CD8 activation-state conversion (reaction 8) so the CD137+ fraction persists to d21 after single-dose urelumab washes out');
addrule(model, 'L_CD137 = k_CD137_latch_on*H_CD137_T*(1 - L_CD137) - k_CD137_latch_off*L_CD137', 'rate');

% CXCL13 organizer — QUASI-STEADY (fast chemokine turnover). An algebraic
% repeatedAssignment parameter, not a species/ODE: CXCL13 equilibrates far faster
% than the niche, so integrating it as a stiff fast state only destabilizes the
% solver at the GVAX dose event. CXCL13 = (vaccine source from depot mature DCs +
% TLA-resident-CD8 feed-forward) / clearance, per V_T volume (units: nanomolarity
% = nanomole/cell/day * cell / (1/day * liter)).
p = addparameter(model, 'CXCL13', 0, 'ValueUnits', 'nanomolarity', ...
    'ConstantValue', false, 'Notes', 'CXCL13 organizer chemokine (quasi-steady; GVAX-depot source + TLA-resident-CD8 feed-forward)');
% CXCL13 has a short tissue diffusion length L = sqrt(D/k_deg) ~ 0.1-0.5 mm
% (10 kDa, GAG-binding) — two-plus orders below the tumor radius (~16 mm) — so it
% is NOT well-mixed in V_T: it concentrates over ~niche scale. The TLA-resident
% feed-forward source (k_CXCL13_ff*CD8_TLA) is therefore diluted in the LOCAL
% niche volume V_niche = TLA*V_T (TLA is the niche volume fraction), with a
% diffusion-length-pocket floor (+1e-6, the nucleation seed; avoids 0/0 at TLA=0
% and adds no new parameter). The GVAX-depot source (mature DCs in V_ID) is a
% recruitment signal, not niche-confined, so it stays diluted in the full V_T.
% NB: at QSS CD8_TLA ~ TLA, so the feed-forward term's TLA cancels — the local
% CXCL13 is a stable level set by bulk CD8 + trafficking, not a bistable switch.
addrule(model, 'CXCL13 = k_CXCL13_sec*V_ID.mAPC/(k_CXCL13_deg*V_T) + k_CXCL13_ff*V_T.CD8_TLA/(k_CXCL13_deg*(TLA + 1e-6)*V_T)', 'repeatedAssignment');

% CXCL13 saturation for niche formation.
p = addparameter(model, 'H_CXCL13', 0, 'ValueUnits', 'dimensionless', ...
    'ConstantValue', false, 'Notes', 'Hill of CXCL13 driving TLA niche formation');
addrule(model, 'H_CXCL13 = CXCL13/(CXCL13 + CXCL13_50_TLA)', 'repeatedAssignment');

% TLA niche fraction (0-1) — an ABSTRACT scalar (not a physical cell/mass pool),
% so it is a dimensionless parameter with a rate rule (an ODE) rather than a
% SimBiology species. Logistic formation toward a CXCL13-set ceiling, slow decay.
p = addparameter(model, 'TLA', 0, 'ValueUnits', 'dimensionless', ...
    'ConstantValue', false, 'Notes', 'TLA niche fraction (0-1); intratumoral tertiary lymphoid aggregate presence/size');
addrule(model, 'TLA = k_TLA_form*H_CXCL13*(1 - TLA) - k_TLA_decay*TLA', 'rate');

%% ---- Total intratumoral CD8 (mass-balance coupling to bulk dynamics) ---------
% The TLA sub-pools (CD8_TLA, CD8_TLA_act) are RELABELED bulk CD8, not new cells.
% But Teff_module's bulk-CD8 logistic source self-limits on V_T.CD8 ALONE, so
% sequestering cells into the niche would relax the carrying-capacity headroom
% and let the bulk source over-refill — inflating total intratumoral CD8 (the
% opposite of the flat-bulk-CD8 constraint). Define the true total and point the
% carrying-capacity DENOMINATOR at it so the niche cells still count against
% capacity. The proliferation PREFACTOR stays V_T.CD8 (only bulk proliferates;
% niche residence is a relabel, not a division source — decision 4).
p = addparameter(model, 'CD8_T_tot', 0, 'ValueUnits', 'cell', ...
    'ConstantValue', false, 'Notes', 'Total intratumoral CD8 (bulk + within-TLA resting + within-TLA activated + egressed vicinity effector); carrying-capacity denominator');
addrule(model, 'CD8_T_tot = V_T.CD8 + V_T.CD8_TLA + V_T.CD8_TLA_act + V_T.CD8_vic', 'repeatedAssignment');

% Rewrite the bulk-CD8 logistic self-limitation to see CD8_T_tot. Modify the
% existing Teff_module reaction in place (it was added earlier in assembly), the
% same assembly-time-modification pattern macrophage_module uses. Targets only
% the denominator substring, leaving the proliferating-population prefactor.
rxnCD8pro = sbioselect(model, 'Type', 'reaction', 'Where', 'ReactionRate', ...
    'regexp', 'k_CD8_T_pro\*V_T\.CD8');
if isempty(rxnCD8pro)
    warning('tla_module:noCD8pro', ...
        'Bulk-CD8 logistic source (k_CD8_T_pro*V_T.CD8...) not found; CD8_T_tot carrying-capacity coupling not applied.');
else
    for i = 1:numel(rxnCD8pro)
        rate = get(rxnCD8pro(i), 'ReactionRate');
        % Point the density brake at the TOTAL CD8. The proliferation PREFACTOR
        % stays bulk V_T.CD8 (only bulk divides), so niche/vicinity sub-pools count
        % against capacity without themselves dividing.
        %
        % NOTE (2026-07-13): the comment that used to live here argued that letting
        % the logistic brake "go negative" past capacity was harmless -- the rate is
        % proportional to bulk, so bulk would just decay toward zero without crossing
        % it. That is true of the TRAJECTORY and false of the BIOLOGY: a negative
        % brake makes the proliferation reaction a quadratic SINK whose strength
        % scales with k_CD8_T_pro, so raising CD8 proliferation DESTROYED CD8, and
        % GVAX routinely pushed CD8_T_tot past capacity. The brake is now saturating
        % (Teff_module: k_pro*CD8/(1 + CD8_T_tot/cap)) and can never go negative, so
        % this concern is gone -- and it is still smooth for CVODE, which was the
        % other reason max() was avoided. See notes/model/tcell_carrying_capacity_inversion.md.
        % Death is handled by separate reactions.
        newrate = strrep(rate, 'V_T.CD8/(rho_CD8_max*V_T + cell)', ...
                               'CD8_T_tot/(rho_CD8_max*V_T + cell)');
        if ~strcmp(rate, newrate)
            set(rxnCD8pro(i), 'ReactionRate', newrate);
        end
    end
end

%% ---- Vaccine-driven in-situ intratumoral CD8 proliferation (cd8_fc) ----------
% GVAX establishes intratumoral TLAs with Ki-67+ proliferating CD8 (Lutz 2014,
% Cancer Immunol Res 2:616 -- the Zheng group, source of the cd8_fc data).
% Vaccine-primed, antigen-experienced CD8 undergo antigen-driven in-situ
% proliferation on tumoral antigen re-encounter -- the treatment-DIFFERENTIAL half
% of tumor CD8 expansion that the constitutive k_CD8_T_pro (vaccine-blind,
% Blessin-anchored, pins baseline cd8_density) cannot supply.
%
% A SEPARATE saturating source on BULK V_T.CD8 (not the CD8_TLA sub-pool): the
% cd8_fc observable is CD8_total_T = V_T.CD8 + exh + pex + vic (the non-LA bulk;
% CD8_TLA/act are excluded), and the within-TLA pool is only ~a_TLA (~0.55%) of
% bulk -- far too small a reservoir to source a bulk-scale fold at any Ki-67-
% anchored division rate. So the aggregate is the histological SIGNATURE of the
% vaccine proliferative program, but numerically the non-LA expansion is bulk cells
% dividing. Gated by the GVAX antigen Hill H_P1_spread (defined in gvax_module;
% identically 0 at diagnosis and in every no-vaccine arm), which decouples cd8_fc
% from the baseline density. Shares the SAME CD8_T_tot density brake as the
% constitutive term so the two jointly saturate at carrying capacity and cannot run
% away. Zero without GVAX -> baseline / Arm-nivo / no-treatment fits untouched.
% See notes/calibration/cd8_fc_vaccine_expansion_deficit.md.
r = addreaction(model, 'null -> V_T.CD8');
    set(r, 'ReactionRate', 'k_CD8_T_vac_pro*V_T.CD8*H_P1_spread/(1 + CD8_T_tot/(rho_CD8_max*V_T + cell))');
    set(r, 'Notes', 'GVAX-vaccine-gated in-situ tumor CD8 proliferation increment (Lutz 2014 Ki-67+ CD8 within intratumoral TLAs); H_P1_spread-gated (0 without vaccine), saturating in CD8_T_tot density. The treatment-differential half of cd8_fc, decoupled from the constitutive k_CD8_T_pro so baseline cd8_density is untouched. Acts on bulk V_T.CD8 (the non-LA countable pool), not the small CD8_TLA sub-pool -- see header note.');

%% ---- Vicinity effector tumor killing (throttle-bypassing, durable) -----------
% The egressed vicinity pool CD8_vic (Heumann Fig 5d 2-4% outside-TLA activated
% subset) is the urelumab antitumor/survival effector. Two facts force a SEPARATE
% killing pathway rather than adding CD8_vic to Teff_module's bulk k_C_Tcell_eff:
%   (1) Quality, not quantity: 4-1BB upgrades effector function (bulk GZMB density
%       is flat across arms; the survival-associated CD8+CD137+GZMB+ subset is the
%       one that egresses to the neoplastic vicinity).
%   (2) The vicinity pool is intrinsically TINY (~0.3% of bulk CD8; the TLA niche
%       is a_TLA~0.55% of tumor area). Added to the bulk-THROTTLED density it would
%       move killing <1%. But bulk CD8 killing in PDAC is throttled to ~16% by ECM
%       pore exclusion (H_pore_T~0.84 at d21) and further by Treg. Vicinity
%       effectors are POSITIONED at the neoplastic interface and 4-1BB-activated,
%       so they BYPASS pore exclusion and Treg suppression — each kills several-
%       fold more per cell than a bulk CD8, which is how a 0.3% pool matters.
% k_C_vic_eff mirrors k_C_Tcell_eff but DROPS (1-H_pore_T) and the Treg-saturation
% factor, keeps PD-1/TGFb relief, and carries the phi_vic_kill quality multiplier.
% Durability (slow k_vic_decay, reaction 10b) keeps CD8_vic alive for months, so as
% the vaccine/nivo-driven bulk CD8 contracts post-treatment the vicinity pool's
% RELATIVE killing weight grows — the long-horizon recurrence control behind the
% urelumab DFS/OS benefit. Urelumab-gated (CD8_vic~0 without urelumab) so Arms A/B
% and no-treatment are untouched.
p = addparameter(model, 'k_C_vic_eff', 0, 'ValueUnits', '1/day', ...
    'ConstantValue', false, 'Notes', 'Effective vicinity-effector (CD8_vic) tumor-killing rate; throttle-bypassing (no pore exclusion / Treg) with phi_vic_kill quality multiplier');
addrule(model, ['k_C_vic_eff = p_T_kill_per_contact * v_T_search_volume * phi_vic_kill ' ...
    '* V_T.CD8_vic / V_T / cell * (1-H_PD1_C1) * (1-Emax_TGFb_Teff*H_TGFb_Teff)'], 'repeatedAssignment');

%% ---- Reactions ---------------------------------------------------------------
% (1,2,3) CXCL13 source/feed-forward/clearance are folded into the quasi-steady
%         CXCL13 repeatedAssignment above, not reactions.
% (4,5) TLA niche formation/decay are the TLA rate rule above (parameter ODE),
%       not reactions — TLA is an abstract dimensionless niche fraction.
% (6) CD8 residence into the niche (scaled by the TLA AREA fraction a_TLA; relabels existing cells).
r = addreaction(model, 'V_T.CD8 -> V_T.CD8_TLA');
    set(r, 'ReactionRate', 'k_TLA_seq*a_TLA*V_T.CD8');
    set(r, 'Notes', 'Intratumoral CD8 sequestration into the TLA niche (scaled by the TLA AREA fraction a_TLA ~0.0055, Hiraoka 2015 — NOT the 0-1 niche-formation dial TLA, which saturates ~0.82 and over-sequestered CD8 ~150x, draining the non-LA bulk pool). Count partition matches the tla_seq_release_partition submodel target.');
% (7) CD8 return from the niche to the bulk pool (residence turnover, RESTING).
r = addreaction(model, 'V_T.CD8_TLA -> V_T.CD8');
    set(r, 'ReactionRate', 'k_TLA_release*V_T.CD8_TLA');
    set(r, 'Notes', 'Within-TLA resting CD8 return to the bulk intratumoral pool (residence turnover)');
% (7b) ACTIVATED cells undergo the SAME residence turnover (state-independent exit).
% Without this, activated cells leave the niche only via the slow k_TLA_egress,
% so converting resting->activated lengthens niche residence and the TOTAL within-
% TLA number becomes coupled to the activation fraction (Arm C/B number ratio drifts
% to ~2.3 in the model vs Heumann's flat ~1.0 — the cd8_intla_number_invariance
% conflict). Adding release to the activated pool makes the residence exit rate
% state-independent: T_total = S/(k_TLA_release + phi*k_TLA_egress), decoupling the
% number from the CD137+ activation fraction. Egress (10) stays the small EXTRA
% activation-specific flux that builds the vicinity effector pool.
r = addreaction(model, 'V_T.CD8_TLA_act -> V_T.CD8');
    set(r, 'ReactionRate', 'k_TLA_release*V_T.CD8_TLA_act');
    set(r, 'Notes', 'Activated within-TLA CD8 residence turnover to bulk (same k_TLA_release as resting; decouples within-TLA NUMBER from the CD137+ activation fraction so Arm C/B number ratio stays ~1, Heumann Suppl Fig 9)');
% (8) CD137 activation-state conversion (the urelumab effect), with baseline floor.
r = addreaction(model, 'V_T.CD8_TLA -> V_T.CD8_TLA_act');
    set(r, 'ReactionRate', 'k_CD137_act*(f_act_base + (1 - f_act_base)*L_CD137)*V_T.CD8_TLA');
    set(r, 'Notes', 'CD8_TLA -> CD137+ activated state; driven by the durable exposure LATCH L_CD137 (not instantaneous H_CD137_T) so the CD137+ fraction persists to d21 after single-dose urelumab washes out (Heumann Arm C ~83%); f_act_base floor keeps Arms A/B low');
% (9) Reversion to the resting state.
r = addreaction(model, 'V_T.CD8_TLA_act -> V_T.CD8_TLA');
    set(r, 'ReactionRate', 'k_CD137_deact*V_T.CD8_TLA_act');
    set(r, 'Notes', 'Loss of CD137+ activation state (sets the ~75% asymptote with k_CD137_act)');
% (10) Egress of activated cells to the neoplastic vicinity as the durable, higher-
% quality effector pool CD8_vic (the 2-4% CD137+ subset outside the TLA, Fig 5d).
% This is the EXTRA activation-specific exit on top of residence turnover (7b);
% it is the mechanistic route by which urelumab exports effector QUALITY to the
% tumor (vs Arms A/B where CD137+ cells stay locked inside the TLA).
r = addreaction(model, 'V_T.CD8_TLA_act -> V_T.CD8_vic');
    set(r, 'ReactionRate', 'k_TLA_egress*V_T.CD8_TLA_act');
    set(r, 'Notes', 'Activated within-TLA CD8 egress to the neoplastic vicinity effector pool CD8_vic (Fig 5d 2-4% subset; the urelumab antitumor effector export)');
% (10b) Vicinity effector durability: slow decay of CD8_vic back to bulk CD8
% (effector-memory-like persistence). Sets the steady vicinity fraction
% (egress flux / k_vic_decay) and the long-horizon persistence of the killing
% advantage that underlies the urelumab DFS/OS benefit.
r = addreaction(model, 'V_T.CD8_vic -> V_T.CD8');
    set(r, 'ReactionRate', 'k_vic_decay*V_T.CD8_vic');
    set(r, 'Notes', 'Vicinity effector (CD8_vic) turnover back to bulk CD8 (effector-memory durability; slow)');
% (13) Vicinity-effector tumor killing (throttle-bypassing; the urelumab antitumor
% pathway). Separate from bulk k_C_Tcell_eff killing so the two are additive, not
% double-counted. Mirrors Teff_module's killing-reaction + R_Tcell antigen-release
% wiring for the single PDAC cancer clone C1.
cancer_vic = {'C1'};
for i = 1:numel(cancer_vic)
    ci = cancer_vic{i};
    if isempty(sbioselect(model, 'Name', ci, 'Type', 'species'))
        warning('tla_module:noCancer', 'Cancer species %s not found; vicinity killing reaction skipped.', ci);
        continue;
    end
    r = addreaction(model, ['V_T.' ci ' -> V_T.C_x']);
        set(r, 'ReactionRate', ['k_C_vic_eff*V_T.' ci]);
        set(r, 'Notes', 'Cancer cell killing by egressed vicinity effectors (CD8_vic; throttle-bypassing, urelumab antitumor pathway)');
    % Append to the model-wide cancer-kill accumulator R_Tcell so vicinity killing
    % feeds antigen release consistently (same convention as Teff/macrophage).
    rR = sbioselect(model, 'Type', 'rule', 'Where', 'Rule', 'regexp', ['R_Tcell\s*=']);
    if isempty(rR)
        warning('tla_module:noRTcell', 'R_Tcell rule not found; vicinity killing not added to antigen release.');
    else
        for j = 1:numel(rR)
            rule = get(rR(j), 'Rule');
            set(rR(j), 'Rule', [rule '+k_C_vic_eff*V_T.' ci]);
        end
    end
end
% (14) Urelumab/4-1BB drive on the progenitor-exhausted (Tpex) reservoir.
% CD137 agonism accelerates Tpex DIFFERENTIATION out of the reservoir, gated by
% tumor urelumab occupancy H_CD137_T. Biologically this is NOT a clean reflux to
% resting functional effectors: 4-1BB agonism drives Tpex proliferation AND
% TERMINAL differentiation (Pichler 2023, Immunity, PMID via Zotero: "TCR-
% independent CD137 (4-1BB) signaling promotes CD8+-exhausted T cell proliferation
% and terminal differentiation"). So the CD137-driven Tpex exit flux is FATE-SPLIT:
% a functional/transitory-effector fraction f_CD137_func goes to bulk CD8, and the
% remaining (1-f_CD137_func) goes to terminal CD8_exh. f_CD137_func is skewed
% terminal (default ~0.35) and is anchorable from the functional:terminal ratio in
% head-to-head 4-1BB Tpex assays (Kim 2020 HCC / Pichler 2023 / Muth 2023). This
% replaces the prior 100%-functional "reflux" (which excluded the terminal fate
% and mislabeled this as the dominant antitumor lever — a 2026-07-15 trace on the
% Heumann arm showed it moves bulk CD8 <0.5% over 21 days; the urelumab antitumor
% effect rides on the vicinity/quality arm, not this one). k_CD137_reinvig is now
% the TOTAL CD137-driven Tpex differentiation rate. Synergizes with the anti-PD-1
% k_reinvig arm on the SAME reservoir. Inert in Arms A/B and no-treatment
% (H_CD137_T=0 without urelumab).
if isempty(sbioselect(model, 'Name', 'CD8_pex', 'Type', 'species'))
    warning('tla_module:noCD8pex', 'CD8_pex (Tpex) not found; urelumab reinvigoration reaction skipped.');
else
    % (14a) CD137-driven Tpex -> functional (transitory-effector) fate.
    r = addreaction(model, 'V_T.CD8_pex -> V_T.CD8');
        set(r, 'ReactionRate', 'k_CD137_reinvig*V_T.CD8_pex*H_CD137_T*f_CD137_func');
        set(r, 'Notes', 'Urelumab/4-1BB reinvigoration, FUNCTIONAL fate: CD137-driven Tpex differentiation to functional (transitory-effector) CD8, fraction f_CD137_func of the total CD137 Tpex exit flux (gated by H_CD137_T; synergizes with anti-PD-1 k_reinvig)');
    % (14b) CD137-driven Tpex -> TERMINAL fate (Pichler 2023: CD137 drives terminal diff).
    r = addreaction(model, 'V_T.CD8_pex -> V_T.CD8_exh');
        set(r, 'ReactionRate', 'k_CD137_reinvig*V_T.CD8_pex*H_CD137_T*(1 - f_CD137_func)');
        set(r, 'Notes', 'Urelumab/4-1BB, TERMINAL fate: CD137-driven Tpex terminal differentiation to CD8_exh, fraction (1-f_CD137_func) of the total CD137 Tpex exit flux (Pichler 2023 — the biologically-correct cost arm the prior pure-functional reflux omitted)');
    % (14c) CD137-driven Tpex PROLIFERATION (self-renewal). 4-1BB agonism expands the
    % Tpex pool (Pichler 2023 / Li 2023 Ki67), not just re-partitions its exit. This is
    % a GENERATIVE null->CD8_pex source — the SAME topology retired at baseline on
    % 2026-06-17 (it self-ignited the reservoir to the niche cap, decoupling terminal
    % exhaustion from the exhaustion routes) — but here GATED by H_CD137_T, so it is
    % identically zero at baseline / Arms A/B (the regime the retirement protected) and
    % active only transiently under urelumab, where CD137-driven Tpex proliferation IS
    % the intended biology. Niche-capped by the SAME logistic as the exhaustion inflow
    % (1 - CD8_pex/(rho_pex_max*V_T)), so proliferation + inflow jointly saturate at the
    % Jansen 2019 reservoir ceiling and cannot run away. k_CD137_prolif is anchorable
    % from CD137-induced Tpex Ki67 / fold-expansion (Li 2023).
    r = addreaction(model, 'null -> V_T.CD8_pex');
        set(r, 'ReactionRate', 'k_CD137_prolif*V_T.CD8_pex*H_CD137_T*(1 - V_T.CD8_pex/(rho_pex_max*V_T + cell))');
        set(r, 'Notes', 'Urelumab/4-1BB CD137-driven Tpex proliferation (self-renewal); generative but H_CD137_T-gated (zero at baseline/Arms A/B) and niche-capped by rho_pex_max (Pichler 2023 proliferation; anchorable from Li 2023 Ki67)');
end

% (11,12) NOTE (2026-06-18, CD8 Proposal 1): the within-TLA CD8 deaths
%   V_T.CD8_TLA -> null     : k_CD8_death*V_T.CD8_TLA
%   V_T.CD8_TLA_act -> null : k_CD8_death*V_T.CD8_TLA_act
% were REMOVED alongside the bulk tumoral CD8 death in Teff_module.m. These pools
% are relabeled bulk CD8 (sequestered/released, not new cells), and every one of
% them routes back to bulk -- CD8_TLA and CD8_TLA_act via k_TLA_release, CD8_vic
% via k_vic_decay -- so they are cleared by the bulk k_CD8_death rather than
% carrying their own. That remains correct and is unaffected by the 2026-07-13
% reformulation.
%
% STALE CLAIM CORRECTED (2026-07-13): this note used to say "tumoral CD8 turnover
% is now entirely in the logistic" and that removing the sub-pool deaths let
% CD8_T_tot "cleanly pin to rho_CD8_max". Neither holds. The bulk k_CD8_death was
% RESTORED on 2026-06-22, and the renewal term is now saturating and strictly
% positive, so it contributes NO turnover at all -- the pool is bounded by
% k_CD8_death alone, and rho_CD8_max is a half-saturation density, not a pin.
% "Pinning to the cap" was the defect, not the design: see
% notes/model/tcell_carrying_capacity_inversion.md. Circulating CD8 death
% (Teff_module V_C/V_P/V_LN) is unchanged.

end

function add_param_if_absent(model, params, nm)
% Add a model parameter from the params struct only if it is not already present.
existing = sbioselect(model, 'Name', nm, 'Type', 'parameter');
if ~isempty(existing); return; end
p = addparameter(model, nm, params.(nm).Value, 'ValueUnits', params.(nm).Units, ...
    'Notes', params.(nm).Notes);
p.UserData.Source = params.(nm).Notes;
if isfield(params.(nm), 'derivedFrom')
    p.UserData.derivedFrom = params.(nm).derivedFrom;
end
end
