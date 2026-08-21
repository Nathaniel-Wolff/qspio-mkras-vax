%% Immune Oncology Model
% Script for setting up and running the immune oncology model in simbiology

%% Add use-defined units into both Simbiology and Symbolic libraries
% Add 'cell' unit to SimBiology and Symbolic Toolboxes
if (isempty(sbioshowunits('cell')))
    cell_unit = sbiounit('cell','molecule');
    sbioaddtolibrary(cell_unit);
end
% Symbolic Unit
u = symunit;
try u.cell;
catch
    newUnit('cell',u.molecule);
end

% Add 'meter_mercury' unit so 'millimeter_mercury' (milli prefix) resolves
if (isempty(sbioshowunits('meter_mercury')))
    mmHg_base = sbiounit('meter_mercury','molecule');
    sbioaddtolibrary(mmHg_base);
end

% Add 'mU' unit to SimBiology and Symbolic Toolboxes
mU_unit = sbiounit('mU','mole/liter');
sbioaddtolibrary(mU_unit);
% Symbolic Unit
u = symunit;
try u.mU;
catch
    newUnit('mU',u.molarity);
end


%% Setup Parameters
% Setup Parameters using hybrid YAML+Legacy approach (includes calculated parameters)
params_in = parameters_PDAC;
params_out = load_parameters(params_in);

%% Create the SimBiology Model
% Model Settings
model_name = 'Immune Oncology Model';
start_time = 0.0; % [days]
time_step = 0.1; % [days] 0.01 days ~ 15 mins
end_time = 365; % [days]
absolute_tolerance = 1e-9;
relative_tolerance = 1e-6;
% solver = 'ode15s';
solver = 'sundials';
% Model Object
time = start_time:time_step:end_time;
model = simbio_init(model_name,time,solver,absolute_tolerance,relative_tolerance,params_out);
%model = simbio_init(model_name,time,params_out);
% Maximal simulation time
config = getconfigset(model);
set(config, 'MaximumWallClock', 120)
% set(config, 'MaximumNumberOfLogs', 1)
set(config.SolverOptions, 'AbsoluteToleranceScaling', false)

%% Add Modules to the Model
% Cancer Modules
% Initial amount is a placeholder; set_baseline_populations sets actual value from density priors
model = cancer_module(model,'C1',params_out,1,'Gompertzian');
% model = cancer_module(model,'C2',params_out,0); % chemotherapy-resistant clone
% T cell Modules
model = Treg_module(model,params_out);
model = CCL5_Treg_module(model,params_out); % CCL5/CCR5-mediated Treg recruitment (Wang 2017 Oncogene + Tan 2009 J Immunol)
model = Teff_module(model,params_out,{'C1'});
% APC Module
model = APC_module(model,params_out);
% Antigen Modules
antigenCP = create_antigen({'C1'},[5.4e-13],'antigenID',0);
model = antigen_module(model,'0',params_out,antigenCP);
antigen   = create_antigen({'C1'},[5.4e-13],'antigenID',1);
model = antigen_module(model,'1',params_out,antigen);
% PK Modules
params_aPD1    = pk_parameters('nivolumab');
params_aPDL1   = pk_parameters('durvalumab');
params_aCTLA4  = pk_parameters('tremelimumab');
params_aCD137  = pk_parameters('urelumab');   % anti-CD137/4-1BB agonist (Heumann Arm C)
model = pk_module(model,'aPD1',params_aPD1);
model = pk_module(model,'aPDL1',params_aPDL1,'n');
model = pk_module(model,'aCTLA4' ,params_aCTLA4);
% Urelumab PK: 4-compartment drug species V_{C,P,T,LN}.aCD137. Always added
% (=0 when undosed) so the TLA-local H_CD137_T CD8-activation Hill in tla_module
% (reads V_T.aCD137) is well-defined for every scenario.
model = pk_module(model,'aCD137',params_aCD137);
% Checkpoint Modules
model = checkpoint_module(model,params_out,'CD8','C1');
%model = checkpoint_module(model,params_out,'CD8','C2');
model = checkpoint_module(model,params_out,'CD8','APC');

% QSPIO-HCC Modules
model = Th_module(model,params_out);
model = cytokine_tcell_module(model, params_out);  % Centralized T-cell cytokine reactions (IL2, IFNg, TGFb)
model = MDSC_module(model,params_out,{'C1'},'inostat',0,'drugName','entinostat');
% model = nabpaclitaxel_module(model,params_out);
% model = cabozantinib_module(model,params_out); % Disabled due to numerical instability
model = macrophage_module(model,params_out,{'C1'},'aCD47',0); % PK module in development for aCD47
% Fibroblast Module
model = fibroblast_module_enhanced(model,params_out);

% Innate immunity + chemokine trafficking (reactivated 2026-05-13)
% NK_module: NK recruitment/killing/IFNg/XCL1 (depends on H_CXCL12_Texcl from fibroblast)
% chemotaxis_module: CXCL9/10 IFNg-induced + CXCR3 recruitment of CD8/Th/NK
% Order: NK first (chemotaxis references V_T.NK / V_P.NK / NK_P_steady),
%        both before gvax_module (baseline immunity before treatment perturbation).
model = NK_module(model, params_out, {'C1'});
model = chemotaxis_module(model, params_out);

model = gvax_module(model, params_out, 'add_cyclophosphamide', true);

% Tertiary lymphoid aggregate (TLA) niche + urelumab/CD137 mechanism.
% After gvax_module (needs V_ID.mAPC), Teff_module (V_T.CD8, k_CD8_death) and
% pk_module('aCD137') (V_T.aCD137). All TLA pools are GVAX-seeded -> zero in the
% baseline / non-vaccinated arms, so the calibrated A/B fit is untouched.
model = tla_module(model, params_out);

% Aggregate observables (derived totals for cal-target denominators).
% Must run after all species-adding modules so nucleated_total_T can
% reference every V_T cell species defined above.
model = aggregate_observables_module(model);

%% Baseline Populations
% Note: set_baseline_populations() is called in initial_conditions() after
% variant parameters are applied. This ensures IC fraction parameters from
% priors/variants are used correctly.
