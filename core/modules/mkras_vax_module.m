function model = mkras_vax_module(model, params, varagin)
%modeled after clinical trial in resected PDAC, https://pmc.ncbi.nlm.nih.gov/articles/PMC12891733/
%varagin might be used later

%mKRAS-VAX Module with subcutaneous injection site compartment (V_SC)
% -All 6 peptides (G12V, G12A, G12R, G12C, G12D, G13D) administered at once to V_SC 
% -No modeling of Poly-ICLC

% ---- Create subcutaneous injection site compartment and species ----
%Creating subcutaneous compartment
V_SC = sbioselect(model, "Name", 'V_SC', 'Type', 'compartment');
if isempty(V_SC)
    V_SC = addcompartment(model, 'V_SC', 'CapacityUnits', 'liter');
    set(V_SC, 'Notes', 'Subcutaneous injection site compartment');
end

%Adding each peptide species separately in injection site(though from dosing, concentration of each will have to calculated)
G12D = sbioselect(model, 'Name', 'G12D', 'Type', 'species');
if isempty(G12D)
    G12D = addspecies(V_SC, 'G12D', 0, 'InitialAmountUnits', 'micromoles');
    set(G12D, 'Notes', 'G12D concentration at subcutaneous injection site');
end

G12V = sbioselect(model, 'Name', 'G12V', 'Type', 'species');
if isempty(G12V)
    G12V = addspecies(V_SC, 'G12V', 0, 'InitialAmountUnits', 'micromoles');
    set(G12V, 'Notes', 'G12V concentration at subcutaneous injection site');
end

G12R = sbioselect(model, 'Name', 'G12R', 'Type', 'species');
if isempty(G12R)
    G12R = addspecies(V_SC, 'G12R', 0, 'InitialAmountUnits', 'micromoles');
    set(G12R, 'Notes', 'G12R concentration at subcutaneous injection site');
end


G12C = sbioselect(model, 'Name', 'G12C', 'Type', 'species');
if isempty(G12C)
    G12C = addspecies(V_SC, 'G12C', 0, 'InitialAmountUnits', 'micromoles');
    set(G12C, 'Notes', 'G12C concentration at subcutaneous injection site');
end

G12A = sbioselect(model, 'Name', 'G12A', 'Type', 'species');
if isempty(G12A)
    G12A = addspecies(V_SC, 'G12A', 0, 'InitialAmountUnits', 'micromoles');
    set(G12A, 'Notes', 'G12A concentration at subcutaneous injection site');
end

G13D = sbioselect(model, 'Name', 'G13D', 'Type', 'species');
if isempty(G13D)
    G13D = addspecies(V_SC, 'G13D', 0, 'InitialAmountUnits', 'micromoles');
    set(G13D, 'Notes', 'G13D concentration at subcutaneous injection site');
end


% Add APC species at injection site
apc_sc = sbioselect(model, 'Name', 'APC', 'Type', 'species', 'Parent', V_SC);
if isempty(apc_sc)
    apc_SC = addspecies(V_SC, 'APC', 0, 'InitialAmountUnits', 'cell');
    set(apc_SC, 'Notes', 'Naive APCs recruited to subcutaneous injection site');
end

% Add mature APC species at injection site
mapc_sc = sbioselect(model, 'Name', 'mAPC', 'Type', 'species', 'Parent', V_SC);
if isempty(mapc_sc)
    mapc_sc = addspecies(V_SC, 'mAPC', 0, 'InitialAmountUnits', 'cell');
    set(mapc_sc, 'Notes', 'Mature APCs with antigen at subcutaneous injection site');
end


% Load parameters from params structure

% Peptide transport from injection site into axillary lymph node
% These may be deprecated later

k_transport_mKRAS = sbioselect(model,'Type','parameter','Name','k_transport_mKRAS');
if isempty(k_transport_mKRAS), k_transport_mKRAS = addparameter(model,'k_transport_mKRAS',params.k_transport_mKRAS.Value,'ValueUnits',params.k_transport_mKRAS.Units); end
set(k_transport_mKRAS, 'Notes', params.k_transport_mKRAS.Notes);

k_transport_G12D = sbioselect(model,'Type','parameter','Name','k_transport_G12D');
if isempty(k_transport_G12D), k_transport_G12D = addparameter(model,'k_transport_G12D',params.k_transport_G12D.Value,'ValueUnits',params.k_transport_G12D.Units); end
set(k_transport_G12D, 'Notes', params.k_transport_G12D.Notes);

k_transport_G12V = sbioselect(model,'Type','parameter','Name','k_transport_G12V');
if isempty(k_transport_G12V), k_transport_G12V = addparameter(model,'k_transport_G12V',params.k_transport_G12V.Value,'ValueUnits',params.k_transport_G12V.Units); end
set(k_transport_G12V, 'Notes', params.k_transport_G12V.Notes);

k_transport_G12R = sbioselect(model,'Type','parameter','Name','k_transport_G12R');
if isempty(k_transport_G12R), k_transport_G12R = addparameter(model,'k_transport_G12R',params.k_transport_G12R.Value,'ValueUnits',params.k_transport_G12R.Units); end
set(k_transport_G12R, 'Notes', params.k_transport_G12R.Notes);

k_transport_G12C = sbioselect(model,'Type','parameter','Name','k_transport_G12C');
if isempty(k_transport_G12C), k_transport_G12C = addparameter(model,'k_transport_G12C',params.k_transport_G12C.Value,'ValueUnits',params.k_transport_G12C.Units); end
set(k_transport_G12C, 'Notes', params.k_transport_G12C.Notes);

k_transport_G12A = sbioselect(model,'Type','parameter','Name','k_transport_G12A');
if isempty(k_transport_G12A), k_transport_G12A = addparameter(model,'k_transport_G12A',params.k_transport_G12A.Value,'ValueUnits',params.k_transport_G12A.Units); end
set(k_transport_G12A, 'Notes', params.k_transport_G12A.Notes);

k_transport_G13D = sbioselect(model,'Type','parameter','Name','k_transport_G13D');
if isempty(k_transport_G13D), k_transport_G13D = addparameter(model,'k_transport_G13D',params.k_transport_G13D.Value,'ValueUnits',params.k_transport_G13D.Units); end
set(k_transport_G13D, 'Notes', params.k_transport_G13D.Notes);

% Injection site volume
V_SC_vol = sbioselect(model,'Type','parameter','Name','V_SC_volume');
if isempty(V_SC_vol), V_SC_vol = addparameter(model,'V_SC_volume',params.V_SC_volume.Value,'ValueUnits',params.V_SC_volume.Units); end
set(V_SC_vol, 'Notes', params.V_SC_volume.Notes);
set(V_SC, 'Capacity', V_SC_vol.Value);


% APC recruitment to injection site (subcutaneous)
k_APC_recruit_SC = sbioselect(model,'Type','parameter','Name','k_APC_recruit_SC');
if isempty(k_APC_recruit_SC), k_APC_recruit_SC = addparameter(model,'k_APC_recruit_SC',params.k_APC_recruit_SC.Value,'ValueUnits',params.k_APC_recruit_SC.Units); end
set(k_APC_recruit_SC, 'Notes', params.k_APC_recruit_SC.Notes);

% APC maturation at injection site (subcutaneous)
k_APC_mature_SC = sbioselect(model,'Type','parameter','Name','k_APC_mature_SC');
if isempty(k_APC_mature_SC), k_APC_mature_SC = addparameter(model,'k_APC_mature_SC',params.k_APC_mature_SC.Value,'ValueUnits',params.k_APC_mature_SC.Units); end
set(k_APC_mature_SC, 'Notes', params.k_APC_mature_SC.Notes);

% APC migration from injection site (subcutaneous) to lymph nodes
k_APC_mig_SC_LN = sbioselect(model,'Type','parameter','Name','k_APC_mig_SC_LN');
if isempty(k_APC_mig_SC_LN), k_APC_mig_SC_LN = addparameter(model,'k_APC_mig_SC_LN',params.k_APC_mig_SC_LN.Value,'ValueUnits',params.k_APC_mig_SC_LN.Units); end
set(k_APC_mig_SC_LN, 'Notes', params.k_APC_mig_SC_LN.Notes);

% APC death at injection site (subcutaneous)
k_APC_death_SC = sbioselect(model,'Type','parameter','Name','k_APC_death_SC');
if isempty(k_APC_death_SC), k_APC_death_SC = addparameter(model,'k_APC_death_SC',params.k_APC_death_SC.Value,'ValueUnits',params.k_APC_death_SC.Units); end
set(k_APC_death_SC, 'Notes', params.k_APC_death_SC.Notes);

% Antigen consumption per APC maturation - assumed to be equal for all mKRAS
sigma_P1_per_mAPC_mKRAS = sbioselect(model,'Type','parameter','Name','sigma_P1_per_mAPC_mKRAS');
if isempty(sigma_P1_per_mAPC_mKRAS), sigma_P1_per_mAPC_mKRAS = addparameter(model,'sigma_P1_per_mAPC_mKRAS',params.sigma_P1_per_mAPC_mKRAS.Value,'ValueUnits',params.sigma_P1_per_mAPC_mKRAS.Units); end
set(sigma_P1_per_mAPC_mKRAS, 'Notes', params.sigma_P1_per_mAPC_mKRAS.Notes);

% EC50 for antigen-dependent APC maturation
% All 6 values for all different mKRAS variants
EC50_G12D_mature = sbioselect(model,'Type','parameter','Name','EC50_G12D_mature');
if isempty(EC50_G12D_mature), EC50__G12D_mature = addparameter(model,'EC50_G12D_mature',params.EC50_G12D_mature.Value,'ValueUnits',params.EC50_G12D_mature.Units); end
set(EC50_G12D_mature, 'Notes', params.EC50_G12D_mature.Notes);

EC50_G12V_mature = sbioselect(model,'Type','parameter','Name','EC50_G12V_mature');
if isempty(EC50_G12V_mature), EC50_G12V_mature = addparameter(model,'EC50_G12V_mature',params.EC50_G12V_mature.Value,'ValueUnits',params.EC50_G12V_mature.Units); end
set(EC50_G12V_mature, 'Notes', params.EC50_G12V_mature.Notes);

EC50_G12R_mature = sbioselect(model,'Type','parameter','Name','EC50_G12R_mature');
if isempty(EC50_G12R_mature), EC50_G12R_mature = addparameter(model,'EC50_G12R_mature',params.EC50_G12R_mature.Value,'ValueUnits',params.EC50_G12R_mature.Units); end
set(EC50_G12R_mature, 'Notes', params.EC50_G12R_mature.Notes);

EC50_G12C_mature = sbioselect(model,'Type','parameter','Name','EC50_G12C_mature');
if isempty(EC50_G12C_mature), EC50_G12C_mature = addparameter(model,'EC50_G12C_mature',params.EC50_G12C_mature.Value,'ValueUnits',params.EC50_G12C_mature.Units); end
set(EC50_G12C_mature, 'Notes', params.EC50_G12C_mature.Notes);

EC50_G12A_mature = sbioselect(model,'Type','parameter','Name','EC50_G12A_mature');
if isempty(EC50_G12A_mature), EC50_G12A_mature = addparameter(model,'EC50_G12A_mature',params.EC50_G12A_mature.Value,'ValueUnits',params.EC50_G12A_mature.Units); end
set(EC50_G12A_mature, 'Notes', params.EC50_G12A_mature.Notes);

EC50_G13D_mature = sbioselect(model,'Type','parameter','Name','EC50_G13D_mature');
if isempty(EC50_G13D_mature), EC50_G13D_mature = addparameter(model,'EC50_G13D_mature',params.EC50_G13D_mature.Value,'ValueUnits',params.EC50_G13D_mature.Units); end
set(EC50_G13D_mature, 'Notes', params.EC50_G13D_mature.Notes);

% Capacity and baseline recruitment to injection site (subcutaneous)
K_APC_SC = sbioselect(model,'Type','parameter','Name','K_APC_SC');
if isempty(K_APC_SC), K_APC_SC = addparameter(model,'K_APC_SC',params.K_APC_SC.Value,'ValueUnits','cell'); end

k_APC_recruit0_SC = sbioselect(model,'Type','parameter','Name','k_APC_recruit0_SC');
if isempty(k_APC_recruit0_SC), k_APC_recruit0_SC = addparameter(model,'k_APC_recruit0_SC',params.k_APC_recruit0_SC.Value,'ValueUnits',params.k_APC_recruit0_SC.Units); end

% Separate mature-DC death
k_mAPC_death_SC = sbioselect(model,'Type','parameter','Name','k_mAPC_death_SC');
if isempty(k_mAPC_death_SC), k_mAPC_death_SC = addparameter(model,'k_mAPC_death_SC',params.k_mAPC_death_SC.Value,'ValueUnits',params.k_mAPC_death_SC.Units); end

% Antigen yield and degradation/drain for each of the mKRAS variants
k_G12D_deg_SC = sbioselect(model,'Type','parameter','Name','k_G12D_deg_SC');
if isempty(k_G12D_deg_SC), k_G12D_deg_SC = addparameter(model,'k_G12D_deg_SC',params.k_G12D_deg_SC.Value,'ValueUnits',params.k_G12D_deg_SC.Units); end

k_G12V_deg_SC = sbioselect(model,'Type','parameter','Name','k_G12V_deg_SC');
if isempty(k_G12V_deg_SC), k_G12V_deg_SC = addparameter(model,'k_G12V_deg_SC',params.k_G12V_deg_SC.Value,'ValueUnits',params.k_G12V_deg_SC.Units); end

k_G12R_deg_SC = sbioselect(model,'Type','parameter','Name','k_G12R_deg_SC');
if isempty(k_G12R_deg_SC), k_G12R_deg_SC = addparameter(model,'k_G12R_deg_SC',params.k_G12R_deg_SC.Value,'ValueUnits',params.k_G12R_deg_SC.Units); end

k_G12C_deg_SC = sbioselect(model,'Type','parameter','Name','k_G12C_deg_SC');
if isempty(k_G12C_deg_SC), k_G12C_deg_SC = addparameter(model,'k_G12C_deg_SC',params.k_G12C_deg_SC.Value,'ValueUnits',params.k_G12C_deg_SC.Units); end

k_G12A_deg_SC = sbioselect(model,'Type','parameter','Name','k_G12A_deg_SC');
if isempty(k_G12A_deg_SC), k_G12A_deg_SC = addparameter(model,'k_G12A_deg_SC',params.k_G12A_deg_SC.Value,'ValueUnits',params.k_G12A_deg_SC.Units); end

k_G13D_deg_SC = sbioselect(model,'Type','parameter','Name','k_G13D_deg_SC');
if isempty(k_G13D_deg_SC), k_G13D_deg_SC = addparameter(model,'k_G13D_deg_SC',params.k_G13D_deg_SC.Value,'ValueUnits',params.k_G13D_deg_SC.Units); end

%degradation rate constants in Lymph Node
k_G12D_deg_ALN = sbioselect(model,'Type','parameter','Name','k_G12D_deg_ALN');
if isempty(k_G12D_deg_ALN), k_G12D_deg_ALN = addparameter(model,'k_G12D_deg_ALN',params.k_G12D_deg_ALN.Value,'ValueUnits',params.k_G12D_deg_ALN.Units); end

k_G12V_deg_ALN = sbioselect(model,'Type','parameter','Name','k_G12V_deg_ALN');
if isempty(k_G12V_deg_ALN), k_G12V_deg_ALN = addparameter(model,'k_G12V_deg_ALN',params.k_G12V_deg_ALN.Value,'ValueUnits',params.k_G12V_deg_ALN.Units); end

k_G12R_deg_ALN = sbioselect(model,'Type','parameter','Name','k_G12R_deg_ALN');
if isempty(k_G12R_deg_ALN), k_G12R_deg_ALN = addparameter(model,'k_G12R_deg_ALN',params.k_G12R_deg_ALN.Value,'ValueUnits',params.k_G12R_deg_ALN.Units); end

k_G12C_deg_ALN = sbioselect(model,'Type','parameter','Name','k_G12C_deg_ALN');
if isempty(k_G12C_deg_ALN), k_G12C_deg_ALN = addparameter(model,'k_G12C_deg_ALN',params.k_G12C_deg_ALN.Value,'ValueUnits',params.k_G12C_deg_ALN.Units); end

k_G12A_deg_ALN = sbioselect(model,'Type','parameter','Name','k_G12A_deg_ALN');
if isempty(k_G12A_deg_ALN), k_G12A_deg_ALN = addparameter(model,'k_G12A_deg_ALN',params.k_G12A_deg_ALN.Value,'ValueUnits',params.k_G12A_deg_ALN.Units); end

k_G13D_deg_ALN = sbioselect(model,'Type','parameter','Name','k_G13D_deg_ALN');
if isempty(k_G13D_deg_ALN), k_G13D_deg_ALN = addparameter(model,'k_G13D_deg_ALN',params.k_G13D_deg_ALN.Value,'ValueUnits',params.k_G13D_deg_ALN.Units); end



% APC migration saturation
K_APC_mig_SC = sbioselect(model,'Type','parameter','Name','K_APC_mig_SC');
if isempty(K_APC_mig_SC), K_APC_mig_SC = addparameter(model,'K_APC_mig_SC',params.K_APC_mig_SC.Value,'ValueUnits',params.K_APC_mig_SC.Units); end

% Antigen-transfer rate: migratory mcDC2 -> resident cDC1 cross-licensing in LN
% Assumed to be the same for mKRAS as GVAX antigen at this time 
k_cDC1_Ag_transfer = sbioselect(model,'Type','parameter','Name','k_cDC1_Ag_transfer');
if isempty(k_cDC1_Ag_transfer), k_cDC1_Ag_transfer = addparameter(model,'k_cDC1_Ag_transfer',params.k_cDC1_Ag_transfer.Value,'ValueUnits',params.k_cDC1_Ag_transfer.Units); end
set(k_cDC1_Ag_transfer, 'Notes', params.k_cDC1_Ag_transfer.Notes);

% (f_mAPC_mcDC1_fate removed in mid-tier #15: the direct depot-moDC->mcDC1
% cross-priming route was retired in favor of the host-cDC1 antigen-handoff
% route -- see the migration block below.)

% THIS SECTION HAS NOT BEEN ADJUSTED YET
% Vaccine-antigen presentation rate: feeds mKRAS antigens into the shared
% endosomal antigen-1 pool (V_e.P1) so the vaccine raises H_P1 (pMHC display),
% not just the DC count. Parallel to tumor-antigen uptake k_up (same 1/day/cell).
k_GVAX_present = sbioselect(model,'Type','parameter','Name','k_GVAX_present');
if isempty(k_GVAX_present), k_GVAX_present = addparameter(model,'k_GVAX_present',params.k_GVAX_present.Value,'ValueUnits',params.k_GVAX_present.Units); end
set(k_GVAX_present, 'Notes', params.k_GVAX_present.Notes);




%creating regional lymph node compartment
V_ALN = sbioselect(model, 'Name', "V_ALN", "Type", "compartment")
if isempty(V_ALN)
    V_ALN = addcompartment(model, 'V_ALN', 'CapacityUnits', 'liter');
    set(V_ALN, 'Notes', 'Regional (Axillary) lymph node compartment');
end

%Adding each peptide species separately at regional lymph node(though from dosing, concentration of each will have to calculated)
G12D_ALN = sbioselect(model, 'Name', 'G12D_ALN', 'Type', 'species');
if isempty(G12D_ALN)
    G12D_ALN = addspecies(V_ALN, 'G12D_ALN', 0, 'InitialAmountUnits', 'micromoles');
    set(G12D_ALN, 'Notes', 'G12D concentration at axillary lymph node.');
end

G12V_ALN = sbioselect(model, 'Name', 'G12V_ALN', 'Type', 'species');
if isempty(G12V_ALN)
    G12V_ALN = addspecies(V_ALN, 'G12V_ALN', 0, 'InitialAmountUnits', 'micromoles');
    set(G12V_ALN, 'Notes', 'G12V concentration at axillary lymph node.');
end

G12R_ALN = sbioselect(model, 'Name', 'G12R_ALN', 'Type', 'species');
if isempty(G12R_ALN)
    G12R_ALN = addspecies(V_ALN, 'G12R_ALN', 0, 'InitialAmountUnits', 'micromoles');
    set(G12R_ALN, 'Notes', 'G12R concentration at axillary lymph node.');
end


G12C_ALN = sbioselect(model, 'Name', 'G12C_ALN', 'Type', 'species');
if isempty(G12C_ALN)
    G12C_ALN = addspecies(V_ALN, 'G12C_ALN', 0, 'InitialAmountUnits', 'micromoles');
    set(G12C_ALN, 'Notes', 'G12C concentration at axillary lymph node.');
end

G12A_ALN = sbioselect(model, 'Name', 'G12A_ALN', 'Type', 'species');
if isempty(G12A_ALN)
    G12A_ALN = addspecies(V_ALN, 'G12A_ALN', 0, 'InitialAmountUnits', 'micromoles');
    set(G12A_ALN, 'Notes', 'G12A concentration at axillary lymph node.');
end

G13D_ALN = sbioselect(model, 'Name', 'G13D_ALN', 'Type', 'species');
if isempty(G13D_ALN)
    G13D_ALN = addspecies(V_ALN, 'G13D_ALN', 0, 'InitialAmountUnits', 'micromoles');
    set(G13D_ALN, 'Notes', 'G13D concentration at axillary lymph node.');
end


%creating blood compartment
V_B = sbioselect(model, 'Name', "V_B", "Type", "compartment")
if isempty(V_B)
    V_B = addcompartment(model, 'V_B', 'CapacityUnits', 'liter');
    set(V_B, 'Notes', 'Blood/circulation compartment');
end 

%creating tumor microenvironment compartment
V_TME = sbioselect(model, 'Name', "V_TME", 'Type', 'compartment');
if isempty(V_TME)
    V_TME = addcompartment(model, 'V_TME', 'CapacityUnits', 'liter');
    set(V_TME, 'Notes', 'Tumor microenvironment compartment');

end
%% --- Subcutaneous Injection Site mKRAS Peptide Dynamics
% Represents the degradation of individual peptides, transport of free
% peptides into lymph nodes

% Peptide degradation and transport from the injection site
% At this time, a single transport rate will be used
peptideNames = {'G12D','G12V','G12R','G12C','G12A','G13D'};
peptideDegradationsSC = {k_G12D_deg_SC, k_G12V_deg_SC, k_G12R_deg_SC, k_G12C_deg_SC, k_G12A_deg_SC, k_G13D_deg_SC}
peptideDegradationsALN = {k_G12D_deg_ALN, k_G12V_deg_ALN, k_G12R_deg_ALN, k_G12C_deg_ALN, k_G12A_deg_ALN, k_G13D_deg_ALN}
%transport_constant = k_transport_mKRAS


for i = 1:numel(peptideNames)
    %form of each peptide transport reaction is G12D_SC -> G12D_ALN
    rxtn_name_sc_deposit = ['null' ' -> ' 'V_SC.' peptideNames{i}]
    rxtn_name_ln_transfer = ['V_SC.' peptideNames{i}  ' -> ' 'V_ALN.' peptideNames{i}]
    rxtn_name_sc_deg = ['V_SC.' peptideNames{i}  ' -> ' 'null']
    rxtn_name_ALN_deg = ['V_ALN.' peptideNames{i}  ' -> ' 'null']

    peptide_transport_ln_rxtn = sbioselect(model, 'Type', 'reaction', 'Name', rxtn_name_ln_transfer);
    if isempty(peptide_transport_ln_rxtn)
        
        peptide = sbioselect(model,'Type','species','Name','peptideNames{i}','Parent',V_SC);
        %deposition of vaccine into SC site: null -> mKRAS_SC
        %discrete jump time deposition is the model here
        input_rxtn = addreaction(model, rxtn_name_sc_deposit)
        mass_action_input = addkineticlaw(input_rxtn, 'MassAction') %placeholder

        %transport in SC: mKRAS_SC -> mKRAS_ALN, transport into ALN 
        % Assuming mass action transport per Yogurtcu et. al 2026
        trans_rxtn = addreaction(model, rxtn_name_ln_transfer); 
        mass_action_transport = addkineticlaw(trans_rxtn, 'MassAction')
        addparameter(mass_action_transport, 'k', k_transport_mKRAS.Value) %k here is equivalent to k_vac,drain in the publication

        %degradation in SC: mKRAS_SC -> null
        deg_SC_rxtn = addreaction(model, rxtn_name_sc_deg);
        mass_action_transport_deg_sc = addkineticlaw(deg_SC_rxtn, 'MassAction')
        addparameter(mass_action_transport_deg_sc, 'k', peptideDegradationsSC{i}.Value) %k_deg_SC here is equivalent to k_vac,deg,depot in the publication


        %degradation in ALN: mKRAS_ALN -> null
        deg_ALN_rxtn = addreaction(model, rxtn_name_ALN_deg);
        mass_action_transport_deg_ALN = addkineticlaw(deg_ALN_rxtn, 'MassAction')
        addparameter(mass_action_transport_deg_ALN, 'k', peptideDegradationsALN{i}.Value) %k_deg_ALN here is equivalent to k_vac,deg,depot in the publication


    end
end
   