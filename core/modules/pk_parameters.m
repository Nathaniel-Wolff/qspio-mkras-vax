% PK parameter database
%
% Inputs: name -- string containing drug name <optional>:
%                 - nivolumab
%                 - pembrolizumab
%                 - urelumab
%
% Output: params_out -- object containing model parameters
%                       - q_P--rate of diffusive transport C<->P
%                       - q_T--rate of diffusive transport C<->T
%                       - q_LN--rate of diffusive transport C<->LN
%                       - q_LD--rate of convective transport T->LN
%                       - k_cl--clearence rate from central
%                       - gamma_C--volume fraction in C
%                       - gamma_P--volume fraction in P
%                       - gamma_T--volume fraction in T
%                       - gamma_LN--volume fraction in LN


function params_out = pk_parameters(name)

% Optional Inputs
if ~exist('name','var')
    name = '';
end

% Drug Database
% - Transport Rates (q_P,q_T,q_LN,q_LD)
% - Clearence Rate (k_cl)
% - Volume Fraction (gamma_C,gamma_P,gamma_T,gamma_LN)
switch (name)

    case 'pembrolizumab'
        % Central to Peripheral
        q_P.Value = 9.16e-6;
        q_P.Units = 'liter/second';
        q_P.Notes = '(Ahamadi 2016, PMID: 27863186)';
        % Central to Tumour
        q_T.Value = 8.52e-5;
        q_T.Units = 'milliliter/second';
        q_T.Notes = '(Wang 2019, PMID: 31218069)';
        % Central to LN
        q_LN.Value = 3.25e-6;
        q_LN.Units = 'milliliter/second';
        q_LN.Notes = '(Meijer 2017, PMID: 28510992)';
        % Tumour to LN
        q_LD.Value = 0.0015;
        q_LD.Units = '1/minute';
        q_LD.Notes = '(Zhu 1996, PMID: 8706023)';
        % Clearence
        k_cl.Value = 0.315;
        k_cl.Units = 'liter/day';
        k_cl.Notes = '(Ahamadi 2016, PMID: 27863186)';
        % Volume Fractions
        gamma_C.Value = 0.698;
        gamma_C.Notes = '(Ahamadi 2016, PMID: 27863186)';
        gamma_P.Value = 0.0789;
        gamma_P.Notes = '(Ahamadi 2016, PMID: 27863186)';
        gamma_T.Value = 0.522;
        gamma_T.Notes = '(Coughlin 2010, PMID: 20125040; Finley 2015, PMID: 26783500)';
        gamma_LN.Value = 0.2;
        gamma_LN.Notes = '(Wang 2019, PMID: 31218069)';

    case 'nivolumab'
        % Central to Peripheral
        q_P.Value = 8.8e-6;
        q_P.Units = 'liter/second';
        q_P.Notes = '(Bajaj 2017, PMID: 28019091)';
        % Central to Tumour
        q_T.Value = 8.52e-5;
        q_T.Units = 'milliliter/second';
        q_T.Notes = '(Wang 2019, PMID: 31218069)';
        % Central to LN
        q_LN.Value = 3.25e-6;
        q_LN.Units = 'milliliter/second';
        q_LN.Notes = '(Meijer 2017, PMID: 28510992)';
        % Tumour to LN
        q_LD.Value = 0.0015;
        q_LD.Units = '1/minute';
        q_LD.Notes = '(Zhu 1996, PMID: 8706023)';
        % Clearence
        k_cl.Value = 0.37;
        k_cl.Units = 'liter/day';
        k_cl.Notes = '(Bajaj 2017, PMID: 28019091)';
        % Volume Fractions
        gamma_C.Value = 0.61;
        gamma_C.Notes = '(Bajaj 2017, PMID: 28019091)';
        gamma_P.Value = 0.0452;
        gamma_P.Notes = '(Bajaj 2017, PMID: 28019091)';
        gamma_T.Value = 0.522;
        gamma_T.Notes = '(Coughlin 2010, PMID: 20125040; Finley 2015, PMID: 26783500)';
        gamma_LN.Value = 0.2;
        gamma_LN.Notes = '(Wang 2019, PMID: 31218069)';

    case 'urelumab'
        % Urelumab (BMS-663513): fully human IgG4 anti-CD137/4-1BB agonist mAb.
        % Tissue distribution (q_T, q_LN, q_LD, gamma_T, gamma_LN) uses the
        % framework-generic IgG values shared by ALL mAbs here (target-agnostic
        % antibody tissue penetration). Clearance (k_cl) is NOT borrowed from
        % nivolumab — it is set to reproduce urelumab's MEASURED serum
        % half-life of ~5.2-5.6 days at the relevant 0.1-0.3 mg/kg dose
        % (Segal 2017 integrated phase-1 analysis, PMID 27756788; the trial dose
        % is 8 mg flat ~ 0.114 mg/kg). The nivolumab-class k_cl=0.36 L/day gives
        % a model serum t_half of ~19 d (nivolumab-like, ~3.6x too long for
        % urelumab); k_cl=1.7 L/day reproduces t_half~5.4 d in-model
        % (scripts/trace_urelumab_pk.m). Central q_P / gamma_C / gamma_P retain
        % IgG-class values (distribution phase; no urelumab-specific Vss data).
        % Linear clearance is adequate here: phase-1 Cmax/AUC are dose-
        % proportional at 0.1-0.3 mg/kg and the CD137 target sink (receptors on
        % a few hundred activated T cells) is negligible vs the ~33 nmol dose,
        % so TMDD is not modeled (see notes/model/pdac_changelog.md Phase 17).
        % Central to Peripheral
        q_P.Value = 8.8e-6;
        q_P.Units = 'liter/second';
        q_P.Notes = '(IgG-class distribution; Bajaj 2017, PMID: 28019091)';
        % Central to Tumour
        q_T.Value = 8.52e-5;
        q_T.Units = 'milliliter/second';
        q_T.Notes = '(Wang 2019, PMID: 31218069)';
        % Central to LN
        q_LN.Value = 3.25e-6;
        q_LN.Units = 'milliliter/second';
        q_LN.Notes = '(Meijer 2017, PMID: 28510992)';
        % Tumour to LN
        q_LD.Value = 0.0015;
        q_LD.Units = '1/minute';
        q_LD.Notes = '(Zhu 1996, PMID: 8706023)';
        % Clearence — set to match urelumab measured serum t_half ~5.4 d
        k_cl.Value = 1.7;
        k_cl.Units = 'liter/day';
        k_cl.Notes = '(fit to urelumab serum t_half ~5.2-5.6 d at 0.1-0.3 mg/kg, Segal 2017, PMID: 27756788; in-model t_half~5.4 d via scripts/trace_urelumab_pk.m)';
        % Volume Fractions
        gamma_C.Value = 0.61;
        gamma_C.Notes = '(IgG-class distribution; Bajaj 2017, PMID: 28019091)';
        gamma_P.Value = 0.0452;
        gamma_P.Notes = '(IgG-class distribution; Bajaj 2017, PMID: 28019091)';
        gamma_T.Value = 0.522;
        gamma_T.Notes = '(Coughlin 2010, PMID: 20125040; Finley 2015, PMID: 26783500)';
        gamma_LN.Value = 0.2;
        gamma_LN.Notes = '(Wang 2019, PMID: 31218069)';

end

% Diffusive Transport: C<->P
params_out.q_P = q_P;
% Diffusive Transport: C<->T
params_out.q_T = q_T;
% Diffusive Transport: C<->LN
params_out.q_LN = q_LN;
% Convective Transport: T->LN || LN->C
params_out.q_LD = q_LD;
% Clearence from Central
params_out.k_cl = k_cl;
% Volume Fractions

gamma_C.Units = 'dimensionless';
gamma_P.Units = 'dimensionless';
gamma_T.Units = 'dimensionless';
gamma_LN.Units = 'dimensionless';
params_out.gamma_C = gamma_C;
params_out.gamma_P = gamma_P;
params_out.gamma_T = gamma_T;
params_out.gamma_LN = gamma_LN;

% Add derivedFrom field to all parameters for provenance tracking
names = fieldnames(params_out);
for i = 1:numel(names)
    params_out.(names{i}).derivedFrom = {names{i}};
end

end
