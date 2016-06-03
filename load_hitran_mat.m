function absorbance_ = load_hitran_mat(hitran_struct, wavenum, ...
		pressure_atm, partial_pressure_atm, ...
		path_length_cm, isotopologues_, molecular_weight_array_amu, T)

%%load_hitran
% hitran_struct = load_hitran(file_name, wavenumber_1_per_cm, ...
%		pressure_atm, partial_pressure_atm, ...
%		path_length_cm, isotopologues_, molecular_weight_array_amu)
%
% Loads 2012 HITRAN format files, and applies a Voigt lineshape to each line 
% within the desired range. Returns data in a structure containing the raw data
% as well as the absorption coefficient spectrum.
%
% Reference for format: HITRAN.org


% Inputs:
%   file_name              - String of HITRAN file path to load
%   wavenumber_1_per_cm    - Array of wavenumber samples on which to 
%                            calculate the lineshape and absorption
%   pressure_atm           - Total pressure of gas.
%   partial_pressure_atm   - Pressure of gas if all gas but the molecular species
%                            of interest were evacuated.
%   path_length_cm         - Path length of the cell.
%   isotopologues_         - Array of isotopologue codes to keep in 
%                            the hitran_struct. Meaning is *loosely* defined in 
%                            L. S. Rothman et al., "AFGL atmospheric absorption 
%                            line parameters compilation: 1982 edition," 
%                            Appl. Opt. 22, 2247 (1983).
%                            As far as I can tell, if the molecule is linear,
%                            you take the ones digit of the nuclear weight of 
%                            each atom in order to get the isotopologue code.
%   molecular_weight_amu   - Array of molecular weights of the isotopologues given 
%                            in isotopologues_.
%
% Outputs:
%   hitran_struct          - A structure containing the following fields:
%     absorption_coefficient_1_per_cm    
%              - An array of absorption coefficients    
%                as a function of wavenumber_1_per_cm    
%     molecule_number    
%              - HITRAN "Molecule number"
%     isotopologue_number    
%              - HITRAN "Isotopologue number"
%     line_center_wavenumber_1_per_cm    
%              - HITRAN "Vacuum wavenumber"
%     line_strength_at_reference_temperature_cm_per_molecule    
%              - HITRAN "Intensity"
%     A_coefficient_1_per_s    
%              - HITRAN "Einstein A coefficient"
%     air_broadened_half_width_1_per_cm_per_atm    
%              - HITRAN "Air-broadened half-width"
%     self_broadened_half_width_1_per_cm_per_atm    
%              - HITRAN "Self-broadened half-width"
%     air_pressure_induced_line_shift_1_per_cm_per_atm    
%              - HITRAN "Air pressured-induce line shift"
%     voigt_half_width_1_per_cm    
%              - The Voigt profile as a function of inputs at 296 K, 
%                as calculated using the formulae in:
%                M. Gharavi and S. Buckley, "Single Diode Laser
%                Sensor for Wide-Range H2O Temperature 
%                Measurements," Appl. Spectrosc. 58, 468 (2004).
%     peak_absorbance_    
%              - The peak absorbance of each line, as calculated using the
%                formulae in:
%                M. Gharavi and S. Buckley, "Single Diode Laser 
%                Sensor for Wide-Range H2O Temperature 
%                Measurements," Appl. Spectrosc. 58, 468 (2004).
%

CH4_Total_Partition_Function = importdata('12CH4_Partition_Extrapolate.txt');
N2O_Total_Partition_Function = importdata('N2O_Partition_Extrapolate.txt');
c = 29979245800;
d2 = 1.4387770;
Q_Tref = [5.0018e+03 5.9052e+02]; % for [N2O CH4], main isotopologue, 296K]
partial_pressure_molecules_1_per_cm3 = partial_pressure_atm * 2.48e19;
T_ref = 296;

molecular_weight_amu = 50;%molecular_weight_array_amu( ...hitran_struct.iso == isotopologues_);
doppler_half_width_1_per_cm = 3.581e-7 * hitran_struct.wnum * sqrt(T./molecular_weight_amu);

line_strength_temperature_dependence_N2O = hitran_struct.int.* (Q_Tref(1)./N2O_Total_Partition_Function(T,2)).*((exp(-d2.*hitran_struct.els./T))./(exp(-d2.*hitran_struct.els./T_ref)))...
        .* ((1-exp(-d2.*hitran_struct.wnum./T))./(1-exp(-d2.*hitran_struct.wnum./T_ref)));
temp_dependent_s = line_strength_temperature_dependence_N2O;

lorentzian_half_width_1_per_cm = (T_ref/T).^(hitran_struct.abcoef).*((hitran_struct.abroad ...
        * (pressure_atm - partial_pressure_atm)) ...
        + (hitran_struct.sbroad * partial_pressure_atm));

% Construct the spectrum
linestrength = temp_dependent_s;
gaussianBroad = doppler_half_width_1_per_cm;
lorentzianBroad = lorentzian_half_width_1_per_cm;

crossSection = zeros(size(wavenum));
for i = 1:numel(wavenum)
    idx = abs(hitran_struct.wnum- wavenum(i)) < 1;
    %crossSection(i) = sum(transitionStrength(idx).*areaNormalizedGaussian(wavenum(i),transitionWavenum(idx),gaussianBroad(idx)));
    %crossSection(i) = sum(transitionStrength(idx).*areaNormalizedLorentzian(wavenum(i),transitionWavenum(idx),lorentzianBroad(idx)));
    %crossSection(i) = sum(transitionStrength(idx).*areaNormVoigt(wavenum(i)-transitionWavenum(idx),gaussianBroad(idx),lorentzianBroad(idx)));
    crossSection(i) = sum(linestrength(idx).*areaNormPseudoVoigt(wavenum(i)-hitran_struct.wnum(idx),gaussianBroad(idx),lorentzianBroad(idx)));
end

absorbance_ = crossSection * path_length_cm * partial_pressure_molecules_1_per_cm3;
end

%---------- SPECTRUM GENERATION FUNCTIONS ----------------

function y = areaNormGaussian(x,FWHM)
% FWHM: Gaussian FWHM

    if FWHM == 0
        y = zeros(size(x));
        return
    end
    
    y = exp(-x.^2./(FWHM./sqrt(2.*log(2))./2).^2./2)./(FWHM./sqrt(2.*log(2))./2)./sqrt(2.*pi);
end

function y = areaNormLorentzian(x,FWHM)
% FWHM: Lorentzian FWHM

    y = FWHM./2./pi./(x.^2+FWHM.^2./4);
	
end

function y = areaNormVoigt(x,FWHMg,FWHMl)
% FWHMg: Gaussian FWHM
% FWHMl: lorentzian FWHM

    z = sqrt(log(2)).*complex(x,FWHMl./2)./(FWHMg./2);
    y = sqrt(log(2)).*cef(z,1000)./sqrt(pi)./(FWHMg./2);
	
end

function y = areaNormPseudoVoigt(x,GammaG,GammaL)
    % GammaG: Gaussian FWHM
    % GammaL: Lorentz FWHM
    
    Gamma = (GammaG.^5 + ...
        2.69269.*GammaG.^4.*GammaL + ...
        2.42843.*GammaG.^3.*GammaL.^2 + ...
        4.47163.*GammaG.^2.*GammaL.^3 + ...
        0.07842.*GammaG.^1.*GammaL.^4 + ...
        GammaL.^5).^(1./5);
    GammaRatio = GammaL./Gamma;
    eta = 1.36603.*GammaRatio - ...
        0.47719.*GammaRatio.^2 + ...
        0.11116.*GammaRatio.^3;

    y = (1-eta).*areaNormGaussian(x,GammaG) + ...
        eta.*areaNormLorentzian(x,GammaL);
		
end

function w = cef(z,N)

	%  Computes the function w(z) = exp(-z^2) erfc(-iz) using a rational 
	%  series with N terms.  It is assumed that Im(z) > 0 or Im(z) = 0.
	%
	%                             Andre Weideman, 1995

	M = 2*N;  M2 = 2*M;  k = [-M+1:1:M-1]';    % M2 = no. of sampling points.
	L = sqrt(N/sqrt(2));                       % Optimal choice of L.
	theta = k*pi/M; t = L*tan(theta/2);        % Define variables theta and t.
	f = exp(-t.^2).*(L^2+t.^2); f = [0; f];    % Function to be transformed.
	a = real(fft(fftshift(f)))/M2;             % Coefficients of transform.
	a = flipud(a(2:N+1));                      % Reorder coefficients.
	Z = (L+1i*z)./(L-1i*z); p = polyval(a,Z);    % Polynomial evaluation.
	w = 2*p./(L-1i*z).^2+(1/sqrt(pi))./(L-1i*z); % Evaluate w(z).

	% abs_carrier = exp(-4*(real(z)).^2*log(2));

end