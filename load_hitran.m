function hitran_struct = load_hitran(file_name, wavenumber_1_per_cm, ...
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

%%%
% Setup


CH4_Total_Partition_Function = importdata('12CH4_Partition_Extrapolate.txt');
N2O_Total_Partition_Function = importdata('N2O_Partition_Extrapolate.txt');
c = 29979245800;
d2 = 1.4387770;
Q_Tref = [5.0018e+03 5.9052e+02]; % for [N2O CH4], main isotopologue, 296K]
N = numel(wavenumber_1_per_cm);
number_of_lines_to_read_per_pass = 100000;
number_of_lines_from_previous_passes = 0;
fid = fopen(file_name);
T_ref = 296;
hitran_struct = struct();
hitran_struct.absorption_coefficient_1_per_cm = zeros(1, N);
hitran_struct.line_strength_temperature_dependence_N2O = zeros(1,N);
hitran_struct.line_strength_temperature_dependence_CH4 = zeros(1,N);
hitran_struct.temp_dependent_s = zeros(1,N);
hitran_struct.molecule_number = [];
hitran_struct.isotopologue_number = [];
hitran_struct.line_center_wavenumber_1_per_cm = [];
hitran_struct.line_strength_at_reference_temperature_cm_per_molecule = [];
hitran_struct.t_ref = [];
hitran_struct.A_coefficient_1_per_s = [];
hitran_struct.air_broadened_half_width_1_per_cm_per_atm = [];
hitran_struct.self_broadened_half_width_1_per_cm_per_atm = [];
hitran_struct.lower_state_energy = [];
hitran_struct.temperature_dependent_width_air = [];
hitran_struct.air_pressure_induced_line_shift_1_per_cm_per_atm = [];
hitran_struct.voigt_half_width_1_per_cm = [];
hitran_struct.peak_absorbance_ = [];
hitran_struct.upper_state_degeneracy = [];
hitran_struct.absorption_cross_section_at_line_center_cm2_per_molecule = zeros(1, N);
hitran_struct.absorption_cross_section_cm2_per_molecule = zeros(1,N);



while(true)
	%%%
	% Load data
	%%%
	HITRAN_data = textscan(fid, ...
			['%2c' '%1f' '%12f' '%10f' '%10f' ...
			'%5f' '%5f' '%10f' '%4f' '%8f' ...
			'%15c' '%15c' '%15c' '%15c' '%6c' ...
			'%12c' '%1c' '%7f' '%7f'], ...
			number_of_lines_to_read_per_pass, 'delimiter', '', ...
			'whitespace', '');
			
	%%%
	% Done condition
	%%%
	HITRAN_data{1} = str2num(HITRAN_data{1});
	number_of_lines_read_this_pass = numel(HITRAN_data{1});
	if(number_of_lines_read_this_pass == 0)
		break;
	end
	
	%%%
	% Add lines into datasets
	%%%
	hitran_struct.molecule_number = [hitran_struct.molecule_number ...
			HITRAN_data{1}.'];
	hitran_struct.isotopologue_number = [hitran_struct.isotopologue_number ...
			HITRAN_data{2}.'];
	hitran_struct.line_center_wavenumber_1_per_cm = ...
			[hitran_struct.line_center_wavenumber_1_per_cm ...
			HITRAN_data{3}.'];
    hitran_struct.line_strength_at_reference_temperature_cm_per_molecule = ...
			[hitran_struct.line_strength_at_reference_temperature_cm_per_molecule ...
			HITRAN_data{4}.'];
	hitran_struct.A_coefficient_1_per_s = ...
			[hitran_struct.A_coefficient_1_per_s ...
			HITRAN_data{5}.'];
	hitran_struct.air_broadened_half_width_1_per_cm_per_atm = ...
			[hitran_struct.air_broadened_half_width_1_per_cm_per_atm ...
			HITRAN_data{6}.'];
	hitran_struct.self_broadened_half_width_1_per_cm_per_atm = ...
			[hitran_struct.self_broadened_half_width_1_per_cm_per_atm ...
			HITRAN_data{7}.'];
    hitran_struct.lower_state_energy = ...
			[hitran_struct.lower_state_energy ...
			HITRAN_data{8}.'];
    hitran_struct.temperature_dependent_width_air = ...
			[hitran_struct.temperature_dependent_width_air ...
			HITRAN_data{9}.'];
	hitran_struct.air_pressure_induced_line_shift_1_per_cm_per_atm = ...
			[hitran_struct.air_pressure_induced_line_shift_1_per_cm_per_atm ...
			HITRAN_data{10}.'];
    hitran_struct.upper_state_degeneracy = [hitran_struct.upper_state_degeneracy ...
			HITRAN_data{18}.'];
    hitran_struct.t_ref = hitran_struct.line_strength_at_reference_temperature_cm_per_molecule;

% Calculate the temperature dependence of S, using HITRAN forms at http://hitran.org/docs/definitions-and-units/
% Total partition function is obtained from extrapolation from data found at http://hitran.org/docs/iso-meta/
% Some data from HITRAN are 0 (E" values), need to remove NaN resulting
% from calculation
   
    hitran_struct.line_strength_temperature_dependence_CH4 = hitran_struct.t_ref.* (Q_Tref(2)./CH4_Total_Partition_Function(T,2)).*((exp(-d2.*hitran_struct.lower_state_energy./T))./(exp(-d2.*hitran_struct.lower_state_energy./T_ref)))...
            .* ((1-exp(-d2.*hitran_struct.line_center_wavenumber_1_per_cm./T))./(1-exp(-d2.*hitran_struct.line_center_wavenumber_1_per_cm./T_ref)));
        hitran_struct.temp_dependent_s =  hitran_struct.line_strength_temperature_dependence_CH4;
        
    
    hitran_struct.line_strength_temperature_dependence_N2O = hitran_struct.t_ref.* (Q_Tref(1)./N2O_Total_Partition_Function(T,2)).*((exp(-d2.*hitran_struct.lower_state_energy./T))./(exp(-d2.*hitran_struct.lower_state_energy./T_ref)))...
            .* ((1-exp(-d2.*hitran_struct.line_center_wavenumber_1_per_cm./T))./(1-exp(-d2.*hitran_struct.line_center_wavenumber_1_per_cm./T_ref)));
        hitran_struct.temp_dependent_s = hitran_struct.line_strength_temperature_dependence_N2O;

	lines_to_delete = [];
    
	
	for(set_line_number = 1:number_of_lines_read_this_pass)
		line_number = set_line_number + number_of_lines_from_previous_passes;
		%%%
		% Skip isotopologues that are not of interest
		%%%
		if(~ismember(hitran_struct.isotopologue_number(line_number), ...
				isotopologues_))
			lines_to_delete = union(lines_to_delete, line_number);
			continue;
		end
		%%%
		% Skip lines out of wavenumber range
		%%%
		if(hitran_struct.line_center_wavenumber_1_per_cm(line_number) ...
				< min(wavenumber_1_per_cm) || ...
				max(wavenumber_1_per_cm) < ...
				hitran_struct.line_center_wavenumber_1_per_cm(line_number))
			lines_to_delete = union(lines_to_delete, line_number);
			continue;
		end
			
		%%%
		% Calculate absorption coefficient
	
		
		% Temperature dependence of width is accounted for

        % Lorentzian width
		lorentzian_half_width_1_per_cm = (T_ref/T).^(hitran_struct.temperature_dependent_width_air(line_number)).*((hitran_struct.air_broadened_half_width_1_per_cm_per_atm(line_number) ...
				* (pressure_atm - partial_pressure_atm)) ...
				+ (hitran_struct.self_broadened_half_width_1_per_cm_per_atm(line_number) * ...
				partial_pressure_atm));

        % Doppler width
		molecular_weight_amu = molecular_weight_array_amu( ...
				hitran_struct.isotopologue_number(line_number) == isotopologues_);
		doppler_half_width_1_per_cm = 3.581e-7 * ...
				hitran_struct.line_center_wavenumber_1_per_cm(line_number) * ...
				sqrt(T/molecular_weight_amu);
				
		hitran_struct.voigt_half_width_1_per_cm ...
				= [hitran_struct.voigt_half_width_1_per_cm ...
				(0.5346 * lorentzian_half_width_1_per_cm) ...
				+ sqrt(0.2166 * (lorentzian_half_width_1_per_cm ^ 2) ...
				+ (doppler_half_width_1_per_cm ^ 2))];
		
		% pV = NkT -> p = nkT -> n = p/kT -> 
		% n (1/cm^3) = n (1/m^3) * 1e-6 (m/cm)^3 = p (Pa) / (1.38e-23 (J/K) * Tref (K)) * 1e-6 (m/cm)^3 
		% = p (atm) * (101325 Pa/atm) / (1.3806e-23 (J/K) * Tref (K)) * 1e-6 (m/cm)^3 
		% = p (atm) * 101325 * 1e-6 / (1.3806e-23 * 296) * (Pa * m^3 / J) * (1/(atm * cm^3))
		% = p (atm) * 2.48e19 * (1/(atm * cm^3))
	
		partial_pressure_molecules_1_per_cm3 = partial_pressure_atm * 2.48e19;
				
		x = (lorentzian_half_width_1_per_cm / ...
				hitran_struct.voigt_half_width_1_per_cm(end));
		y = (abs(wavenumber_1_per_cm - hitran_struct.line_center_wavenumber_1_per_cm(line_number) ...
				- hitran_struct.air_pressure_induced_line_shift_1_per_cm_per_atm(line_number) * pressure_atm) ...
				/ hitran_struct.voigt_half_width_1_per_cm(end));
		
        
        % Temperature dependence linestrength, S
%         
%          hitran_struct.line_strength_temperature_dependence_N2O = (hitran_struct.A_coefficient_1_per_s./(8*pi*c*hitran_struct.line_center_wavenumber_1_per_cm.^2))...
%             .*hitran_struct.upper_state_degeneracy.*((exp(-c2./hitran_struct.lower_state_energy)./T)./(1-exp(c2./hitran_struct.line_center_wavenumber_1_per_cm)./T))...
%             ./total_partition_function_N2O(T);
%          
%      
%         hitran_struct.absorption_coefficient_1_per_cm_temperature_dependence_CH4 = hitran_struct.line_strength_at_reference_temperature_cm_per_molecule...
%             .* (Q_Tref(2)./total_partition_function_CH4(T)) .*((exp(-c2./hitran_struct.lower_state_energy)./T)./(exp(-c2./hitran_struct.lower_state_energy)./T_ref))...
%             .* ((1-exp(-c2.*hitran_struct.line_center_wavenumber_1_per_cm./T))./(1-exp(-c2.*hitran_struct.line_center_wavenumber_1_per_cm./T_ref)));
%           
% 		absorption_cross_section_at_line_center_cm2_per_molecule = ...
% 				hitran_struct.line_strength_at_reference_temperature_cm_per_molecule(line_number) ...
% 				./ (2 * hitran_struct.voigt_half_width_1_per_cm(end) * (1.065 + (0.447 * x) + (0.058 * (x .^ 2))));

    
	absorption_cross_section_at_line_center_cm2_per_molecule = ...
				   hitran_struct.temp_dependent_s(line_number) ...
				./ (2 .* hitran_struct.voigt_half_width_1_per_cm(end) .* (1.065 + (0.447 .* x) + (0.058 .* (x .^ 2))));
            
            
    hitran_struct.absorption_cross_section_at_line_center_cm2_per_molecule = absorption_cross_section_at_line_center_cm2_per_molecule;	
    

		absorption_cross_section_cm2_per_molecule = ...
				absorption_cross_section_at_line_center_cm2_per_molecule ...
				.* ( ((1 - x) .* exp(-0.693 .* (y .^ 2))) + ...
				x ./ (1 + (y .^ 2)) + (0.016 * ...
				(1 - x) .* x .* ...
				(exp( -0.0841 .* (y .^ 2.25)) - 1 ./ (1 + 0.021 .* (y .^ 2.25)))) );

        % alternative Voigt calculation, but the lineshape looks too
        % Lorentzian even at low pressures, something is wrong


% 		absorption_cross_section_cm2_per_molecule = ...
% 				absorption_cross_section_at_line_center_cm2_per_molecule.*voigtf(x,y,1);


        % end of alternative Voigt calculation

            
        hitran_struct.absorption_cross_section_cm2_per_molecule = absorption_cross_section_cm2_per_molecule;
                  
		hitran_struct.absorption_coefficient_1_per_cm = hitran_struct.absorption_coefficient_1_per_cm + partial_pressure_molecules_1_per_cm3 .* absorption_cross_section_cm2_per_molecule;

	end
	%%%
	% Delete lines not in set of desired isotopologues
	%%%
	hitran_struct.molecule_number(lines_to_delete) = [];
	hitran_struct.isotopologue_number(lines_to_delete) = [];
	hitran_struct.line_center_wavenumber_1_per_cm(lines_to_delete) = [];
    hitran_struct.t_ref(lines_to_delete) = [];
	hitran_struct.line_strength_at_reference_temperature_cm_per_molecule(lines_to_delete) = [];
    hitran_struct.line_strength_temperature_dependence_N2O(lines_to_delete) = [];
    hitran_struct.line_strength_temperature_dependence_CH4(lines_to_delete) = [];
	hitran_struct.A_coefficient_1_per_s(lines_to_delete) = [];
    hitran_struct.air_broadened_half_width_1_per_cm_per_atm(lines_to_delete) = [];
	hitran_struct.self_broadened_half_width_1_per_cm_per_atm(lines_to_delete) = [];
	hitran_struct.air_pressure_induced_line_shift_1_per_cm_per_atm(lines_to_delete) = [];
	number_of_lines_from_previous_passes = numel(hitran_struct.molecule_number);

end

hitran_struct.absorbance_ = ...
		hitran_struct.absorption_coefficient_1_per_cm * path_length_cm;
hitran_struct.transmission_ = exp(-hitran_struct.absorbance_);

fclose(fid);
end