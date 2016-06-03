%%load_hitran

T = 10;  % set temperature of gas for simulation
pressure_torr = 1; % input pressure in Torr for simulation
pressure_atm = pressure_torr / 760;
mole_fraction = 1; 
partial_pressure_atm = pressure_atm * mole_fraction;
path_length_cm = 5; % absorption pathlength
colors = discrete_plot_colors_2();
load_file_name = 'N2O_2725_2850.par'; % input file name from HITRAN
const = phys_const_2();
N = 50000; % set this number higher to increase spectral resolution
wavenumber_max = 2850;
wavenumber_min = 2725;
wavelength_max_nm = 1e7/wavenumber_min ;
wavelength_min_nm = 1e7/wavenumber_max;

frequency_samples_wavenumber = 1e7./([wavelength_max_nm wavelength_min_nm]);

df_wavenumber = (frequency_samples_wavenumber(2) - frequency_samples_wavenumber(1)) / (N - 1);

frequency_samples_wavenumber = frequency_samples_wavenumber(1): df_wavenumber:frequency_samples_wavenumber(2);

wavelength_samples_nm = 1e7 ./ frequency_samples_wavenumber;

wavenumber_samples_1_per_cm = 1 ./ wavelength_samples_nm * 1e7;

f0_wavenumber = min(frequency_samples_wavenumber) + (df_wavenumber * N / 2);


% if(1) % 14N2O
% 	isotopologues_array_ = [1];
% 	molecular_weight_array_amu = [14 + 14 + 16];
% else % 15N2O
% 	isotopologues_array_ = [2]; 
% 	molecular_weight_array_amu = [14 + 15 + 16];
% end

% 14N2O mass

	isotopologues_array_ = [1];
	molecular_weight_array_amu = [14 + 14 + 16];

%%%
% Execute
%%%

if(1)
	hitran_struct = load_hitran(load_file_name, wavenumber_samples_1_per_cm, ...
			pressure_atm, partial_pressure_atm, ...
			path_length_cm, isotopologues_array_, molecular_weight_array_amu, T);
end	
	
%%%
% Plot vs. wavelength

absorbance = hitran_struct.absorbance_;

if(1)
	figure;
	plot(frequency_samples_wavenumber, absorbance);

    	xlabel('Wavenumber (cm^-^1)');

	ylabel('Absorbance');
	title('N_2O Absorbance');
else
	figure;
    	plot(frequency_samples_wavenumber, absorbance);

end
