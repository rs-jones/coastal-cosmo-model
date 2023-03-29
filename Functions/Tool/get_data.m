%
% sample_data = get_data(input_name)
% sample_data = get_data(input_name,inheritance)
% 
% Loads input data (either from an .xlsx file, or headerless .txt or .csv), 
% extracts sample names and concentrations, sorts data for functions, and 
% exports to single output.
%
% inheritance is an optional two-part (struct) input containing 
% inheritance.conc_mean and inheritance.conc_uncert, which specify the 
% nuclide concentration of an inheritance sample (atoms/g)
%
% For CRONUScalc calculations, the erosion rate and inheritance is set to 
% zero. Attenuation lengths are automatically determined using the Sato 
% model, based on the location information of each sample.
%
% Input data must include the following information:
% 1. Sample name
% 2. Latitude (decimal degrees)
% 3. Longitude (decimal degrees)
% 4. Elevation (m asl)
% 5. Pressure (hPa) (zero if not known)
% 6. Distance along transect from cliff (m)
% 7. Sample thickness (cm)
% 8. Bulk density (g/cm^3)
% 9. Shielding factor for terrain, snow, etc. (unitless)
% 10. Sample 10-Be concentration, mean (atoms of 10-Be/g)
% 11. Sample 10-Be concentration, 1 sigma uncertainty (atoms of 10-Be/g)
% 12. Year the sample was collected (calendar year)
%
% Written by Richard Selwyn Jones, Monash University
% richard.s.jones@monash.edu
% Modification of function within the iceTEA tools suite, which is built 
% on versions of CRONUS-Earth and CRONUScalc code.
%
%
%%

function sample_data = get_data(input_name,inheritance)

  if nargin<2
      inheritance.conc_mean = 0;
      inheritance.conc_uncert = 0;
  end

  if ispc
      input_path = strcat(pwd,'\',input_name);
  else
      input_path = strcat(pwd,'/',input_name);
  end
  [~,~,ext] = fileparts(input_path);
  if isempty(ext)
      error('File name of sample data does not include the file type extension (e.g. .xlsx)')
  end

  % Load input data
  if strcmpi(ext,'.xlsx') || strcmpi(ext,'.xls')
      [in_data,in_txt,in_raw] = xlsread(input_name);
      for a = 1:length(in_raw(1,:))
          this_header = in_raw{1,a};
          strrow(a) = ~isnumeric(this_header);
      end
      n_raw_columns = numel(in_raw(1,:));
      n_data_columns = numel(in_data(1,:));
      % Remove any empty columns
      if (n_raw_columns > 12) && all(~strrow(13:end))
          in_raw = in_raw(:,1:12);
      elseif (n_raw_columns > 12) && (n_data_columns < 12)
          in_raw = in_raw(:,1:12);
      end
      for a = 1:length(in_raw(1,:))
          this_header = in_raw{1,a};
          new_strrow(a) = ~isnumeric(this_header);
      end
      if all(new_strrow) % Remove header
          in_raw = in_raw(2:end,:); 
          in_txt = in_txt(2:end,:);
      end
      
  elseif strcmpi(ext,'.txt') || strcmpi(ext,'.csv')
      warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
      in_raw = table2cell(readtable(input_name));
      in_txt = cell(size(in_raw));
      in_data = zeros(size(in_raw));
      for a = 1:length(in_raw(:,1))
          for b = 1:length(in_raw(1,:))
              this_raw = in_raw{a,b};
              this_num_log = isnumeric(this_raw);
              if this_num_log
                  in_data(a,b) = this_raw;
              else
                  in_txt{a,b} = this_raw;
              end
          end
      end
      in_data = in_data(:,2:end);
  end
  
  n_raw_columns = numel(in_raw(1,:));
  n_data_columns = numel(in_data(1,:));
  
  % Check inputs
  if (n_raw_columns > 12)
      error('sample input data has too many fields!');
  end
  
  % Check sample names
  col1_nan_log = isnan(in_data(:,1)); % Check for NaNs in column 1
  if any(col1_nan_log) && (n_data_columns == 12)
      miss_names = cellstr(num2str(in_data(~col1_nan_log,1))); % Get and convert numeric names from in_data
      in_txt_names = in_raw(:,1);
      in_txt_names(~col1_nan_log,1) = miss_names; % Add missed non-numeric names to data
      in_txt = cell(size(in_raw));
      in_txt(1,:) = in_raw(1,:); in_txt(:,1) = in_txt_names; % Combine for new text data
      in_data = in_data(:,2:end); % Adjust in_data (remove names column)
      warning('some sample names are numeric... fixed.');
  end
  
  % Set data columns
  name_txt_c = 1;
  lat_c = 1;
  lon_c = 2;
  elev_c = 3;
  pressure_c = 4;
  distance_c = 5;
  thk_c = 6;
  density_c = 7;
  shield_c = 8;
  Be_c = 9;
  Be_sig_c = 10;
  year_c = 11;
  
  
  % Get number of samples (and samples with 10-Be and 26-Al concs.)
  n_samples = numel(in_data(:,lat_c));
  
  % Get sample names
  names = in_txt(:,name_txt_c)';
  
  
  % If concentrations are NaN, make zero
  for a = 1:n_samples
      if (isnan(in_data(a,Be_c)))
          in_data(a,Be_c) = 0;
          in_data(a,Be_sig_c) = 0;
      end
  end
 
  % Sort details for each sample and nuclide (currently only 10Be and 26Al)
  logical_10 = any(in_data(:,Be_c),2)';
  
  
  % Get raw measurement data
  meas_raw_mean = in_data(:,Be_c);
  meas_raw_uncert = in_data(:,Be_sig_c);
  
  % Correct raw data for inheritance
  meas_inhcorr_mean = meas_raw_mean - inheritance.conc_mean;
  meas_inhcorr_uncert = sqrt(meas_raw_uncert.^2 + inheritance.conc_uncert.^2);
  

  % Collate sample details
  n_data = numel(in_data(:,1));
  for c = 1:n_data
      
      sample_data.s{c}.name = names(c);
      
      if logical_10(c)
          sample_data.s{c}.nuclide10 = 1;
          sample_data.s{c}.N10 = meas_inhcorr_mean(c);
          sample_data.s{c}.dN10 = meas_inhcorr_uncert(c);  
      else
          sample_data.s{c}.nuclide10 = 0;
      end           
  end
   
  
  % Assume certain sample details are unknown or zero
  e_rate = zeros(n_data,1); % Erosion rate (mm/kyr)
  inh10 = zeros(n_data,1);  % 10Be inheritance (atoms/g)
  inh26 = zeros(n_data,1);  % 26Al inheritance (atoms/g)
  top_depth_gcm2 = zeros(n_data,1); % Top depth (g/cm^2) - i.e. the surface
  
  % Initially set attenuation length to zero
  init_L = zeros(n_data,1);
  
  % Combine necessary data for CronusCalc
  Al_N = zeros(size(meas_raw_mean));
  sample_data.CC = [in_data(:,[lat_c:pressure_c,thk_c:shield_c]),e_rate,meas_raw_mean,Al_N,inh10,inh26,init_L,top_depth_gcm2,in_data(:,year_c)];
  
  % Determine atmospheric pressure (if unknown)
  sample_data.CC = atm_pressure(in_data(:,lat_c),in_data(:,lon_c),sample_data.CC);
  

  % Determine attenuation lengths (g/cm^2) based on modified Sato model
  for d = 1:length(sample_data.CC(:,1))
      L(d,1) = attenuationlength(sample_data.CC(d,1),sample_data.CC(d,2),sample_data.CC(d,3),sample_data.CC(d,4));
  end
  sample_data.CC(:,13) = L; % Add attenuation lengths
  
  
  % Create uncertainties struct
  
  % Create a standard elevation uncertainty
  elev_uncert = 0.5; % 0.5 metres
  uncert = zeros(size(sample_data.CC)); % Assume zero uncertainty for now
  Al_N_uncert = zeros(size(meas_raw_uncert));
  uncert(:,[9,10]) = [meas_raw_uncert,Al_N_uncert]; % Add nuclide concentration uncertainties
  uncert(:,3) = elev_uncert; % Add elevation uncertainties
  
  % Determine pressure uncertainties from elevation
  upelev_uncert = uncert;  upelev_uncert(:,3) = sample_data.CC(:,elev_c) + elev_uncert';
  upelev_uncert = atm_pressure(in_data(:,lat_c),in_data(:,lon_c),upelev_uncert);
  uppress_diff = sample_data.CC(:,4) - upelev_uncert(:,4);
  loelev_uncert = uncert;  loelev_uncert(:,3) = sample_data.CC(:,elev_c) - elev_uncert';
  loelev_uncert = atm_pressure(in_data(:,lat_c),in_data(:,lon_c),loelev_uncert);
  lopress_diff = sample_data.CC(:,4) - loelev_uncert(:,4);
  press_uncert = mean([abs(uppress_diff),abs(lopress_diff)],2); % Get average of upper and lower pressure uncertainties
  uncert(:,4) = press_uncert;
  
  sample_data.CC_uncert = uncert;
  
  
  % Export logical of nuclides
  sample_data.logical_10 = logical_10;
  
  % Export sample position
  sample_data.position = in_data(:,distance_c)';
  
  % Create default cover depth
  sample_data.cover.z = 0;
  
  % Export raw and inheritance-corrected concentrations
  inh_log = inheritance.conc_mean==meas_raw_mean;
  if any(inh_log)
      inhcorr_pos = sample_data.position(~inh_log)';
      meas_inhcorr_mean = meas_inhcorr_mean(~inh_log);
      meas_inhcorr_uncert = meas_inhcorr_uncert(~inh_log);
  else
      inhcorr_pos = sample_data.position';
  end
  sample_data.measured_raw = [sample_data.position',meas_raw_mean,meas_raw_uncert];
  sample_data.measured_inh = [inhcorr_pos,meas_inhcorr_mean,meas_inhcorr_uncert];
  

end
