%
% [param_array]=param_rand_gen(min_value,max_value,rand_number)
% 
% Generates an array of random values given a range for the parameter.
%
% Output is to be used to create combinations of parameters to model 
% nuclide concentrations.
%

%%

function [param_array]=param_rand_gen(min_value,max_value,rand_number)

  % Check inputs
  if (nargin ~= 3)
    error('param_rand_gen has wrong number of inputs!');
  end

  % Load random number generator
  %rng(0,'twister'); % Produces the same output within a Matlab session (good for debugging)
  rng('shuffle','twister'); % Produces a different output every time

  % Determine range
  range = max_value-min_value;

  % Generate values
  param_array = range.*rand(rand_number,1) + min_value;

end
