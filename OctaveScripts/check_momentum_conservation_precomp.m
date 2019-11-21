% Check the momentum conservation for the precomputation output.
%
% Momentum convervation requires the condition
%
% \f[
% 0 = R_{D,m}^{\alpha\alpha'} + R_{I,m}^{\alpha'\alpha}
% \f]
% where R_{D,m}^{\alpha\alpha'} and R_{I,m}^{\alpha'\alpha} are related
% to the differential and integral part of the collision operator,
% respectively.
% They can be computed as
% \f[
% R_{D,m}^{\alpha\alpha'} = \frac{1}{T_\alpha} \sum\limits_{m'}^{} - V_{m'm}^{\alpha\alpha'} + D_{m'm}^{\alpha\alpha'}
% \f]
% and
% \f[
% R_{I,m}^{\alpha\alpha'} = \frac{2}{3} \frac{1}{T_\alpha} \sum\limits_{m'}^{} I_{m'm}^{\alpha\alpha'}
% \f]
% The three matrices V_{m'm}^{\alpha\alpha'}, D_{m'm}^{\alpha\alpha'}
% and I_{m'm}^{\alpha\alpha'} are connected with the matrices anumm_aa,
% denmm_aa and ailmm_aa, respectively. They are not equal, as the output
% has been multiplied with M_transform_inv. Note the M_transform_ is not
% the inverse of M_transform_inv, as for the former just an identiy
% matrix in the output, for unknown reasons.
%
% input
% -----
% precomp_filename: name of the file (+path) from which to read the
%   precomp data to check.
% temperatures: vector with the temperatures of the species. It does
%   not matter if it is a row or column vector, as it is transformed
%   to row vector.
%
% output
% ------
% res: R_D + R_I for the different lag values and species.
% res_rel: (R_D + R_I) / (|R_D| + |R_I|) for the different lag values
%   and species.
function [res, res_rel] = check_momentum_conservation_precomp(precomp_filename, temperatures)
  temperatures = temperatures(:);

  h = load(precomp_filename);

  num_species = size(temperatures, 1);


  for k = 1:num_species
    for l = 1:num_species
      dmm(:,:, k,l) = h.M_transform_inv \ h.denmm_aa(:,:,k,l);
      vmm(:,:, k,l) = h.M_transform_inv \ h.anumm_aa(:,:,k,l);
      % We need the (1) component, and need to remember that index 1 is
      % the (0) component
      imm(:,:, k,l) = h.M_transform_inv \ h.ailmm_aa(:,:,2,k,l);
    end
  end

  for k = 1:num_species
    for l = 1:num_species
      rd(:,k,l) =     sum(-vmm(:,:,k,l) + dmm(:,:,k,l), 1)/temperatures(k);
      ri(:,k,l) = 2/3*sum(imm(:,:,k,l), 1)/temperatures(k);
    end
  end

  for k = 1:num_species
    for l = 1:num_species
      norm=sum(rd(:,k,l));
      res(:,k,l) = (rd(:,k,l) + ri(:,l,k));
      %~ res_rel(:,k,l) = res(:,k,l) ./ (abs(rd(:,k,l)) + abs(ri(:,l,k)));
      res_rel(:,k,l) = res(:,k,l) ./ norm;
    end
  end
end
