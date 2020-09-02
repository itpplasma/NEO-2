% Igichine model: del_s = 4 sqrt(s_res) ka / e
% ka=0.06
function del_s = Igichine_model_del_s(s_res, ka)
  del_s = 4*sqrt(s_res)*ka/exp(1);
end
