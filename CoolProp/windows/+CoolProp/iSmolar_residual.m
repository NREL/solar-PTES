function v = iSmolar_residual()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 36);
  end
  v = vInitialized;
end
