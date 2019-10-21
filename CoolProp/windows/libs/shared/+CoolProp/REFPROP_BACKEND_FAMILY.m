function v = REFPROP_BACKEND_FAMILY()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 136);
  end
  v = vInitialized;
end
