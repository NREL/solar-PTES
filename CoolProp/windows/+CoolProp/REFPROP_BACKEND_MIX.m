function v = REFPROP_BACKEND_MIX()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 149);
  end
  v = vInitialized;
end
