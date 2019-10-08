function v = REFPROP_USE_GERG()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 166);
  end
  v = vInitialized;
end
