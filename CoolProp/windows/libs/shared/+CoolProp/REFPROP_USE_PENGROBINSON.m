function v = REFPROP_USE_PENGROBINSON()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 167);
  end
  v = vInitialized;
end
