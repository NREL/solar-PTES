function v = ALTERNATIVE_REFPROP_HMX_BNC_PATH()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 138);
  end
  v = vInitialized;
end
