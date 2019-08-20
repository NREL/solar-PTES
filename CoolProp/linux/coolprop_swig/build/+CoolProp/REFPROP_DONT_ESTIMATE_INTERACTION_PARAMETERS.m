function v = REFPROP_DONT_ESTIMATE_INTERACTION_PARAMETERS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 139);
  end
  v = vInitialized;
end
