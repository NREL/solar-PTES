function v = REFPROP_DONT_ESTIMATE_INTERACTION_PARAMETERS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 165);
  end
  v = vInitialized;
end
