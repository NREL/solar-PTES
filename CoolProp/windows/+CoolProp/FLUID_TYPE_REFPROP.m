function v = FLUID_TYPE_REFPROP()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 94);
  end
  v = vInitialized;
end
