function v = FLUID_TYPE_REFPROP()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 92);
  end
  v = vInitialized;
end
