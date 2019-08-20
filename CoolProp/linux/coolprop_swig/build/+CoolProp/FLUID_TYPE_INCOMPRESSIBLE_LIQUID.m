function v = FLUID_TYPE_INCOMPRESSIBLE_LIQUID()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 93);
  end
  v = vInitialized;
end
