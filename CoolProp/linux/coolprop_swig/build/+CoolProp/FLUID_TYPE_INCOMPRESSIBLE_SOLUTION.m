function v = FLUID_TYPE_INCOMPRESSIBLE_SOLUTION()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 94);
  end
  v = vInitialized;
end
