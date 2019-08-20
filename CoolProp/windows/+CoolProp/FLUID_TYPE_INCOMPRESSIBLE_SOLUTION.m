function v = FLUID_TYPE_INCOMPRESSIBLE_SOLUTION()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 96);
  end
  v = vInitialized;
end
