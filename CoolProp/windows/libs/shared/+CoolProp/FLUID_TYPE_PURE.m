function v = FLUID_TYPE_PURE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 92);
  end
  v = vInitialized;
end
