function v = FLUID_TYPE_PURE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 90);
  end
  v = vInitialized;
end
