function v = FLUID_TYPE_UNDEFINED()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 95);
  end
  v = vInitialized;
end
