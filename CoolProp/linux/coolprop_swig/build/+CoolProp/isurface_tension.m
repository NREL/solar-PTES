function v = isurface_tension()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 46);
  end
  v = vInitialized;
end
