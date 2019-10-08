function v = isurface_tension()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 48);
  end
  v = vInitialized;
end
