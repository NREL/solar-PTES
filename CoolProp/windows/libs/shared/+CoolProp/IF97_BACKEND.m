function v = IF97_BACKEND()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 151);
  end
  v = vInitialized;
end
