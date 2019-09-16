function v = BICUBIC_BACKEND_FAMILY()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 141);
  end
  v = vInitialized;
end
