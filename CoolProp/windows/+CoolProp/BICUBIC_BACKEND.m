function v = BICUBIC_BACKEND()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 154);
  end
  v = vInitialized;
end
