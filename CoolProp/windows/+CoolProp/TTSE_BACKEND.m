function v = TTSE_BACKEND()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 153);
  end
  v = vInitialized;
end
