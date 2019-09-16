function v = SRK_BACKEND()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 155);
  end
  v = vInitialized;
end
