function v = DmolarQ_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 105);
  end
  v = vInitialized;
end
