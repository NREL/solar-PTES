function v = iphase_liquid()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 78);
  end
  v = vInitialized;
end
