function v = iHH()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 73);
  end
  v = vInitialized;
end
