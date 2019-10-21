function v = iPhase()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 76);
  end
  v = vInitialized;
end
