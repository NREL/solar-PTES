function v = iPH()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 74);
  end
  v = vInitialized;
end
