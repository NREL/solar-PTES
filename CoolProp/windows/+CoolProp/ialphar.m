function v = ialphar()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 54);
  end
  v = vInitialized;
end
