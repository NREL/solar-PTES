function v = iCpmass()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 40);
  end
  v = vInitialized;
end
