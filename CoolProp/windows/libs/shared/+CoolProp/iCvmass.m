function v = iCvmass()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 42);
  end
  v = vInitialized;
end
