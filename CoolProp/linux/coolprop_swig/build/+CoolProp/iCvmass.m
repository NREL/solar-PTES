function v = iCvmass()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 41);
  end
  v = vInitialized;
end
