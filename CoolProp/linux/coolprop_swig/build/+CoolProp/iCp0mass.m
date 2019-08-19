function v = iCp0mass()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 40);
  end
  v = vInitialized;
end
