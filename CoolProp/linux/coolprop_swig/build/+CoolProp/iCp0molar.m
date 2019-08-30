function v = iCp0molar()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 31);
  end
  v = vInitialized;
end