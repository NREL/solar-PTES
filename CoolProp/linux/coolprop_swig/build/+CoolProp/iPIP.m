function v = iPIP()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 63);
  end
  v = vInitialized;
end
