function v = IFRAC_VOLUME()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 87);
  end
  v = vInitialized;
end
