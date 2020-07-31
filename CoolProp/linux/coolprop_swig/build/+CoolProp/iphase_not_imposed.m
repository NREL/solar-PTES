function v = iphase_not_imposed()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 84);
  end
  v = vInitialized;
end
