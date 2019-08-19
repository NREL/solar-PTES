function v = IFRAC_UNDEFINED()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 88);
  end
  v = vInitialized;
end
