function v = SmolarT_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 110);
  end
  v = vInitialized;
end
