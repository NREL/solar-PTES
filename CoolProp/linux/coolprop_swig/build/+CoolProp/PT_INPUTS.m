function v = PT_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 105);
  end
  v = vInitialized;
end
