function v = HmolarP_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 117);
  end
  v = vInitialized;
end
