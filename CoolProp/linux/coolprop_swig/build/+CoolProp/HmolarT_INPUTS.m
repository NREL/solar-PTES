function v = HmolarT_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 108);
  end
  v = vInitialized;
end
