function v = HmassP_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 116);
  end
  v = vInitialized;
end
