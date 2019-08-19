function v = DmassP_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 114);
  end
  v = vInitialized;
end
