function v = DmolarP_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 115);
  end
  v = vInitialized;
end
