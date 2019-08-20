function v = HmassT_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 109);
  end
  v = vInitialized;
end
