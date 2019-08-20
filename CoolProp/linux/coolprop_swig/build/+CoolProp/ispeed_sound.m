function v = ispeed_sound()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 48);
  end
  v = vInitialized;
end
