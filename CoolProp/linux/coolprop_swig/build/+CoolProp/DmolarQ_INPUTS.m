function v = DmolarQ_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 103);
  end
  v = vInitialized;
end
