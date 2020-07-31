function v = HmolarQ_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 101);
  end
  v = vInitialized;
end
