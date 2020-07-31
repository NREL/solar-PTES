function v = SmassT_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 111);
  end
  v = vInitialized;
end
