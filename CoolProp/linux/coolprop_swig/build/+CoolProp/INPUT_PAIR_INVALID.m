function v = INPUT_PAIR_INVALID()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 96);
  end
  v = vInitialized;
end
