function v = iPhase()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 74);
  end
  v = vInitialized;
end
