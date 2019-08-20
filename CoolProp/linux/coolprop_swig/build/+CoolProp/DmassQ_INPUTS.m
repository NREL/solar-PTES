function v = DmassQ_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 104);
  end
  v = vInitialized;
end
