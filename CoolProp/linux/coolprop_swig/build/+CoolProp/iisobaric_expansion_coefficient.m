function v = iisobaric_expansion_coefficient()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 50);
  end
  v = vInitialized;
end
