function v = IFRAC_VOLUME()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 89);
  end
  v = vInitialized;
end
