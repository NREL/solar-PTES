function v = IFRAC_UNDEFINED()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 90);
  end
  v = vInitialized;
end
