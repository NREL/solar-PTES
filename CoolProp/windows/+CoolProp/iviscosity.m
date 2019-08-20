function v = iviscosity()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 46);
  end
  v = vInitialized;
end
