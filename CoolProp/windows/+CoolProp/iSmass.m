function v = iSmass()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 39);
  end
  v = vInitialized;
end
