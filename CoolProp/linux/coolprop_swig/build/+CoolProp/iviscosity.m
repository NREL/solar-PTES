function v = iviscosity()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 44);
  end
  v = vInitialized;
end
