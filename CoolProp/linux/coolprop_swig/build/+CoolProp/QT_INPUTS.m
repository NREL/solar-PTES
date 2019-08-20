function v = QT_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 97);
  end
  v = vInitialized;
end
