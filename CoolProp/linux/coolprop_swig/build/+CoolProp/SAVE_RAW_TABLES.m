function v = SAVE_RAW_TABLES()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 135);
  end
  v = vInitialized;
end
