function v = MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 168);
  end
  v = vInitialized;
end
