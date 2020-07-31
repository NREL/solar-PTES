function v = PHASE_ENVELOPE_STARTING_PRESSURE_PA()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 144);
  end
  v = vInitialized;
end
