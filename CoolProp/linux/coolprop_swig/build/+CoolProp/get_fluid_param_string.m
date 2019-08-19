function varargout = get_fluid_param_string(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(324,varargin{:});
end
